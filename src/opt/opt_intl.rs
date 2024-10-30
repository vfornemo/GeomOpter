use std::iter::Cycle;

use log::{debug, info, trace};
use serde::de;

use crate::geom::coord::{format_coord, CoordVec};
use crate::matrix::mat_blas_lapack::{mat_dgemm, mat_dsyev};
use crate::matrix::MatFull;
use crate::mol::Molecule;
use crate::data::MAX_STEP_SIZE;
use crate::utils::{formated_output_vec, rms};

#[derive(Clone, Debug)]
pub struct OptIntl{
    pub mol: Molecule,
    pub intl_coords: MatFull<f64>,
    pub inv_hess: MatFull<f64>,
    pub inv_G: MatFull<f64>,
    pub grad_q: MatFull<f64>,
    pub p_qk: MatFull<f64>,
    pub q_k1: MatFull<f64>,


}


impl OptIntl {
    pub fn from(mol: Molecule) -> Self {
        OptIntl {
            mol,
            intl_coords: MatFull::new(),
            inv_hess: MatFull::new(),
            inv_G: MatFull::new(),
            grad_q: MatFull::new(),
            p_qk: MatFull::new(),
            q_k1: MatFull::new(),
        }
    }

    pub fn build(&mut self) {
        self.intl_coords = self.mol.geom.to_intl_coord();
        self.get_inv_G();
        self.get_grad_q();
        self.init_inv_hess();


    }

    pub fn kernel(&mut self) {
        self.geom_opt();
    }


    //[nx, nq]T * [nx, nq] = [nq, nq]
    pub fn get_inv_G(&mut self){
        let G_mat = mat_dgemm(&self.mol.B_mat, &self.mol.B_mat, 't', 'n', 1.0, 0.0);
        
        debug!("G matrix, the product of B with its own transpose (square, dimension {})", self.inv_G.size[0]);
        G_mat.formated_output_f64('t', 5, "debug");
        
        let (evec, mut eval, _) = mat_dsyev(&G_mat, 'v');
        let mut evec = evec.unwrap();

        // eval = eval.into_iter().filter(|x| x.abs() > 1.0e-10).collect();
        // eval.iter_mut().for_each(|x| if x.abs() < 1.0e-10 {*x = 0.0});
        debug!("Eigenvalues of G:");
        formated_output_vec(&eval, 5, "debug");

        // eval = eval.into_iter().rev().collect();
        // evec.data = evec.data.into_iter().rev().collect();
        // let mut zero_vec = vec![0.0;self.mol.geom.nintl-3*self.mol.geom.natm+6];
        // eval.append(&mut zero_vec);

        let eval_mat = MatFull::to_diagonal(eval);
        let mut G_inv = mat_dgemm(&evec, &eval_mat, 'n', 'n', 1.0, 0.0);
        G_inv = mat_dgemm(&G_inv, &evec, 'n', 't', 1.0, 0.0);
        debug!("Inverse G Matrix (square matrix with dimension {}) at the initial structure:", G_inv.size[0]);
        G_inv.formated_output_f64('t', 5, "debug");

        self.inv_G = G_inv;

    }

    pub fn get_grad_q(&mut self) {
        let mut grad_q = mat_dgemm(&self.inv_G, &self.mol.B_mat, 'n', 't', 1.0, 0.0); //[nq, nq] [nx, nq]T -> [nq, nx]
        let mut grad_x = self.mol.grad_tot.clone();
        grad_x.reshape([self.mol.geom.natm*3, 1]); //[nx, 1]
        grad_q = mat_dgemm(&grad_q, &grad_x, 'n', 'n', 1.0, 0.0); //[nq, nx] [nx, 1] -> [nq, 1]
        self.grad_q = grad_q;

        debug!("Initial gradient in terms of the internal coordinates (kcal/mol/Angstrom or kcal/mol/radian):");
        self.grad_q.formated_output_f64('t', 6, "debug");
    }

    pub fn init_inv_hess(&mut self) {
        let mut vec_b = vec![1.0/600.0; self.mol.geom.bond_num];
        let mut vec_a = vec![1.0/150.0; self.mol.geom.angle_num];
        let mut vec_d = vec![1.0/80.0; self.mol.geom.dihedral_num];
        vec_b.append(&mut vec_a);
        vec_b.append(&mut vec_d);
        self.inv_hess = MatFull::to_diagonal(vec_b);
        debug!("Initial guess for the inverse Hessian M in internal coordinates:");
        self.inv_hess.formated_output_f64('t', 6, "debug");
    }


    pub fn geom_opt(&mut self) {

        // let mut cycle = 0;

        // loop {
        //     cycle += 1;
        //     debug!("***** Geometry optimization cycle number   {} ******\n", cycle);
        //     *self = self.optimizer();
        //     self.update_hessian();
        //     info!("Cycle number {}, GRMS = {:-16.7}", cycle, self.grms);
        //     if self.grms <= self.mol.geom.input.rms_tol || cycle > self.mol.geom.input.max_cycle {
        //         break
        //     }
            
        // }

        // if self.grms <= self.mol.geom.input.rms_tol {
        //     info!("##########################\n# Optimization converged #\n##########################");
        //     info!("Final energy at mimimum:{:-16.8} kcal/mol", self.mol.geom.e_tot);
        //     info!("Final coordinates:");
        //     self.mol.geom.formated_output_coord("info");
            
        // } else {
        //     warn!("Optimization did not converge after {} cycles", cycle);
        //     warn!("Please check the input and try again");
        // }
        

    }

    /// Give the predicted new structure
    fn optimizer(&self) -> OptIntl {
        let mut opt = self.clone();
        opt.p_qk = mat_dgemm(&self.inv_hess, &self.grad_q, 'n', 'n', -1.0, 0.0); //[nq,1]
        debug!("Predicted update step in internal coordinates s_k (prior to possible scaling):");
        opt.p_qk.formated_output_f64('t', 6, "debug");
        let mut lambda = 1.0;
        let prms = rms(&opt.p_qk.data);
        if prms > MAX_STEP_SIZE {
            debug!("Predicted step is too long: RMS length:{:-10.6}", prms);
            lambda = MAX_STEP_SIZE/prms;
            // loop {
            //     lambda *= MAX_STEP_SIZE/prms;
            //     prms *= lambda;
            //     if prms <= MAX_STEP_SIZE || cycle > self.mol.geom.input.max_cycle {
            //         break
            //     }
            //     debug!("Predicted step is too long: RMS length:{:-10.6}", prms);
            //     cycle += 1;
            // }
            opt.p_qk *= lambda;
            debug!("Scaled update step in internal coordinates s_k:");
            opt.p_qk.formated_output_f64('t', 6, "debug");
        }

        opt.q_k1 = &self.intl_coords + &opt.p_qk;


        // let p_k = format_coord(opt.p_k.data.clone());
        // let new_coord = self.mol.geom.coord.iter().zip(p_k.iter()).map(|(x, y)| x + y*alpha).collect::<Vec<CoordVec>>();
        // opt.mol.geom.coord = new_coord;
        // opt.mol.geom.update();
        // opt.mol.update();
        // opt.y_k = opt.mol.grad_tot.clone() - self.mol.grad_tot.clone();
        // opt.y_k.reshape([self.mol.geom.natm*3, 1]);
        // opt.grms = rms(&opt.mol.grad_tot.data);


        // debug!("\nNew structure r_k+1 = r_k + s_k = r_k + alpha_k p_k");
        // opt.mol.geom.formated_output_coord("debug");
        // debug!("\nEnergy before update and after:{:-16.8} {:-16.8} And GRMS: {:-16.7}", self.mol.geom.e_tot, opt.mol.geom.e_tot, opt.grms);
        // debug!("\nNew gradient g_k+1:"); 
        // opt.mol.grad_tot.formated_output_f64('t', 6, "debug");	

        opt

    }

    fn opt_car_search(&self) -> Vec<CoordVec> {
        let a = vec![];
        let mut x_k = self.mol.geom.to_car_coord();
        let mut s_qk = self.p_qk.clone();
        let mut cycle = 0;
        loop {
            let mut BGs = mat_dgemm(&self.mol.B_mat, &self.inv_G, 'n', 'n', 1.0, 0.0);
            BGs = mat_dgemm(&BGs, &s_qk, 'n', 'n', 1.0, 0.0);
            x_k += &BGs;

            let step_max = BGs.data.iter().map(|x| x.abs()).max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
            if step_max <= self.mol.geom.input.car_conver || cycle > self.mol.geom.input.max_cycle {
                break
            }

            let x_k = format_coord(x_k.data.clone());
            let mut new_geom = self.mol.geom.clone();
            new_geom.coord = x_k;
            let qkj = new_geom.calc_intl_coord();
            s_qk = &self.q_k1 - &qkj;
            cycle += 1;

        }



        a

    }

    fn update_hessian(&mut self) {

        // let mut s_k = self.p_k.clone() * self.alpha;  
        // s_k.reshape([self.mol.geom.natm*3, 1]);                // s_k [cart, 1]
        // // let v_k = self.inv_hess.mat_mul(&self.y_k, 'n', 'n');    // v_k [cart, 1]
        // let v_k = mat_dgemm(&self.inv_hess, &self.y_k, 'n', 'n', 1.0, 0.0);
        
        // let sy = self.y_k.data.iter().zip(s_k.data.iter()).map(|(x,y)| x*y).sum::<f64>();
        // let vy = self.y_k.data.iter().zip(v_k.data.iter()).map(|(x,y)| x*y).sum::<f64>();
        // // let mut ss = s_k.mat_mul(&s_k, 'n', 't');
        // let mut ss = mat_dgemm(&s_k, &s_k, 'n', 't', 1.0, 0.0);
        // // let mut vs = v_k.mat_mul(&s_k, 'n', 't');
        // let mut vs = mat_dgemm(&v_k, &s_k, 'n', 't', 1.0, 0.0);
        // // let sv = s_k.mat_mul(&v_k, 'n', 't');
        // let sv = mat_dgemm(&s_k, &v_k, 'n', 't', 1.0, 0.0);
        // ss *= sy+vy;
        // ss /= sy*sy;
        // vs += &sv;
        // vs /= sy;

        // self.inv_hess += &ss;
        // self.inv_hess -= &vs;

        // trace!("Updated guess for the Inverse Hessian (equation 12)");
        // self.inv_hess.formated_output_f64('t', 6, "debug");

    }


    pub fn logger(&self) {
        info!("Wilson B Matrix at the initial structure: {} rows of {} elements", self.mol.B_mat.size[1], self.mol.B_mat.size[0]);
        self.mol.B_mat.formated_output_f64('t', 5, "debug");

        // info!("Geometry optimization gradient RMS threshold:{:-16.8}", self.mol.geom.input.rms_tol);

        // info!("\nGradient of the energy for initial structure (in kcal/mol/Angstrom)");
        // self.mol.grad_tot.formated_output_f64('t', 6, "info");

        // debug!("Initial Guess for the Inverse Hessian");
        // self.inv_hess.formated_output_f64('t', 6, "debug");

    }

}