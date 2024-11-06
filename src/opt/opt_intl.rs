use std::iter::Cycle;

use log::{debug, info, trace, warn};
use serde::de;

use crate::geom::coord::{format_coord, CoordVec};
use crate::matrix::mat_blas_lapack::{mat_dgemm, mat_dsyev};
use crate::matrix::MatFull;
use crate::mol::Molecule;
use crate::data::MAX_STEP_SIZE;
use crate::utils::{formated_output_vec, rms, dihedral_diff_check_rad};

#[derive(Clone, Debug)]
pub struct OptIntl{
    pub mol: Molecule,
    pub intl_coords: MatFull<f64>,
    pub inv_hess: MatFull<f64>,
    pub inv_G: MatFull<f64>,
    pub grad_q: MatFull<f64>,
    pub p_qk: MatFull<f64>,
    pub q_k1: MatFull<f64>,
    pub y_k: MatFull<f64>,
    pub s_k: MatFull<f64>,
    pub grms: f64,

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
            y_k: MatFull::new(),
            s_k: MatFull::new(),
            grms: 0.0,
        }
    }

    pub fn build(&mut self) {
        self.intl_coords = self.mol.geom.to_intl_coord();
        self.get_inv_G();
        self.get_grad_q();
        self.grms = rms(&self.mol.grad_tot.data);
        debug!("Initial gradient in terms of the internal coordinates (kcal/mol/Angstrom or kcal/mol/radian):");
        self.grad_q.formated_output_f64('t', 6, "debug");
        self.init_inv_hess();
        debug!("Initial guess for the inverse Hessian M in internal coordinates:");
        self.inv_hess.formated_output_f64('t', 5, "debug");
    }

    pub fn update(&mut self) {
        self.inv_G.reset();
        self.grad_q.reset();
        // self.inv_hess.reset();

        self.get_inv_G1();
        self.get_grad_q();
        self.grms = rms(&self.mol.grad_tot.data);
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
        let evec = evec.unwrap();

        debug!("Eigenvalues of G:");
        formated_output_vec(&eval, 5, "debug");

        eval.iter_mut().for_each(|x| if x.abs() < 1.0e-14 {*x = 0.0} else {*x = 1.0/ *x });

        let eval_mat = MatFull::to_diagonal(eval);
        let mut G_inv = mat_dgemm(&evec, &eval_mat, 'n', 'n', 1.0, 0.0);
        G_inv = mat_dgemm(&G_inv, &evec, 'n', 't', 1.0, 0.0);
        debug!("Inverse G Matrix (square matrix with dimension {}) at the initial structure:", G_inv.size[0]);
        G_inv.formated_output_f64('t', 5, "debug");

        self.inv_G = G_inv;

    }

    //[nx, nq]T * [nx, nq] = [nq, nq]
    /// equivalent to get_inv_G, but more brief
    pub fn get_inv_G1(&mut self){
        let G_mat = mat_dgemm(&self.mol.B_mat, &self.mol.B_mat, 't', 'n', 1.0, 0.0);
        let (evec, mut eval, _) = mat_dsyev(&G_mat, 'v');
        let evec = evec.unwrap();
        eval.iter_mut().for_each(|x| if x.abs() < 1.0e-14 {*x = 0.0} else {*x = 1.0/ *x });
        let eval_mat = MatFull::to_diagonal(eval);
        let mut G_inv = mat_dgemm(&evec, &eval_mat, 'n', 'n', 1.0, 0.0);
        G_inv = mat_dgemm(&G_inv, &evec, 'n', 't', 1.0, 0.0);
        self.inv_G = G_inv;
    }

    pub fn get_grad_q(&mut self) {
        let mut grad_q = mat_dgemm(&self.inv_G, &self.mol.B_mat, 'n', 't', 1.0, 0.0); //[nq, nq] [nx, nq]T -> [nq, nx]
        let mut grad_x = self.mol.grad_tot.clone();
        grad_x.reshape([self.mol.geom.natm*3, 1]); //[nx, 1]
        grad_q = mat_dgemm(&grad_q, &grad_x, 'n', 'n', 1.0, 0.0); //[nq, nx] [nx, 1] -> [nq, 1]
        self.grad_q = grad_q;
    }


    pub fn init_inv_hess(&mut self) {
        let mut vec_b = vec![1.0/600.0; self.mol.geom.nbond];
        let mut vec_a = vec![1.0/150.0; self.mol.geom.nangle];
        let mut vec_d = vec![1.0/80.0; self.mol.geom.ndihedral];
        vec_b.append(&mut vec_a);
        vec_b.append(&mut vec_d);
        self.inv_hess = MatFull::to_diagonal(vec_b);

    }


    pub fn geom_opt(&mut self) {

        info!("##################################\n# Start of Geometry Optimization #\n##################################\n");
        info!("Reminder: maximum step size is {:-8.3}", MAX_STEP_SIZE);

        let mut cycle = 0;

        loop {
            cycle += 1;
            debug!("\n*** Optimization Cycle  {} ***\n", cycle);
            let e_old = self.mol.geom.e_tot;
            *self = self.optimizer();
            
            debug!("New coordinates:");
            self.mol.geom.formated_output_coord("debug");
            debug!("New set of internals q (note these are the ones that correspond to the best fit Cartesians):");
            self.intl_coords.formated_output_f64('t', 4, "debug");
            debug!("Wilson B Matrix at the new structure:");
            self.mol.B_mat.formated_output_f64('t', 5, "debug");
            debug!("Inverse G Matrix at the new structure:");
            self.inv_G.formated_output_f64('t', 5, "debug");
            debug!("Gradient in terms of the internal coordinates:");
            self.grad_q.formated_output_f64('t', 6, "debug");
            debug!("Old and new energies:{:-10.5}{:-10.5} And GRMS:{:-8.4}", e_old, self.mol.geom.e_tot, self.grms);
            info!("Cycle number {}, new energy = {:-10.6}, GRMS = {:-8.4}", cycle, self.mol.geom.e_tot, self.grms);

            if self.grms <= self.mol.geom.input.rms_tol || cycle >= self.mol.geom.input.max_cycle {
                break
            }

            self.update_hessian();
            
        }

        if self.grms <= self.mol.geom.input.rms_tol {
            info!("##########################\n# Optimization converged #\n##########################");
            info!("Final energy at mimimum:{:-16.8} kcal/mol", self.mol.geom.e_tot);
            info!("Final coordinates:");
            self.mol.geom.formated_output_coord("info");
            
        } else {
            warn!("Optimization did not converge after {} cycles", cycle);
            warn!("Please check the input and try again");
        }
        

    }

    /// Give the predicted new structure
    fn optimizer(&self) -> OptIntl {
        let mut opt = self.clone();
        opt.p_qk = mat_dgemm(&self.inv_hess, &self.grad_q, 'n', 'n', -1.0, 0.0); //[nq,1]
        debug!("Predicted update step in internal coordinates s_k (prior to possible scaling):");
        opt.p_qk.formated_output_f64('t', 6, "debug");
        let prms = rms(&opt.p_qk.data);
        if prms > MAX_STEP_SIZE {
            debug!("Predicted step is too long: RMS length:{:-10.6}", prms);
            let lambda = MAX_STEP_SIZE/prms;
            opt.p_qk *= lambda;
            debug!("Scaled update step in internal coordinates s_k:");
            opt.p_qk.formated_output_f64('t', 6, "debug");
        }

        opt.q_k1 = &self.intl_coords + &opt.p_qk; // q_k+1 = q_k + s_qk (p_qk)

        let new_coord = opt.opt_car_search();
        opt.mol.geom.coord = new_coord;
        opt.mol.geom.update();
        opt.mol.update();
        opt.intl_coords = opt.mol.geom.calc_intl_coord();
        opt.update();
        opt.s_k = &opt.intl_coords - &self.intl_coords; // s_k = q_k+1 - q_k
        opt.y_k = &opt.grad_q - &self.grad_q; // y_k = grad_q_k+1 - grad_q_k

        opt

    }

    fn opt_car_search(&mut self) -> Vec<CoordVec> {

        debug!("    %% Iterative determination of optimal Cartesian coordinates:  %%");
        debug!("    tolerated diff for Cartsian coordinates:{:-16.8}", self.mol.geom.input.car_conver);

        let mut s_qk = self.p_qk.clone();    // s_qk^(0) = p_qk
        let mut x_k = self.mol.geom.to_car_coord();  // x_k+1^(0) = x_k

        let BG = mat_dgemm(&self.mol.B_mat, &self.inv_G, 'n', 'n', 1.0, 0.0); 
        let mut BGs = mat_dgemm(&BG, &s_qk, 'n', 'n', 1.0, 0.0); //B^T G^- s_qk^(j) 
        x_k += &BGs; // x_k+1^(j+1) = x_k+1^(j) + BGs

        debug!("Initially predicted dx = BTG-sk:");
        BGs.formated_output_f64('t', 4, "debug");

        let mut cycle = 1;

        loop {

            debug!("Cartesian fitting iteration number {}", cycle);

            let new_coord = format_coord(x_k.data.clone());  // x_(k+1)^(j+1)
            let mut new_geom = self.mol.geom.clone();
            new_geom.coord = new_coord;
            let qkj = new_geom.calc_intl_coord();  // q_(k+1)^(j+1)
            
            debug!("current set of internals q_(k+1)^(j+1):");
            qkj.formated_output_f64('t', 4, "debug");

            s_qk = &self.q_k1 - &qkj;  // s_qk^(j+1) = q_k+1 - q_(k+1)^(j+1)

            // check dihedral angle difference. If they are out of range, correct them
            s_qk.data[self.mol.geom.nbond + self.mol.geom.nangle..].iter_mut().for_each(|x| dihedral_diff_check_rad(x)); 

            debug!("difference between these internals and the desired internals, s_q,k^j+1:");
            s_qk.formated_output_f64('t', 4, "debug");

            BGs = mat_dgemm(&BG, &s_qk, 'n', 'n', 1.0, 0.0);
            x_k = new_geom.to_car_coord();  // x_k+1^(j+1)
            x_k += &BGs;

            debug!("Corresponding Cartesians x+k+1^j:");
            x_k.formated_output_f64('t', 4, "debug");

            let step_max = BGs.data.iter().map(|x| x.abs()).max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
            debug!("Maximum change in x from previous iteration:{:-10.6}", step_max);
            
            if step_max.abs() <= self.mol.geom.input.car_conver || cycle >= self.mol.geom.input.max_cycle {
                debug!("fitting of Cartesian coordinates considered converged.");
                break
            }
            cycle += 1;
        }

        self.p_qk = s_qk.clone(); // here we need to confirm
        format_coord(x_k.data)

    }

    fn update_hessian(&mut self) {

        // debug!("y_k is:");
        // self.y_k.formated_output_f64('t', 6, "debug");

        // debug!("s_k is:");
        // self.s_k.formated_output_f64('t', 6, "debug");

        let v_k = mat_dgemm(&self.inv_hess, &self.y_k, 'n', 'n', 1.0, 0.0);
        
        // debug!("v_k is:");
        // v_k.formated_output_f64('t', 6, "debug");

        let sy = self.y_k.data.iter().zip(self.s_k.data.iter()).map(|(x,y)| x*y).sum::<f64>();
        let vy = self.y_k.data.iter().zip(v_k.data.iter()).map(|(x,y)| x*y).sum::<f64>();
        let mut ss = mat_dgemm(&self.s_k, &self.s_k, 'n', 't', 1.0, 0.0);
        let mut vs = mat_dgemm(&v_k, &self.s_k, 'n', 't', 1.0, 0.0);
        let sv = mat_dgemm(&self.s_k, &v_k, 'n', 't', 1.0, 0.0);
        ss *= sy+vy;
        ss /= sy*sy;
        vs += &sv;
        vs /= sy;

        self.inv_hess += &ss;
        self.inv_hess -= &vs;

        debug!("Updated approximate inverse Hessian M in internal coordinates:");
        self.inv_hess.formated_output_f64('t', 6, "debug");

    }


    pub fn logger(&self) {

        info!("Geometry optimization gradient RMS threshold:{:-16.8}", self.mol.geom.input.rms_tol);

        info!("Initial potential energy:");
        info!("{:-10.6}", self.mol.geom.e_tot);

        info!("Wilson B Matrix at the initial structure: {} rows of {} elements", self.mol.B_mat.size[1], self.mol.B_mat.size[0]);
        self.mol.B_mat.formated_output_f64('t', 5, "debug");

        // info!("\nGradient of the energy for initial structure (in kcal/mol/Angstrom)");
        // self.mol.grad_tot.formated_output_f64('t', 6, "info");

        // debug!("Initial Guess for the Inverse Hessian");
        // self.inv_hess.formated_output_f64('t', 6, "debug");

    }

}