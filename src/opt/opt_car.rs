//! Cartesian coordinate optimization module

#![allow(non_snake_case)]

use crate::geom::coord::{format_coord, CoordVec};
use crate::utils::rms;
use crate::{matrix::MatFull, grad::Gradient};
use crate::utils::constant::WOLFE_C1;
use log::{debug, info, trace, warn};
use crate::matrix::mat_blas_lapack::mat_dgemm;
use crate::io::Result;

/// Cartesian optimization structure
/// # Fields:
/// * `grad`: Gradient
/// * `inv_hess`: Inverse Hessian
/// * `p_k`: p_k
/// * `alpha`: alpha
/// * `y_k`: y_k
/// * `grms`: gradient RMS
#[derive(Clone, Debug)]
pub struct OptCar {
    pub grad: Gradient,
    pub inv_hess: MatFull<f64>,
    pub p_k: MatFull<f64>,
    pub alpha: f64,
    pub y_k: MatFull<f64>,
    pub grms: f64,

}

impl OptCar {
    pub fn from(grad: Gradient) -> Self {
        OptCar {
            grad,
            inv_hess: MatFull::new(),
            p_k: MatFull::new(),  // p_k
            alpha: 0.0,
            y_k: MatFull::new(), // y_k
            grms: 0.0,
        }
    }

    pub fn build(&mut self) {
        self.inv_hess = MatFull::eye(self.grad.geom.natm*3, 1.0/300.0);

    }

    pub fn kernel(&mut self, res: &mut Result) {
        self.geom_opt(res);
    }

    pub fn geom_opt(&mut self, res: &mut Result) {

        let mut cycle = 0;

        loop {
            cycle += 1;
            debug!("***** Geometry optimization cycle number   {} ******\n", cycle);
            *self = self.optimizer();
            self.update_hessian();
            info!("Cycle number {}, GRMS = {:-16.7}", cycle, self.grms);
            if self.grms <= self.grad.geom.input.rms_tol || cycle >= self.grad.geom.input.max_cycle {
                break
            }
            
        }

        if self.grms <= self.grad.geom.input.rms_tol {
            info!("\n##########################\n# Optimization converged #\n##########################");
            info!("Final energy at mimimum:{:-16.8} kcal/mol", self.grad.geom.e_tot);
            info!("Final coordinates:");
            self.grad.geom.formated_output_coord("info");
            res.is_converged = true;
            res.e_tot = Some(self.grad.geom.e_tot);
            res.grms = Some(self.grms);
            res.geom = Some(self.grad.geom.clone());
            res.conv_cyc = Some(cycle);
        } else {
            warn!("Warning: Optimization did not converge after {} cycles", cycle);
            warn!("Please check the input and try again");
            res.is_converged = false;
            res.grms = Some(self.grms);
        }
        

    }

    /// Give the predicted new structure
    fn optimizer(&self) -> OptCar {
        let mut opt = self.clone();
        let mut grad_V = self.grad.grad_tot.clone();
        grad_V.reshape([self.grad.geom.natm*3, 1]);
        // opt.p_k = -self.inv_hess.mat_mul(&grad_V, 'n', 'n');
        opt.p_k = mat_dgemm(&self.inv_hess, &grad_V, 'n', 'n', -1.0, 0.0);
        opt.p_k.reshape([3,self.grad.geom.natm]);

        trace!("Predicted structure change p_k = -M_k grad V(r_k)");
        opt.p_k.formated_output_f64('t', 6, "debug");

        let mut alpha = 0.8;
        let mut cycle = 0;
        while !opt.line_search(alpha) && cycle < self.grad.geom.input.max_cycle {
            alpha *= 0.8;
            cycle += 1;
        }

        opt.alpha = alpha;

        let p_k = format_coord(opt.p_k.data.clone());
        let new_coord = self.grad.geom.coord.iter().zip(p_k.iter()).map(|(x, y)| x + y*alpha).collect::<Vec<CoordVec>>();
        opt.grad.geom.coord = new_coord;
        opt.grad.geom.update();
        opt.grad.update();
        opt.y_k = opt.grad.grad_tot.clone() - self.grad.grad_tot.clone();
        opt.y_k.reshape([self.grad.geom.natm*3, 1]);
        opt.grms = rms(&opt.grad.grad_tot.data);


        debug!("\nNew structure r_k+1 = r_k + s_k = r_k + alpha_k p_k");
        opt.grad.geom.formated_output_coord("debug");
        debug!("\nEnergy before update and after:{:-16.8} {:-16.8} And GRMS: {:-16.7}", self.grad.geom.e_tot, opt.grad.geom.e_tot, opt.grms);
        debug!("\nNew gradient g_k+1:"); 
        opt.grad.grad_tot.formated_output_f64('t', 6, "debug");	

        opt

    }

    fn update_hessian(&mut self) {

        let mut s_k = self.p_k.clone() * self.alpha;  
        s_k.reshape([self.grad.geom.natm*3, 1]);                // s_k [cart, 1]
        let v_k = mat_dgemm(&self.inv_hess, &self.y_k, 'n', 'n', 1.0, 0.0);
        
        let sy = self.y_k.data.iter().zip(s_k.data.iter()).map(|(x,y)| x*y).sum::<f64>();
        let vy = self.y_k.data.iter().zip(v_k.data.iter()).map(|(x,y)| x*y).sum::<f64>();
        let mut ss = mat_dgemm(&s_k, &s_k, 'n', 't', 1.0, 0.0);
        let mut vs = mat_dgemm(&v_k, &s_k, 'n', 't', 1.0, 0.0);
        let sv = mat_dgemm(&s_k, &v_k, 'n', 't', 1.0, 0.0);
        ss *= sy+vy;
        ss /= sy*sy;
        vs += &sv;
        vs /= sy;

        self.inv_hess += &ss;
        self.inv_hess -= &vs;

        trace!("Updated guess for the Inverse Hessian (equation 12)");
        self.inv_hess.formated_output_f64('t', 6, "trace");

    }

    fn line_search(&self, alpha: f64) -> bool {
        let e_old = self.grad.geom.e_tot;
        let p_k = format_coord(self.p_k.data.clone());
        let new_coord = self.grad.geom.coord.iter().zip(p_k.iter()).map(|(x, y)| x + y*alpha).collect::<Vec<CoordVec>>();
        let mut new_geom = self.grad.geom.clone();
        new_geom.coord = new_coord;
        new_geom.update();
        let lhs = new_geom.e_tot;
        let mut rhs = self.p_k.data.iter().zip(self.grad.grad_tot.data.iter()).map(|(x,y)| x*y*WOLFE_C1*alpha).sum::<f64>();
        rhs += e_old;

        trace!("Line search: alpha and energy: {:-11.6} {:-11.6}", alpha, lhs);

        if lhs <= rhs {
            return true
        } else {
            return false
        }
        
    }


    pub fn logger(&self) {
        info!("Geometry optimization gradient RMS threshold:{:-16.8}", self.grad.geom.input.rms_tol);

        info!("\nGradient of the energy for initial structure (in kcal/mol/Angstrom)");
        self.grad.grad_tot.formated_output_f64('t', 6, "info");

        debug!("Initial Guess for the Inverse Hessian");
        self.inv_hess.formated_output_f64('t', 6, "debug");

    }

}
