
#![allow(non_snake_case)]

use crate::{data, geom::{Angle, AngleType, BondType, Geom}, matrix::MatFull};
use data::{A0_CCC, A0_XCX, Kb_CC, Kb_CH, Ka_CCC, Ka_XCX, N_XCCX, A_XCCX, R0_CC, R0_CH, SIGMA_C, SIGMA_H, EPSILON_C, EPSILON_H};

#[derive(Debug)]
pub struct Molecule {
    pub geom: Geom,
    pub B_mat: MatFull<f64>,
    pub grad_str: MatFull<f64>,
    pub grad_bend: MatFull<f64>,
    pub grad_tors: MatFull<f64>,
    pub grad_vdw: MatFull<f64>,
    pub grad_tol: MatFull<f64>,

    // ...
}

impl Molecule {
    pub fn new() -> Molecule {
        Molecule {
            geom: Geom::new(),
            B_mat: MatFull::new(),
            grad_str: MatFull::new(),
            grad_bend: MatFull::new(),
            grad_tors: MatFull::new(),
            grad_vdw: MatFull::new(),
            grad_tol: MatFull::new(),

        }
    }

    pub fn from(geom: Geom) -> Molecule {
        let natm = geom.atm_num;
        let n_intl = geom.bond_num + geom.angle_num + geom.dihedral_num;
        Molecule {
            geom,
            B_mat: MatFull::<f64>::from_vec([natm*3, n_intl], vec![0.0; 3*natm*n_intl]),
            grad_str: MatFull::<f64>::from_vec([3, natm], vec![0.0; 3*natm]),
            grad_bend: MatFull::<f64>::from_vec([3, natm], vec![0.0; 3*natm]),
            grad_tors: MatFull::<f64>::from_vec([3, natm], vec![0.0; 3*natm]),
            grad_vdw: MatFull::<f64>::from_vec([3, natm], vec![0.0; 3*natm]),
            grad_tol: MatFull::<f64>::from_vec([3, natm], vec![0.0; 3*natm]),
        }
    }

    pub fn build(&mut self) {
        self.get_B_mat();

        self.get_g_str(); 

        self.get_g_bend();

        self.get_g_tors();

        self.get_g_vdw();

        self.get_g_tol();


    }

    pub fn logger(&self) {
        println!("Wilson B Matrix at the initial structure:  {} rows of  {} elements", self.B_mat.size[1], self.B_mat.size[0]);
        self.B_mat.formated_output_f64('t',5);

        println!("Analytical gradient of overall energy: (in kcal/mol/angstrom)");
        self.grad_tol.formated_output_f64('t',6);

        println!("Analytical gradient of stretching energy:");
        self.grad_str.formated_output_f64('t',6);

        println!("Analytical gradient of bending energy:");
        self.grad_bend.formated_output_f64('t',6);

        println!("Analytical gradient of torsional energy:");
        self.grad_tors.formated_output_f64('t',6);
        
        println!("Analytical gradient of VDW energy:");
        self.grad_vdw.formated_output_f64('t',6);




    }

    //               nq
    //     b1/x1 ... an/x1 ... dn/x1
    // nx  b1/y1 .
    //     b1/z1   .
    //     b1/x2     .
    //      ...
    //     b1/zn    ...        dn/zn
    /// Calculate the Wilson B matrix [nx, nq]
    /// nx: number of cartesian coordinates
    /// nq: number of internal coordinates
    fn get_B_mat(&mut self) {
        // partial bond to partial cartesian
        for bi in 0..self.geom.bond_num {
            let bond = &self.geom.bonds[bi];
            let coord_a = self.geom.coord[bond.atms[0]-1];
            let coord_b = self.geom.coord[bond.atms[1]-1];
            let r_AB = coord_a - coord_b;
            r_AB.to_array().iter().zip(self.B_mat.iter_col_range_mut(bi, 3*bond.atms[0]-3, 3*bond.atms[0])).for_each(|(x,y)| *y += x/bond.bond_length);
            r_AB.to_array().iter().zip(self.B_mat.iter_col_range_mut(bi, 3*bond.atms[1]-3, 3*bond.atms[1])).for_each(|(x,y)| *y -= x/bond.bond_length);

        }

        let mut offset = self.geom.bond_num;

        // partial angle to partial cartesian
        for ai in 0..self.geom.angle_num {
            let angle = &self.geom.angles[ai];
            let coord_a = self.geom.coord[angle.atms[0]-1];
            let coord_b = self.geom.coord[angle.atms[1]-1];
            let coord_c = self.geom.coord[angle.atms[2]-1];
            let r_BA = coord_a - coord_b;
            let r_BC = coord_c - coord_b;
            let p = r_BA.cross(&r_BC);
            let grad_a = r_BA.cross(&p)/(r_BA.square()*p.norm());
            let grad_c = -r_BC.cross(&p)/(r_BC.square()*p.norm());
            let grad_b = -grad_a - grad_c;
            grad_a.to_array().iter().zip(self.B_mat.iter_col_range_mut(ai+offset, 3*angle.atms[0]-3, 3*angle.atms[0])).for_each(|(x,y)| *y += x);
            grad_b.to_array().iter().zip(self.B_mat.iter_col_range_mut(ai+offset, 3*angle.atms[1]-3, 3*angle.atms[1])).for_each(|(x,y)| *y += x);
            grad_c.to_array().iter().zip(self.B_mat.iter_col_range_mut(ai+offset, 3*angle.atms[2]-3, 3*angle.atms[2])).for_each(|(x,y)| *y += x);

        }

        offset += self.geom.angle_num;

        // partial dihedral to partial cartesian 
        for di in 0..self.geom.dihedral_num {
            let dihedral = &self.geom.dihedrals[di];
            let coord_a = self.geom.coord[dihedral.atms[0]-1];
            let coord_b = self.geom.coord[dihedral.atms[1]-1];
            let coord_c = self.geom.coord[dihedral.atms[2]-1];
            let coord_d = self.geom.coord[dihedral.atms[3]-1];
            let r_AB = coord_b - coord_a;
            let r_AC = coord_c - coord_a;
            let r_BC = coord_c - coord_b;
            let r_CD = coord_d - coord_c;
            let r_BD = coord_d - coord_b;
            let t = r_AB.cross(&r_BC);
            let u = r_BC.cross(&r_CD);
            let grad_a = t.cross(&r_BC).cross(&r_BC)/(t.square()*r_BC.norm());
            let grad_d = -u.cross(&r_BC).cross(&r_BC)/(u.square()*r_BC.norm());
            let grad_b = r_AC.cross(&t.cross(&r_BC))/(t.square()*r_BC.norm()) - u.cross(&r_BC).cross(&r_CD)/(u.square()*r_BC.norm());
            let grad_c = t.cross(&r_BC).cross(&r_AB)/(t.square()*r_BC.norm()) - r_BD.cross(&u.cross(&r_BC))/(u.square()*r_BC.norm());
            grad_a.to_array().iter().zip(self.B_mat.iter_col_range_mut(di+offset, 3*dihedral.atms[0]-3, 3*dihedral.atms[0])).for_each(|(x,y)| *y += x);
            grad_b.to_array().iter().zip(self.B_mat.iter_col_range_mut(di+offset, 3*dihedral.atms[1]-3, 3*dihedral.atms[1])).for_each(|(x,y)| *y += x);
            grad_c.to_array().iter().zip(self.B_mat.iter_col_range_mut(di+offset, 3*dihedral.atms[2]-3, 3*dihedral.atms[2])).for_each(|(x,y)| *y += x);
            grad_d.to_array().iter().zip(self.B_mat.iter_col_range_mut(di+offset, 3*dihedral.atms[3]-3, 3*dihedral.atms[3])).for_each(|(x,y)| *y += x);
            
        }


    }


    /// Get the gradient of the stretch energy wrt the coordinates
    fn get_g_str(&mut self) {
        for bi in 0..self.geom.bond_num {
            let bond = &self.geom.bonds[bi];
            // rBA/rAB
            let grad_r = self.B_mat.iter_col_range(bi, 3*bond.atms[0]-3, 3*bond.atms[0]); 
            match bond.bond_type {
                BondType::CC => {
                                grad_r.clone().zip(self.grad_str.iter_col_mut(bond.atms[0]-1)).for_each(|(x,y)| *y += 2.0*Kb_CC*(bond.bond_length - R0_CC)*x);
                                grad_r.zip(self.grad_str.iter_col_mut(bond.atms[1]-1)).for_each(|(x,y)| *y -= 2.0*Kb_CC*(bond.bond_length - R0_CC)*x);
                                },
                       
                BondType::CH => {
                                grad_r.clone().zip(self.grad_str.iter_col_mut(bond.atms[0]-1)).for_each(|(x,y)| *y += 2.0*Kb_CH*(bond.bond_length - R0_CH)*x);
                                grad_r.zip(self.grad_str.iter_col_mut(bond.atms[1]-1)).for_each(|(x,y)| *y -= 2.0*Kb_CH*(bond.bond_length - R0_CH)*x);
                                },
                BondType::HH => (),
                BondType::XX => (),
            }
    
        }

    }

    fn get_g_bend(&mut self) {
        for ai in 0..self.geom.angle_num {
            let angle = &self.geom.angles[ai];
            // rBA/rAB
            let grad_a_a = self.B_mat.iter_col_range(ai+self.geom.bond_num, 3*angle.atms[0]-3, 3*angle.atms[0]); 
            let grad_a_b = self.B_mat.iter_col_range(ai+self.geom.bond_num, 3*angle.atms[1]-3, 3*angle.atms[1]);
            let grad_a_c = self.B_mat.iter_col_range(ai+self.geom.bond_num, 3*angle.atms[2]-3, 3*angle.atms[2]);
            match angle.angle_type {
                AngleType::CCC => {
                                grad_a_a.zip(self.grad_bend.iter_col_mut(angle.atms[0]-1)).for_each(|(x,y)| *y += 2.0*Ka_CCC*(angle.angle - A0_CCC).to_radians()*x);
                                grad_a_b.zip(self.grad_bend.iter_col_mut(angle.atms[1]-1)).for_each(|(x,y)| *y += 2.0*Ka_CCC*(angle.angle - A0_CCC).to_radians()*x);
                                grad_a_c.zip(self.grad_bend.iter_col_mut(angle.atms[2]-1)).for_each(|(x,y)| *y += 2.0*Ka_CCC*(angle.angle - A0_CCC).to_radians()*x);
                                
                                },
                       
                AngleType::XCX => {
                                grad_a_a.zip(self.grad_bend.iter_col_mut(angle.atms[0]-1)).for_each(|(x,y)| *y += 2.0*Ka_XCX*(angle.angle - A0_XCX).to_radians()*x);
                                grad_a_b.zip(self.grad_bend.iter_col_mut(angle.atms[1]-1)).for_each(|(x,y)| *y += 2.0*Ka_XCX*(angle.angle - A0_XCX).to_radians()*x);
                                grad_a_c.zip(self.grad_bend.iter_col_mut(angle.atms[2]-1)).for_each(|(x,y)| *y += 2.0*Ka_XCX*(angle.angle - A0_XCX).to_radians()*x);
                                },
                AngleType::XXX => (),
            }
    
        }

    }
    

    fn get_g_tors(&mut self) {
        for di in 0..self.geom.dihedral_num {
            let dihedral = &self.geom.dihedrals[di];
            let offset = self.geom.bond_num + self.geom.angle_num;
            // rBA/rAB
            let grad_d_a = self.B_mat.iter_col_range(di+offset, 3*dihedral.atms[0]-3, 3*dihedral.atms[0]); 
            let grad_d_b = self.B_mat.iter_col_range(di+offset, 3*dihedral.atms[1]-3, 3*dihedral.atms[1]);
            let grad_d_c = self.B_mat.iter_col_range(di+offset, 3*dihedral.atms[2]-3, 3*dihedral.atms[2]);
            let grad_d_d = self.B_mat.iter_col_range(di+offset, 3*dihedral.atms[3]-3, 3*dihedral.atms[3]);
            grad_d_a.zip(self.grad_tors.iter_col_mut(dihedral.atms[0]-1)).for_each(|(x,y)| *y -= N_XCCX*A_XCCX*(N_XCCX*dihedral.dihedral.to_radians()).sin()*x);	
            grad_d_b.zip(self.grad_tors.iter_col_mut(dihedral.atms[1]-1)).for_each(|(x,y)| *y -= N_XCCX*A_XCCX*(N_XCCX*dihedral.dihedral.to_radians()).sin()*x);
            grad_d_c.zip(self.grad_tors.iter_col_mut(dihedral.atms[2]-1)).for_each(|(x,y)| *y -= N_XCCX*A_XCCX*(N_XCCX*dihedral.dihedral.to_radians()).sin()*x);
            grad_d_d.zip(self.grad_tors.iter_col_mut(dihedral.atms[3]-1)).for_each(|(x,y)| *y -= N_XCCX*A_XCCX*(N_XCCX*dihedral.dihedral.to_radians()).sin()*x);
    
        }
        
    }

    
    fn get_g_vdw(&mut self) {
        for ap in self.geom.atm_pair.iter(){
            if ap.is_LJ_pair {
                let vij = self.geom.coord[ap.atms[0]-1] - self.geom.coord[ap.atms[1]-1];
                let rij2 = vij.square();
                match ap.bond_type {
                    BondType::CC => { let sigma_CC = 2.0*SIGMA_C;
                                      let Aij = 4.0*sigma_CC.powi(12);
                                      let Bij = 4.0*sigma_CC.powi(6);
                                      vij.clone().to_array().iter().zip(self.grad_vdw.iter_col_mut(ap.atms[0]-1)).for_each(|(x,y)| *y += EPSILON_C*x*(-12.0*Aij/rij2.powi(7) + 6.0*Bij/rij2.powi(4)));
                                      vij.to_array().iter().zip(self.grad_vdw.iter_col_mut(ap.atms[1]-1)).for_each(|(x,y)| *y += EPSILON_C*x*(-12.0*Aij/rij2.powi(7) + 6.0*Bij/rij2.powi(4)));
                                    },
                    BondType::CH => { let sigma_CH = 2.0*(SIGMA_C*SIGMA_H).sqrt();
                                      let Aij = 4.0*sigma_CH.powi(12);
                                      let Bij = 4.0*sigma_CH.powi(6);
                                      let epsilon_CH = (EPSILON_C*EPSILON_H).sqrt();
                                      vij.clone().to_array().iter().zip(self.grad_vdw.iter_col_mut(ap.atms[0]-1)).for_each(|(x,y)| *y += epsilon_CH*x*(-12.0*Aij/rij2.powi(7) + 6.0*Bij/rij2.powi(4)));
                                      vij.to_array().iter().zip(self.grad_vdw.iter_col_mut(ap.atms[1]-1)).for_each(|(x,y)| *y -= epsilon_CH*x*(-12.0*Aij/rij2.powi(7) + 6.0*Bij/rij2.powi(4)));
                                    },
                    BondType::HH => { let sigma_HH = 2.0*SIGMA_H;
                                      let Aij = 4.0*sigma_HH.powi(12);
                                      let Bij = 4.0*sigma_HH.powi(6);
                                      vij.clone().to_array().iter().zip(self.grad_vdw.iter_col_mut(ap.atms[0]-1)).for_each(|(x,y)| *y += EPSILON_H*x*(-12.0*Aij/rij2.powi(7) + 6.0*Bij/rij2.powi(4)));
                                      vij.to_array().iter().zip(self.grad_vdw.iter_col_mut(ap.atms[1]-1)).for_each(|(x,y)| *y -= EPSILON_H*x*(-12.0*Aij/rij2.powi(7) + 6.0*Bij/rij2.powi(4)));
                                    },
                    BondType::XX => (),
                }


            }
        }
        
        
    }

    fn get_g_tol(&mut self) {
        self.grad_tol = self.grad_str.clone();
        self.grad_tol += &self.grad_bend;
        self.grad_tol += &self.grad_tors;
        self.grad_tol += &self.grad_vdw;

    }

}

