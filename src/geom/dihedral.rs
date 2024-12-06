//! Dihedral structure and methods

use crate::utils::constant;
use crate::geom::coord::CoordVec;
use crate::geom::coord;

/// Structure of a dihedral angle
#[derive(Debug, Clone)]
pub struct Dihedral {
    /// atom indices of the four bonding atoms
    pub atms: [usize;4],
    /// dihedral angle in degrees
    pub dihedral: f64,
    /// torsion energy
    pub e_tor: f64, 
}

impl Dihedral {
    /// Create a new dihedral angle from four atom indices
    pub fn from(atms: [usize;4]) -> Dihedral {
        Dihedral {
            atms,
            dihedral: 0.0f64,
            e_tor: 0.0f64,
        }
    }

    /// Calculate dihedral angle
    pub fn get_dihedral(&mut self, coord: &Vec<CoordVec>) {
        self.dihedral = coord::get_dihedral_deg(&coord[self.atms[0]-1], &coord[self.atms[1]-1], &coord[self.atms[2]-1], &coord[self.atms[3]-1]);
    }

    /// Calculate torsion energy
    pub fn get_e_tor(&mut self) {
        self.e_tor = constant::A_XCCX*(1.0 + (constant::N_XCCX*self.dihedral.to_radians()).cos());
    }
    
}
