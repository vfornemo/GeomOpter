//! Bond structure and methods

use crate::utils::constant;
use crate::geom::coord::CoordVec;

#[derive(Debug, Clone)]
/// Enumeration of bond types
///  `CH`, `CC`, `HH`, `XX` \
/// `X = H, C`
pub enum BondType {
    CH,
    CC,
    HH,
    XX,
}

#[derive(Debug, Clone)]
/// Structure of a bond
pub struct Bond {
    /// atom indices of the bonding atoms
    pub atms: [usize;2],
    /// bond type
    pub bond_type: BondType,
    /// bond length
    pub bond_length: f64,
    /// stretching energy or vdw energy, depending on `is_LJ_pair`
    pub e_str: f64, 
    /// whether the bond is counted in vdw energy calculation
    pub is_LJ_pair: bool,
}


impl Bond {

    /// Create a new bond from two atom indices
    pub fn from(atm1: usize, atm2: usize) -> Bond {
        Bond {
            atms: [atm1, atm2],
            bond_type: BondType::XX,
            bond_length: 0.0f64,
            e_str: 0.0f64,
            is_LJ_pair: false,
        }
    }

    /// Get bond type of a bond according to the bonding atoms
    pub fn get_bond_type(&mut self, atoms_idx: &Vec<usize>) {
        let mut c_num = 0;
        for i in 0..2 {
            if atoms_idx[self.atms[i]-1] == 6 {
                c_num += 1;
            }
        }

        if c_num == 2 {
            self.bond_type = BondType::CC;
        } else if c_num == 1 {
            self.bond_type = BondType::CH;
        } else {
            self.bond_type = BondType::HH;
        }
    }

    /// Calculate bond length of a bond
    pub fn get_bond_length(&mut self, coord: &Vec<CoordVec>) {
        self.bond_length = coord[self.atms[0]-1].distance(&coord[self.atms[1]-1]);
    }

    /// Calculate stretching energy of a bond
    pub fn get_e_str(&mut self) {
        self.e_str = match self.bond_type {
            BondType::CC => {constant::Kb_CC*(self.bond_length - constant::R0_CC).powi(2)},
            BondType::CH => {constant::Kb_CH*(self.bond_length - constant::R0_CH).powi(2)},
            BondType::HH => 0.0,
            BondType::XX => 0.0,
        }

    }

    /// Calculate vdw energy of a bond
    pub fn get_e_vdw(&mut self) {
        let epsilon_CH = (constant::EPSILON_C*constant::EPSILON_H).sqrt();
        let sigma_CH = 2.0*(constant::SIGMA_C*constant::SIGMA_H).sqrt();
        let sigma_HH = 2.0*constant::SIGMA_H;
        let sigma_CC = 2.0*constant::SIGMA_C;
        if self.is_LJ_pair {
            self.e_str = match self.bond_type {
                BondType::CH => {4.0*epsilon_CH*((sigma_CH/self.bond_length).powi(12) - (sigma_CH/self.bond_length).powi(6))},
                BondType::HH => {4.0*constant::EPSILON_H*((sigma_HH/self.bond_length).powi(12) - (sigma_HH/self.bond_length).powi(6))},
                BondType::CC => {4.0*constant::EPSILON_C*((sigma_CC/self.bond_length).powi(12) - (sigma_CC/self.bond_length).powi(6))},
                BondType::XX => 0.0,
            }
        }

    }

    
}
