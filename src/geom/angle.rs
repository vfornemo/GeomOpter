//! Angle structure and methods

use crate::geom::coord;
use crate::geom::coord::CoordVec;
use crate::utils::constant;


/// Enumeration of angle types
/// `CCC`, `XCX`, `XXX` \
/// `X = H, C`
#[derive(Debug, Clone)]
pub enum AngleType {
    CCC,
    XCX,
    XXX,
}

/// Structure of a bond angle
#[derive(Debug, Clone)]
pub struct Angle {
    /// atom indices of the three bonding atoms
    pub atms: [usize;3],
    /// angle type
    pub angle_type: AngleType,
    /// angle in degrees
    pub angle: f64,
    /// bending energy
    pub e_bend: f64, 
}

impl Angle {

    /// Create a new angle from three atom indices
    pub fn from(atms: [usize;3]) -> Angle {
        Angle {
            atms,
            angle_type: AngleType::XXX,
            angle: 0.0f64,
            e_bend: 0.0f64,
        }
    }

    /// Calculate bond angle
    pub fn get_angle(&mut self, coord: &Vec<CoordVec>) {
        self.angle = coord::get_angle(&coord[self.atms[0]-1], &coord[self.atms[1]-1], &coord[self.atms[2]-1]);
    }

    /// Get angle type according to the bonding atoms
    pub fn get_angle_type(&mut self, atoms_idx: &Vec<usize>) {
        if atoms_idx[self.atms[0]-1] == 6 && atoms_idx[self.atms[1]-1] == 6 && atoms_idx[self.atms[2]-1] == 6 {
            self.angle_type = AngleType::CCC;
        } else {
            self.angle_type = AngleType::XCX;
        }
    }

    /// Calculate bending energy
    pub fn get_e_bend(&mut self) {
        self.e_bend = match self.angle_type {
            AngleType::CCC => {constant::Ka_CCC*(self.angle - constant::A0_CCC).to_radians().powi(2)},
            AngleType::XCX => {constant::Ka_XCX*(self.angle - constant::A0_XCX).to_radians().powi(2)},
            AngleType::XXX => 0.0,
        }
    }

    
}
