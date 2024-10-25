use core::f64;

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct CoordVec {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl CoordVec {
    pub fn new() -> CoordVec {
        CoordVec {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    pub fn from(vec: &[f64;3]) -> CoordVec {
        CoordVec {
            x: vec[0],
            y: vec[1],
            z: vec[2],
        }
    }

    pub fn distance(&self, other: &CoordVec) -> f64 {
        ((self.x-other.x).powi(2) + (self.y-other.y).powi(2) + (self.z-other.z).powi(2)).sqrt()
    }

    /// Return the cross product of two vectors
    /// Note: Check if the CoordVec is a vector or a point
    pub fn cross(&self, other: &CoordVec) -> CoordVec {
        CoordVec {
            x: self.y*other.z - self.z*other.y,
            y: self.z*other.x - self.x*other.z,
            z: self.x*other.y - self.y*other.x,
        }
    }

    /// Return the dot product of two vectors
    /// Note: Check if the CoordVec is a vector or a point
    pub fn dot(&self, other: &CoordVec) -> f64 {
        self.x*other.x + self.y*other.y + self.z*other.z
    }

    /// Return the norm of a vector
    pub fn norm(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }
    
    /// Return the square of the norm of a vector v**2
    pub fn square(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }

    /// Return the angle between two vectors in degrees
    pub fn angle_deg(&self, other: &CoordVec) -> f64 {
        let cos_theta = self.dot(other) / (self.norm() * other.norm());
        cos_theta.acos()*180.0/f64::consts::PI
    }

    /// Return the angle between two vectors in radians
    pub fn angle_rad(&self, other: &CoordVec) -> f64 {
        let cos_theta = self.dot(other) / (self.norm() * other.norm());
        cos_theta.acos()
    }

    pub fn to_array(&self) -> [f64;3] {
        [self.x, self.y, self.z]
    }

    
}

impl<'a> Sub for &'a CoordVec
{
    type Output = CoordVec;
    fn sub(self, other: Self) -> CoordVec {
        CoordVec {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Sub for CoordVec
{
    type Output = Self;
    fn sub(self, other: Self) -> Self::Output {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl SubAssign for CoordVec
{
    fn sub_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        };
    }
}


impl<'a> Add for &'a CoordVec
{
    type Output = CoordVec;
    fn add(self, other: Self) -> CoordVec {
         CoordVec{
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl AddAssign for CoordVec
{
    fn add_assign(&mut self, other: Self) {
        *self = Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        };
    }
}

impl MulAssign<f64> for CoordVec
{
    fn mul_assign(&mut self, a: f64) {
        self.x *= a;
        self.y *= a;
        self.z *= a;
    }
}

impl Mul<f64> for CoordVec
{
    type Output = Self;
    fn mul(self, a: f64) -> Self::Output {
        Self {
            x: self.x * a,
            y: self.y * a,
            z: self.z * a,
        }
    }
}

impl DivAssign<f64> for CoordVec
{
    fn div_assign(&mut self, a: f64) {
        self.x /= a;
        self.y /= a;
        self.z /= a;
    }
}

impl Div<f64> for CoordVec
{
    type Output = Self;
    fn div(self, a: f64) -> Self::Output {
        Self {
            x: self.x / a,
            y: self.y / a,
            z: self.z / a,
        }
    }
}

impl Neg for CoordVec
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}


/// Get dihedral angle for atom 1, 2, 3, 4
/// 
///           3 - 4
///         /
///   1 - 2
/// 
// pub fn get_dihedral_v1(atm1: &[f64;3], atm2: &[f64;3], atm3: &[f64;3], atm4: &[f64;3]) -> f64 {
//     let vec1 = [atm2[0]-atm1[0], atm2[1]-atm1[1], atm2[2]-atm1[2]];
//     let vec2 = [atm3[0]-atm2[0], atm3[1]-atm2[1], atm3[2]-atm2[2]];
//     let vec3 = [atm4[0]-atm3[0], atm4[1]-atm3[1], atm4[2]-atm3[2]];
//     let t = utils::cross_product(&vec1, &vec2);
//     let u = utils::cross_product(&vec2, &vec3);
//     let v = utils::cross_product(&t, &u);
//     let cosa = utils::dot_product(&t, &u)/(utils::norm(&t)*utils::norm(&u));
//     let sina = utils::dot_product(&vec2, &v)/(utils::norm(&vec2)*utils::norm(&t)*utils::norm(&u));
//     (sina).atan2(cosa)*180.0/f64::consts::PI
// }

// pub fn get_dihedral_v2(atm1: &[f64;3], atm2: &[f64;3], atm3: &[f64;3], atm4: &[f64;3]) -> f64 {
//     let vec1 = [atm2[0]-atm1[0], atm2[1]-atm1[1], atm2[2]-atm1[2]];
//     let vec2 = [atm3[0]-atm2[0], atm3[1]-atm2[1], atm3[2]-atm2[2]];
//     let vec3 = [atm4[0]-atm3[0], atm4[1]-atm3[1], atm4[2]-atm3[2]];
//     let n1 = utils::cross_product(&vec1, &vec2);
//     let n2 = utils::cross_product(&vec2, &vec3);
//     let d = utils::cross_product(&n1, &vec2);
//     if utils::dot_product(&d, &vec2) < 0.0 {
//         return -utils::angle(&n1, &n2);
//     } else {
//         return utils::angle(&n1, &n2);
//     }
// }

// pub fn get_dihedral(atm1: &[f64;3], atm2: &[f64;3], atm3: &[f64;3], atm4: &[f64;3]) -> f64 {
//     get_dihedral_v1(atm1, atm2, atm3, atm4)
// }


/// Get dihedral angle for atom 1, 2, 3, 4 in degrees
/// 
///           3 - 4
///         /
///   1 - 2
/// 
pub fn get_dihedral_deg(atm1: &CoordVec, atm2: &CoordVec, atm3: &CoordVec, atm4: &CoordVec) -> f64 {
    let vec1 = atm2 - atm1;
    let vec2 = atm3 - atm2;
    let vec3 = atm4 - atm3;
    let t = vec1.cross(&vec2);
    let u = vec2.cross(&vec3);
    let v = t.cross(&u);
    let cosa = t.dot(&u)/(t.norm()*u.norm());
    let sina = vec2.dot(&v)/(vec2.norm()*t.norm()*u.norm());
    (sina).atan2(cosa)*180.0/f64::consts::PI

}

/// Get dihedral angle for atom 1, 2, 3, 4 in radians
/// 
///           3 - 4
///         /
///   1 - 2
/// 
pub fn get_dihedral_rad(atm1: &CoordVec, atm2: &CoordVec, atm3: &CoordVec, atm4: &CoordVec) -> f64 {
    let vec1 = atm2 - atm1;
    let vec2 = atm3 - atm2;
    let vec3 = atm4 - atm3;
    let t = vec1.cross(&vec2);
    let u = vec2.cross(&vec3);
    let v = t.cross(&u);
    let cosa = t.dot(&u)/(t.norm()*u.norm());
    let sina = vec2.dot(&v)/(vec2.norm()*t.norm()*u.norm());
    (sina).atan2(cosa)
}



/// Get bond angle for atom 1, 2, 3
/// 
/// Note: atom 2 is the central atom
pub fn get_angle(atm1: &CoordVec, atm2: &CoordVec, atm3: &CoordVec) -> f64 {

    let vec1 = atm1 - atm2;
    let vec2 = atm3 - atm2;
    return vec1.angle_deg(&vec2);
}

// pub fn get_distance(atm1: &[f64;3], atm2: &[f64;3]) -> f64 {
//     return ((atm1[0]-atm2[0]).powi(2) + (atm1[1]-atm2[1]).powi(2) + (atm1[2]-atm2[2]).powi(2)).sqrt();
// }


#[test]
fn test_dihedral_angle() {
    let atm1 = CoordVec {x: 0.0, y: 0.0, z: 0.0};
    let atm2 = CoordVec {x: 1.0, y: 0.0, z: 0.0};
    let atm3 = CoordVec {x: 1.0, y: 1.0, z: 0.0};
    let atm4 = CoordVec {x: 1.0, y: 1.0, z: 1.0};

    assert_eq!(get_dihedral_deg(&atm1, &atm2, &atm3, &atm4), 90.0);
}

