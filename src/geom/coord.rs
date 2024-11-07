//! Cartesian coordinate structure,
//! including functions for coordinate operations

use core::f64;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use forward_ref::{forward_ref_binop, forward_ref_op_assign};

/// Coordinate vector for coordinate calculations
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct CoordVec {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl CoordVec {
    /// Create a zero `CoordVec`
    pub fn new() -> CoordVec {
        CoordVec {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Create a `CoordVec` from an array
    pub fn from(vec: &[f64;3]) -> CoordVec {
        CoordVec {
            x: vec[0],
            y: vec[1],
            z: vec[2],
        }
    }

    #[inline]
    /// Calculate the distance between two `CoordVec`
    pub fn distance(&self, other: &CoordVec) -> f64 {
        ((self.x-other.x).powi(2) + (self.y-other.y).powi(2) + (self.z-other.z).powi(2)).sqrt()
    }

    /// Return the cross product of two vectors \
    /// Note: Check if the CoordVec is a vector or a point
    #[inline]
    pub fn cross(&self, other: &CoordVec) -> CoordVec {
        CoordVec {
            x: self.y*other.z - self.z*other.y,
            y: self.z*other.x - self.x*other.z,
            z: self.x*other.y - self.y*other.x,
        }
    }

    /// Return the dot product of two vectors \
    /// Note: Check if the CoordVec is a vector or a point
    #[inline]
    pub fn dot(&self, other: &CoordVec) -> f64 {
        self.x*other.x + self.y*other.y + self.z*other.z
    }

    /// Return the norm of a vector |v|
    #[inline]
    pub fn norm(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }
    
    /// Return the square of the norm of a vector v<sup>2</sup>
    #[inline]
    pub fn square(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }

    /// Return the angle between two vectors in degrees
    #[inline]
    pub fn angle_deg(&self, other: &CoordVec) -> f64 {
        let cos_theta = self.dot(other) / (self.norm() * other.norm());
        cos_theta.acos()*180.0/f64::consts::PI
    }

    /// Return the angle between two vectors in radians
    #[inline]
    pub fn angle_rad(&self, other: &CoordVec) -> f64 {
        let cos_theta = self.dot(other) / (self.norm() * other.norm());
        cos_theta.acos()
    }

    #[inline]
    /// Convert `CoordVec` to an array
    pub fn to_array(&self) -> [f64;3] {
        [self.x, self.y, self.z]
    }

}

// Add
impl Add for CoordVec {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

forward_ref_binop!(impl Add, add for CoordVec, CoordVec);

impl AddAssign for CoordVec {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

forward_ref_op_assign!(impl AddAssign, add_assign for CoordVec, CoordVec);

// Sub
impl Sub for CoordVec {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

forward_ref_binop!(impl Sub, sub for CoordVec, CoordVec);

impl SubAssign for CoordVec {
    #[inline]
    fn sub_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

forward_ref_op_assign!(impl SubAssign,sub_assign for CoordVec, CoordVec);

// Mul

impl MulAssign<f64> for CoordVec
{
    #[inline]
    fn mul_assign(&mut self, a: f64) {
        self.x *= a;
        self.y *= a;
        self.z *= a;
    }
}


impl<'a> Mul<f64> for &'a CoordVec
{
    type Output = CoordVec;
    #[inline]
    fn mul(self, a: f64) -> CoordVec {
         CoordVec{
            x: self.x * a,
            y: self.y * a,
            z: self.z * a,
        }
    }
}

impl Mul<f64> for CoordVec
{
    type Output = Self;
    #[inline]
    fn mul(self, a: f64) -> Self::Output {
        Self {
            x: self.x * a,
            y: self.y * a,
            z: self.z * a,
        }
    }
}

// Div

impl DivAssign<f64> for CoordVec
{
    #[inline]
    fn div_assign(&mut self, a: f64) {
        self.x /= a;
        self.y /= a;
        self.z /= a;
    }
}

impl Div<f64> for CoordVec
{
    type Output = Self;
    #[inline]
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
    #[inline]
    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}


/// Get dihedral angle for atom 1, 2, 3, 4 in degrees
/// ```
///           3 - 4
///         /
///   1 - 2
/// ```
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

// /// Get dihedral angle for atom 1, 2, 3, 4 in radians
// /// ```
// ///           3 - 4
// ///         /
// ///   1 - 2
// /// ```
// pub fn get_dihedral_rad(atm1: &CoordVec, atm2: &CoordVec, atm3: &CoordVec, atm4: &CoordVec) -> f64 {
//     let vec1 = atm2 - atm1;
//     let vec2 = atm3 - atm2;
//     let vec3 = atm4 - atm3;
//     let t = vec1.cross(&vec2);
//     let u = vec2.cross(&vec3);
//     let v = t.cross(&u);
//     let cosa = t.dot(&u)/(t.norm()*u.norm());
//     let sina = vec2.dot(&v)/(vec2.norm()*t.norm()*u.norm());
//     (sina).atan2(cosa)
// }


/// Get bond angle for atom 1, 2, 3 \
/// Note: atom 2 is the central atom
pub fn get_angle(atm1: &CoordVec, atm2: &CoordVec, atm3: &CoordVec) -> f64 {

    let vec1 = atm1 - atm2;
    let vec2 = atm3 - atm2;
    return vec1.angle_deg(&vec2);
}


/// format `Vec<CoordVec>` from `Vec<f64>`
pub fn format_coord(coord: Vec<f64>) -> Vec<CoordVec> {
    if coord.len() % 3 != 0 {
        panic!("The length of the coordinate should be a multiple of 3");
    }

    let mut coord_vec = Vec::new();
    let _ = &coord[..].chunks(3).into_iter().for_each(|x| {
        coord_vec.push(CoordVec::from(&[x[0], x[1], x[2]]));
    });

    coord_vec
}



#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_dihedral_angle() {
        let atm1 = CoordVec {x: 0.0, y: 0.0, z: 0.0};
        let atm2 = CoordVec {x: 1.0, y: 0.0, z: 0.0};
        let atm3 = CoordVec {x: 1.0, y: 1.0, z: 0.0};
        let atm4 = CoordVec {x: 1.0, y: 1.0, z: 1.0};   

        assert_eq!(get_dihedral_deg(&atm1, &atm2, &atm3, &atm4), 90.0);
    }   
    

    #[test]
    fn test_format_coord() {
        let coord = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0];
        let coord_vec = format_coord(coord);
        assert_eq!(coord_vec[0], CoordVec {x: 0.0, y: 0.0, z: 0.0});
        assert_eq!(coord_vec[1], CoordVec {x: 1.0, y: 0.0, z: 0.0});
        assert_eq!(coord_vec[2], CoordVec {x: 1.0, y: 1.0, z: 0.0});
    }

}