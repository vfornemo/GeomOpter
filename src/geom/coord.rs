use core::f64;

use crate::utils::{self, angle};


/// Get dihedral angle for atom 1, 2, 3, 4
/// 
///           3 - 4
///         /
///   1 - 2
/// 
pub fn get_dihedral_v1(atm1: &[f64;3], atm2: &[f64;3], atm3: &[f64;3], atm4: &[f64;3]) -> f64 {
    let vec1 = [atm2[0]-atm1[0], atm2[1]-atm1[1], atm2[2]-atm1[2]];
    let vec2 = [atm3[0]-atm2[0], atm3[1]-atm2[1], atm3[2]-atm2[2]];
    let vec3 = [atm4[0]-atm3[0], atm4[1]-atm3[1], atm4[2]-atm3[2]];
    let t = utils::cross_product(&vec1, &vec2);
    let u = utils::cross_product(&vec2, &vec3);
    let v = utils::cross_product(&t, &u);
    let cosa = utils::dot_product(&t, &u)/(utils::norm(&t)*utils::norm(&u));
    let sina = utils::dot_product(&vec2, &v)/(utils::norm(&vec2)*utils::norm(&t)*utils::norm(&u));
    (sina).atan2(cosa)*180.0/f64::consts::PI
}

pub fn get_dihedral_v2(atm1: &[f64;3], atm2: &[f64;3], atm3: &[f64;3], atm4: &[f64;3]) -> f64 {
    let vec1 = [atm2[0]-atm1[0], atm2[1]-atm1[1], atm2[2]-atm1[2]];
    let vec2 = [atm3[0]-atm2[0], atm3[1]-atm2[1], atm3[2]-atm2[2]];
    let vec3 = [atm4[0]-atm3[0], atm4[1]-atm3[1], atm4[2]-atm3[2]];
    let n1 = utils::cross_product(&vec1, &vec2);
    let n2 = utils::cross_product(&vec2, &vec3);
    let d = utils::cross_product(&n1, &vec2);
    if utils::dot_product(&d, &vec2) < 0.0 {
        return -utils::angle(&n1, &n2);
    } else {
        return utils::angle(&n1, &n2);
    }
}

pub fn get_dihedral(atm1: &[f64;3], atm2: &[f64;3], atm3: &[f64;3], atm4: &[f64;3]) -> f64 {
    get_dihedral_v1(atm1, atm2, atm3, atm4)
}



/// Get bond angle for atom 1, 2, 3
/// 
/// Note: atom 2 is the central atom
pub fn get_angle(atm1: &[f64;3], atm2: &[f64;3], atm3: &[f64;3]) -> f64 {

    let vec1 = [atm1[0]-atm2[0], atm1[1]-atm2[1], atm1[2]-atm2[2]];
    let vec2 = [atm3[0]-atm2[0], atm3[1]-atm2[1], atm3[2]-atm2[2]];
    return utils::angle(&vec1, &vec2);
}

pub fn get_distance(atm1: &[f64;3], atm2: &[f64;3]) -> f64 {
    return ((atm1[0]-atm2[0]).powi(2) + (atm1[1]-atm2[1]).powi(2) + (atm1[2]-atm2[2]).powi(2)).sqrt();
}


#[test]
fn test_dihedral_angle() {
    let atm1 = [0.0, 0.0, 0.0];
    let atm2 = [1.0, 0.0, 0.0];
    let atm3 = [1.0, 1.0, 0.0];
    let atm4 = [1.0, 1.0, 1.0];
    assert_eq!(get_dihedral(&atm1, &atm2, &atm3, &atm4), 90.0);
}

