#![allow(non_upper_case_globals)]

use std::collections::HashMap;

// Hashmap of all element and their atomic number
pub fn elem_to_idx(elem: &str) -> usize {
    let elem_list:HashMap<String, usize> = HashMap::from([
        (String::from("H"),1),
        (String::from("He"),2),
        (String::from("Li"),3),
        (String::from("Be"),4),
        (String::from("B"),5),
        (String::from("C"),6),
        (String::from("N"),7),
        (String::from("O"),8)]);

        elem_list.get(elem).unwrap().clone()
    
}

pub fn elems_to_idx(elems: &Vec<String>) -> Vec<usize> {
    let mut idxs = Vec::new();
    for elem in elems {
        idxs.push(elem_to_idx(elem));
    }
    idxs
}

// [charge, mass, name]
const ELEM_INFO: [[f64; 2]; 9] = [
    [0.0, 0.0],         // 0
    [1.0, 1.00784],     // 1
    [2.0, 4.002602],   // 2
    [3.0, 6.938],      // 3
    [4.0, 9.0121831],  // 4
    [5.0, 10.81],       // 5
    [6.0, 12.011],      // 6
    [7.0, 14.007],      // 7
    [8.0, 15.999]];    // 8
    
pub const R0_CC: f64 = 1.5300;
pub const R0_CH: f64 = 1.1100;
pub const Kb_CC: f64 = 300.0;
pub const Kb_CH: f64 = 350.0;
pub const A0_CCC: f64 = 109.50;
pub const A0_XCX: f64 = 109.50;
pub const Ka_CCC: f64 = 60.00;
pub const Ka_XCX: f64 = 35.00;
pub const N_XCCX: f64 = 3.0;
pub const A_XCCX: f64 = 0.300;
pub const EPSILON_H: f64 = 0.0300;
pub const EPSILON_C: f64 = 0.0700;
pub const SIGMA_H: f64 = 1.20000;
pub const SIGMA_C: f64 = 1.75000;



