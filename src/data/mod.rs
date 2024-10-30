#![allow(non_upper_case_globals)]

// Hashmap of all element and their atomic number
pub fn elem_to_idx(elem: &str) -> usize {
    match elem {
        "H" => 1,
        "C" => 6,
        _ => 0,
    }
    
}

pub fn elems_to_idx(elems: &Vec<String>) -> Vec<usize> {
    let mut idxs = Vec::new();
    for elem in elems {
        idxs.push(elem_to_idx(elem));
    }
    idxs
}

    
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
pub const WOLFE_C1: f64 = 0.1;
pub const MAX_STEP_SIZE: f64 = 0.02;