use std::f64;


pub fn cross_product(v1: &[f64;3], v2: &[f64;3]) -> [f64;3] {
    [v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]]
}

pub fn dot_product(v1: &[f64;3], v2: &[f64;3]) -> f64 {
    v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
}

pub fn norm(v: &[f64;3]) -> f64 {
    (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt()
}

pub fn angle(v1: &[f64;3], v2: &[f64;3]) -> f64 {
    let cos_theta = dot_product(v1, v2) / (norm(v1) * norm(v2));
    cos_theta.acos()*180.0/f64::consts::PI
}