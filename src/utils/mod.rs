use std::f64;

/// Return the cross product of two vectors
/// v1 x v2
pub fn cross_product(v1: &[f64;3], v2: &[f64;3]) -> [f64;3] {
    [v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]]
}

/// Return the dot product of two vectors 
/// v1 . v2
pub fn dot_product(v1: &[f64;3], v2: &[f64;3]) -> f64 {
    v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
}

/// Return the norm of a vector 
/// |v|
pub fn norm(v: &[f64;3]) -> f64 {
    (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt()
}

/// Return the angle between two vectors in degrees
pub fn angle(v1: &[f64;3], v2: &[f64;3]) -> f64 {
    let cos_theta = dot_product(v1, v2) / (norm(v1) * norm(v2));
    cos_theta.acos()*180.0/f64::consts::PI
}