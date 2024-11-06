use log::{debug, info, trace};
use std::f64::consts::PI;


pub fn rms(vec: &Vec<f64>) -> f64 {
    (vec.iter().map(|x| x.powi(2)).sum::<f64>()/vec.len() as f64).sqrt()
}

pub fn formated_output_vec(vec: &Vec<f64>,  precision: usize, print_level: &str) {
    let mut output = String::new();
    for i in 0..vec.len() {
        output.push_str(&format!("{:-11.*}", precision, vec[i]));
    }
    output.push_str("\n");

    match print_level {
        "debug" => debug!("{}", output),
        "info" => info!("{}", output),
        "trace" => trace!("{}", output),
        "off" => (),
        _ => println!("{}", output),
    }
}



pub fn dihedral_diff_check_rad(di: &mut f64) {
    let _2pi = 2.0*PI;
    if *di > _2pi {
        *di -= _2pi;
        dihedral_diff_check_rad(di);
    } else if *di > PI {
        *di -= _2pi;
        return
    } else if *di < -_2pi {
        *di += _2pi;
        dihedral_diff_check_rad(di);
    } else if *di < -PI {
        *di += _2pi;
        return
    }  else {
        return
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rms() {
        let vec = vec![-5.597918, -14.720191, -10.062227, 36.805617, -19.415319, -2.191495, 6.600580, 2.679783, 0.403894, -2.949445, 0.294711, -6.247580, 5.066267, 3.039453, -1.049026, 0.635902, 0.151318, -1.282693, -1.808231, -2.982609, -4.473819, -2.097598, -0.182994, -1.674205, 0.702017, 4.417356, 2.926145, 5.302366];
        println!("{:?}", rms(&vec));
        let vec2 = vec![  -10.464441,   -31.812504,   -21.746560,    79.437265,   -41.923997,    -4.770806,    14.261180,     3.283227,     0.220859,    -3.730432,     0.301442,    -7.150790,     6.192046,     3.742520,    -1.358443,     0.844431,     0.102716,    -1.700123,    -2.066855,    -3.503042,    -5.197075,    -2.466898,    -0.252375,    -1.946408,     0.783769,     5.053091,     3.359058,     6.089235 
        ];
        println!("{:?}", rms(&vec2));
    }

    #[test]
    fn test_rad() {
        let mut di = 6.894;
        dihedral_diff_check_rad(&mut di);
        println!("{:?}", di);
        let mut di = 15.26;
        dihedral_diff_check_rad(&mut di);
        println!("{:?}", di);
        let mut di = 4.25;
        dihedral_diff_check_rad(&mut di);
        println!("{:?}", di);
        let mut di = -15.26;
        dihedral_diff_check_rad(&mut di);
        println!("{:?}", di);

    }

}