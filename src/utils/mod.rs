use log::{debug, info, trace};


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
