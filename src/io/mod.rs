
use std::fs;
use toml;
use serde::Deserialize;
use log::{info, debug};

#[derive(Deserialize, Debug, Clone)]
pub struct Input {
    pub path: String,  // file location
    pub max_cycle: usize, // maximum number of cycles
    pub rms_tol: f64, // rms tolerance
    pub car_conver: f64, // cartesian convergence
}


impl Input {

    pub fn from_file(file_path: &str) -> Self {
        let file = fs::read_to_string(file_path).expect("Unable to read input file");

        let input= toml::from_str(&file).unwrap();

        input
    }

    pub fn logger(&self) {
        info!("Input file path: {}", self.path);
        info!("Maximum number of cycles: {}", self.max_cycle);
        info!("RMS tolerance: {}", self.rms_tol);
        info!("Cartesian convergence: {}", self.car_conver);

    }
}
