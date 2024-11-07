//! Input file parsing functionality

use std::fs;
use toml;
use serde::Deserialize;
use log::info;

/// Input structure parsed from input file
/// # Fields:
/// * `path`: file location
/// * `max_cycle`: maximum number of cycles
/// * `rms_tol`: rms tolerance
/// * `car_conver`: cartesian convergence
#[derive(Deserialize, Debug, Clone)]
pub struct Input {
    /// file location
    pub path: String, 
    /// maximum number of cycles
    pub max_cycle: usize, 
    /// rms tolerance
    pub rms_tol: f64, 
    /// cartesian convergence
    pub car_conver: f64, 
}


impl Input {
    /// Parse input file
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
