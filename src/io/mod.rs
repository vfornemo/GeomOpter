//! Input file parsing functionality

use std::fs;
use toml;
use serde::Deserialize;
use log::{error, info};

/// Input structure parsed from input file
/// # Fields:
/// * `path`: file location
/// * `calc_type`: calculation type
/// * `max_cycle`: maximum number of cycles
/// * `rms_tol`: rms tolerance
/// * `car_conver`: cartesian convergence
#[derive(Deserialize, Debug, Clone)]
pub struct Input {
    /// file location
    pub path: String, 
    /// calculation type
    pub calc_type: String,
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
        let input: Input = toml::from_str(&file).unwrap();
        input.input_check();
        input
    }

    pub fn logger(&self) {
        info!("Input file path: {}", self.path);
        info!("Calculation type: {}", self.calc_type);
        info!("Maximum number of cycles: {}", self.max_cycle);
        info!("RMS tolerance: {}", self.rms_tol);
        info!("Cartesian convergence: {}", self.car_conver);
    }

    fn input_check(&self) {
        // check path
        if !std::path::Path::new(&self.path).exists() {
            error!("Error: The mol file does not exist.");
            panic!("The mol file does not exist.");
        }

        // check calc_type
        let calc_types = ["energy", "car", "intl"];
        if !calc_types.contains(&self.calc_type.as_str()) {
            error!("Error: Invalid calculation type, please check input file.");
            panic!("Invalid calculation type, please check input file.");
        }


    }
}
