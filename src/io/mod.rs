//! Input file parsing functionality

use std::fs;
use toml;
use serde::Deserialize;
use log::{error, info, warn, LevelFilter};

use crate::geom::Geom;


/// Config structure parsed from .toml file
#[derive(Deserialize, Debug, Clone)]
pub struct Config {
    /// name of the job
    pub name: Option<String>,
    /// file location
    pub path: Option<String>, 
    /// calculation type
    pub calc_type: Option<String>,
    /// maximum number of cycles
    pub max_cycle: Option<usize>,
    /// rms tolerance
    pub rms_tol: Option<f64>,
    /// cartesian convergence
    pub car_conver: Option<f64>,
    /// log level
    pub log_level: Option<String>,
    /// output file path
    pub output: Option<String>,

}

impl Config {
    pub fn from_file(file_path: &str) -> Self {
        let file = fs::read_to_string(file_path).expect("Unable to read input file");
        let config: Config = toml::from_str(&file).unwrap();
        config
    }

    pub fn to_input(self) -> Input {
        let res = Input {
            name: self.name.unwrap_or(String::from("my_job")),
            path: self.path.expect("Error: No mol file provided in input file"),
            calc_type: self.calc_type.expect("Error: No calculation type provided in input file"),
            max_cycle: self.max_cycle.unwrap_or(200),
            rms_tol: self.rms_tol.unwrap_or(0.001),
            car_conver: self.car_conver.unwrap_or(0.00001),
            log_level: match self.log_level.as_deref() {
                Some("info") => LevelFilter::Info,
                Some("warn") => LevelFilter::Warn,
                Some("error") => LevelFilter::Error,
                Some("debug") => LevelFilter::Debug,
                Some("trace") => LevelFilter::Trace,
                Some("off") => LevelFilter::Off,
                _ => LevelFilter::Info,
            },
            output: self.output.unwrap_or(String::from("./output.log")),
        };
        res.input_check();
        res
    }

}

/// Input structure for the calculation
/// # Fields:
/// * `name`: name of the job
/// * `path`: file location
/// * `calc_type`: calculation type
/// * `max_cycle`: maximum number of cycles
/// * `rms_tol`: rms tolerance
/// * `car_conver`: cartesian convergence
/// * `log_level`: log level
/// * `output`: output file path
#[derive(Deserialize, Debug, Clone)]
pub struct Input {
    /// name of the job
    pub name: String,
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
    /// log level
    pub log_level: LevelFilter,
    /// output file path
    pub output: String,

}

impl Input {
    /// check the validity of input file
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

    pub fn logger(&self) {
        info!("Job name: {}", self.name);
        info!("Input file path: {}", self.path);
        info!("Calculation type: {}", self.calc_type);
        info!("Maximum number of cycles: {}", self.max_cycle);
        info!("RMS tolerance: {}", self.rms_tol);
        info!("Cartesian convergence: {}", self.car_conver);
        info!("Log level: {:?}", self.log_level);
        info!("Output file path: {}", self.output);
    }


}


pub struct Result {
    pub input: Input,
    pub is_converged: bool,
    pub conv_cyc: Option<usize>,
    pub grms: Option<f64>,
    pub e_tot: Option<f64>,
    pub geom: Option<Geom>,
}

impl Result {
    pub fn new(input: Input) -> Self {
        Result {
            input,
            is_converged: false,
            conv_cyc: None,
            grms: None,
            e_tot: None,
            geom: None,
        }
    }

    pub fn logger(&self) {
        let print_digits = (-self.input.rms_tol.log10()).ceil() as usize + 4;

        match self.input.calc_type.as_str() {
            "energy" => {
                info!("Energy calculation of {}", self.input.name);
                info!("Energy of current structure: {:-.8}", self.e_tot.unwrap());
                return;
        },
            "car" => info!("Cartesian optimization of {}", self.input.name),
            "intl" => info!("Internal optimization of {}", self.input.name),
            _ => (),
        }

        match self.is_converged {
            true => {
                info!("Optimization converged within {} cycles", self.conv_cyc.unwrap());
                match self.e_tot {
                    Some(e) => info!("Final energy at mimimum: {:-.1$}", e, print_digits),
                    None => (),
                }
                match &self.geom {
                    Some(g) => {
                        info!("Final coordinates:");
                        g.formated_output_coord("info");
                    },
                    None => (),
                }
            }
            false => {
                warn!("Warning: Optimization did not converge within {} cycles", self.input.max_cycle);
                warn!("The final GRMS is {:-.6}, while the tolerance is {}", self.grms.unwrap(), self.input.rms_tol);
                warn!("Please check the input and try again");
            }
        }


    }
    
}



#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn check_config() {
        let config = Config::from_file("config/ctrl.toml");
        println!("{:?}", config);

    }


}