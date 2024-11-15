mod grad;
mod geom;
mod utils;
mod matrix;
mod opt;
mod io;

use geom::Geom;
use log::info;
use grad::Gradient;
use opt::opt_car::OptCar;
use log4rs;
use crate::opt::opt_intl::OptIntl;

fn main() {

    let cmd_line = std::env::args();
    let input_file = cmd_line.skip(1).next().expect("No input file provided. Please use `GeomOpter [input_file]`");
    let input = io::Input::from_file(&input_file);

    // initialize logger
    log4rs::init_file("./config/log4rs.yaml", Default::default()).unwrap();

    info!("\n**************** Input ****************");
    input.logger();

    let time = std::time::Instant::now();

    match input.calc_type.as_str() {
        "energy" => {
            info!("\n**************** Energy Calculation ****************");
            let mut geom = Geom::from_mol2(input);
            geom.build();
            geom.logger();
        }
        "car" => {
            info!("\n**************** Cartesian Optimization ****************");
            let mut geom = Geom::from_mol2(input);
            geom.build();
            geom.logger();
            let mut grad = Gradient::from(geom);
            grad.build();
            let mut opt = OptCar::from(grad.clone());
            opt.build();
            opt.logger();
            opt.kernel();
        }
        "intl" => {
            info!("\n**************** Internal Optimization ****************");
            let mut geom = Geom::from_mol2(input);
            geom.build();
            geom.logger();
            let mut grad = Gradient::from(geom);
            grad.build();
            let mut opt_intl = OptIntl::from(grad);
            opt_intl.logger();
            opt_intl.build();
            opt_intl.kernel();
        }
        _ => (),
    }

    let t_tot = time.elapsed().as_secs_f64();

    info!("\n**************** Summary ****************");

    info!("Time Usage: {:.4} s", t_tot);

}
