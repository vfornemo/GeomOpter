mod mol;
mod geom;
mod utils;
mod matrix;
mod opt;
mod io;

use geom::Geom;
use log::info;
use mol::Molecule;
use opt::opt_car::OptCar;
use log4rs;
use crate::opt::opt_intl::OptIntl;

fn main() {
    let time = std::time::Instant::now();
    let mut t_tot = 0.0;

    let cmd_line = std::env::args();
    let input_file = cmd_line.skip(1).next().unwrap();
    let input = io::Input::from_file(&input_file);

    // initialize logger
    log4rs::init_file("./config/log4rs.yaml", Default::default()).unwrap();

    info!("\n**************** Input ****************");
    input.logger();

    let mut geom = Geom::from_mol2(input);
    geom.build();
    info!("\n**************** Initialize ****************");
    geom.logger();
    let mut mol = Molecule::from(geom);
    mol.build();
    mol.logger();
    let t_init = time.elapsed().as_secs_f64();
    t_tot += t_init;

    // info!("\n**************** Cartesian Optimization ****************");
    // mol.geom.logger();
    // let mut opt = OptCar::from(mol.clone());
    // opt.build();
    // opt.logger();
    // opt.kernel();
    // let t_opt = time.elapsed().as_secs_f64();
    // t_tot += t_opt;

    info!("\n**************** Internal Optimization ****************");
    mol.geom.logger();
    let mut opt_intl = OptIntl::from(mol);
    opt_intl.logger();
    opt_intl.build();
    opt_intl.kernel();
    let t_opt_intl = time.elapsed().as_secs_f64();
    t_tot += t_opt_intl;


    info!("\n**************** Summary ****************");

    info!("Time Usage: {:.4} s", t_tot);
    info!("Initialization time: {:.4} s", t_init);
    // info!("Cartesian Optimization time: {:.4} s", t_opt);
    info!("Internal Optimization time: {:.4} s", t_opt_intl);


}
