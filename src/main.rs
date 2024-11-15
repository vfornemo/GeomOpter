mod grad;
mod geom;
mod utils;
mod matrix;
mod opt;
mod io;

use geom::Geom;
use grad::Gradient;
use opt::opt_car::OptCar;
use opt::opt_intl::OptIntl;

use log::info;
use chrono;
use log4rs;
use log4rs::config::{Appender, Config, Root};
use log4rs::append::file::FileAppender;
use log4rs::encode::pattern::PatternEncoder;

fn main() {
    let cmd_line = std::env::args();
    let input_file = cmd_line.skip(1).next().expect("No input file provided. Please use `GeomOpter [input_file]`");
    let config = io::Config::from_file(&input_file);
    let input = config.to_input();
    let mut res = io::Result::new(input.clone());

    // initialize logger
    // log4rs::init_file("./config/log4rs.yaml", Default::default()).unwrap();

    let logfile = FileAppender::builder()
    .encoder(Box::new(PatternEncoder::new("{m}{n}"))).append(false)
    .build(&input.output).unwrap();

    let log_config = Config::builder()
        .appender(Appender::builder().build("logfile", Box::new(logfile)))
        .build(Root::builder()
                   .appender("logfile")
                   .build(input.log_level)).unwrap();

    log4rs::init_config(log_config).unwrap();

    info!("\n**************** Input File ****************");
    input.logger();

    let time = std::time::Instant::now();
    let start = chrono::Local::now();

    match input.calc_type.as_str() {
        "energy" => {
            info!("\n**************** Energy Calculation ****************");
            let mut geom = Geom::from_mol2(input);
            geom.build();
            geom.logger();
            res.e_tot = Some(geom.e_tot);
        }
        "car" => {
            info!("\n**************** Cartesian Optimization ****************");
            let mut geom = Geom::from_mol2(input);
            geom.build();
            geom.logger();
            let mut grad = Gradient::from(geom);
            grad.build();
            grad.logger();
            let mut opt = OptCar::from(grad.clone());
            opt.build();
            opt.logger();
            opt.kernel(&mut res);
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
            opt_intl.kernel(&mut res);
        }
        _ => (),
    }

    let t_tot = time.elapsed().as_millis();
    let end = chrono::Local::now();
    info!("\n**************** Summary ****************");
    res.logger();
    info!("Job start at: {:?}", start);
    info!("Job end at: {:?}", end);
    info!("Time Usage: {} ms", t_tot);

}
