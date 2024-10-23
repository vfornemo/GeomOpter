
use crate::geom::Geom;



#[derive(Debug)]
pub struct Molecule {
    pub geom: Geom,
    pub E_potential: f64,


    // ...
}

impl Molecule {
    fn new() -> Molecule {
        Molecule {
            geom: Geom::new(),
            E_potential: 0.0f64,
        }
    }

    fn from(geom: Geom) -> Molecule {
        Molecule {
            geom,
            E_potential: 0.0f64,
        }
    }

    fn build(&mut self) {
    }


    fn get_potential_energy(&self) {
    }


}

