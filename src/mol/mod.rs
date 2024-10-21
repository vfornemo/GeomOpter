
use crate::geom::Geom;



#[derive(Debug)]
pub struct Molecule {
    pub geom: Geom,

    // ...
}

impl Molecule {
    fn new() -> Molecule {
        Molecule {
            geom: Geom::new(),


        }
    }

}

