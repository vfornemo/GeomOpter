mod mol;
mod data;
mod geom;
mod utils;
mod matrix;

use geom::Geom;
use mol::Molecule;

fn main() {
    let file_path = "examples/ethane.mol2";
    let mut geom = Geom::from_mol2(file_path);
    geom.build();
    // geom.logger();
    let mut mol = Molecule::from(geom);
    mol.build();
    mol.logger();

}
