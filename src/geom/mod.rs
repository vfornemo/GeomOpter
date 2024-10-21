mod coord;
use std::backtrace;
use std::fs::File;
use std::io::{BufReader, BufRead};
use crate::data;

#[derive(Debug)]
pub enum BondType {
    CH,
    CC,
    UNK,
}

#[derive(Debug)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub bond_type: BondType,
}


impl Bond {
    fn new() -> Bond {
        Bond {
            atom1: 0usize,
            atom2: 0usize,
            bond_type: BondType::UNK,
        }
    }
    
    fn from(atm1: usize, atm2: usize) -> Bond {
        Bond {
            atom1: atm1,
            atom2: atm2,
            bond_type: BondType::UNK,
        }
    }

    fn get_bond_type(&mut self, atoms_idx: &Vec<usize>) {
        if atoms_idx[self.atom1-1] == 6 && atoms_idx[self.atom2-1] == 6 {
            self.bond_type = BondType::CC;
        } else {
            self.bond_type = BondType::CH;
        }
    }
    
}

#[derive(Debug)]
pub struct Geom {
    pub atm_num: usize,  // Number of atoms
    pub atoms: Vec<String>,  // Atom sequence with element symbols
    pub atoms_idx: Vec<usize>,  // Atom sequence with atom numbers
    pub coord: Vec<[f64;3]>,  // Cartesian coordinates
    pub bond_num: usize,  // Number of bonds
    pub bonds: Vec<Bond>,  // Bonds
    pub c_number: usize,  // Number of carbons
    pub c_neighbor: Vec<Vec<usize>>,  // Neighboring atoms of each carbon
    pub coord_intnl: Vec<[f64;3]>, // Internal coordinates
}

impl Geom {
    pub fn new() -> Geom {
        Geom {
            atm_num: 0usize,
            atoms: Vec::new(),
            atoms_idx: Vec::new(),
            coord: Vec::new(),
            bond_num: 0usize,
            bonds: Vec::new(),
            c_number: 0usize,
            c_neighbor: Vec::new(),
            coord_intnl: Vec::new(),
        }
    }

    
    pub fn from_mol2(mol2: &str) -> Geom {
        let mut geom = Geom::new();
        // open file
        let file = File::open(mol2).expect("File not found");
        let mut reader = BufReader::new(file).lines();
        let line1 = reader.next().unwrap().unwrap();
        
        // first line of mol2 file:
        // number of atoms, number of bonds, number of substructures, number of features, number of sets
        // we only consider the first two numbers
        let mut iter = line1.split_whitespace();
        geom.atm_num = iter.next().unwrap().parse::<usize>().unwrap();
        geom.bond_num = iter.next().unwrap().parse::<usize>().unwrap();
        print!("atm_num: {:?}\n", geom.atm_num);
        print!("bond_num: {:?}\n", geom.bond_num);

        // read atoms
        for _ in 0..geom.atm_num {
            let line = reader.next().unwrap().unwrap();
            let mut iter = line.split_whitespace();
            let x = iter.next().unwrap().parse::<f64>().unwrap();
            let y = iter.next().unwrap().parse::<f64>().unwrap();
            let z = iter.next().unwrap().parse::<f64>().unwrap();
            let atom = iter.next().unwrap();
            geom.atoms.push(atom.to_string());
            geom.coord.push([x, y, z]);
        }
        geom.atoms_idx = data::elems_to_idx(&geom.atoms);

        for i in 0..geom.atm_num {
            if geom.atoms_idx[i] == 6 {
                geom.c_number += 1;
            }
        }

        print!("atoms: {:?}\n", geom.atoms);
        print!("atoms_idx: {:?}\n", geom.atoms_idx);
        print!("coord: {:?}\n", geom.coord);

        // read bonds
        for _ in 0..geom.bond_num {
            let line = reader.next().unwrap().unwrap();
            let mut iter = line.split_whitespace();
            let atm1 = iter.next().unwrap().parse::<usize>().unwrap();
            let atm2 = iter.next().unwrap().parse::<usize>().unwrap();
            let mut bond = Bond::from(atm1, atm2);
            bond.get_bond_type(&geom.atoms_idx);
            geom.bonds.push(bond);
        }
        print!("number of bonds: {:?}\n", geom.bond_num);
        print!("bonds: {:?}\n", geom.bonds);

        geom.get_neighbor();

        println!("c_number: {:?}", geom.c_number);
        println!("c_neighbor: {:?}", geom.c_neighbor);
        
        geom
    }

    fn get_neighbor(&mut self) {
        // construct neighbor table for saturated hydrocarbons
        //        H   H
        //    H - C - C - H
        //        H   H
        self.c_neighbor = vec![vec![]; self.c_number];

        for bond in &self.bonds {
            match bond.bond_type {
                BondType::CC => { 
                    self.c_neighbor[bond.atom1-1].push(bond.atom2); self.c_neighbor[bond.atom2-1].push(bond.atom1);
                },
                BondType::CH => {
                    if self.atoms_idx[bond.atom1-1] == 6 {
                        self.c_neighbor[bond.atom1-1].push(bond.atom2);
                    } else {
                        self.c_neighbor[bond.atom2-1].push(bond.atom1);
                    }
                },
                BondType::UNK => (),
            }
        }

    }

    // construct internal coordinate table for saturated hydrocarbons
    //        H   H
    //    H - C - C - H
    //        H   H
    // pub fn intnl_table(&self) -> Vec<[usize;4]> {

    
    // }

}