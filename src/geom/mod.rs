mod coord;
use std::f64::consts::SQRT_2;
use std::fs::File;
use std::io::{BufReader, BufRead};
use crate::data::{self, EPSILON_C};
use std::collections::VecDeque;

#[derive(Debug)]
pub enum BondType {
    CH,
    CC,
    HH,
    XX,
}

#[derive(Debug)]
pub struct Bond {
    pub atms: [usize;2],
    pub bond_type: BondType,
    pub bond_length: f64,
    pub e_str: f64,  // stretching energy (or vdw energy for atom pairs)
}


impl Bond {
    fn new() -> Bond {
        Bond {
            atms: [0usize;2],
            bond_type: BondType::XX,
            bond_length: 0.0f64,
            e_str: 0.0f64,
        }
    }
    
    fn from(atm1: usize, atm2: usize) -> Bond {
        Bond {
            atms: [atm1, atm2],
            bond_type: BondType::XX,
            bond_length: 0.0f64,
            e_str: 0.0f64,
        }
    }

    fn get_bond_type(&mut self, atoms_idx: &Vec<usize>) {
        let mut c_num = 0;
        for i in 0..2 {
            if atoms_idx[self.atms[i]-1] == 6 {
                c_num += 1;
            }
        }

        if c_num == 2 {
            self.bond_type = BondType::CC;
        } else if c_num == 1 {
            self.bond_type = BondType::CH;
        } else {
            self.bond_type = BondType::HH;
        }
    }

    fn get_bond_length(&mut self, coord: &Vec<[f64;3]>) {
        self.bond_length = coord::get_distance(&coord[self.atms[0]-1], &coord[self.atms[1]-1]);
    }

    /// Calculate stretching energy of a bond
    fn get_e_str(&mut self) {
        self.e_str = match self.bond_type {
            BondType::CC => {data::Kb_CC*(self.bond_length - data::R0_CC).powi(2)},
            BondType::CH => {data::Kb_CH*(self.bond_length - data::R0_CH).powi(2)},
            BondType::HH => 0.0,
            BondType::XX => 0.0,
        }

    }

    fn get_e_vdw(&mut self) {
        let epsilon_CH = (data::EPSILON_C*data::EPSILON_H).sqrt();
        let sigma_CH = 2.0*(data::SIGMA_C*data::SIGMA_H).sqrt();
        let sigma_HH = 2.0*data::SIGMA_H;
        let sigma_CC = 2.0*data::SIGMA_C;

        self.e_str = match self.bond_type {
            BondType::CH => {4.0*epsilon_CH*((sigma_CH/self.bond_length).powi(12) - (sigma_CH/self.bond_length).powi(6))},
            BondType::HH => {4.0*data::EPSILON_H*((sigma_HH/self.bond_length).powi(12) - (sigma_HH/self.bond_length).powi(6))},
            BondType::CC => {4.0*data::EPSILON_C*((sigma_CC/self.bond_length).powi(12) - (sigma_CC/self.bond_length).powi(6))},
            BondType::XX => 0.0,
        }

    }
    
}


#[derive(Debug)]
pub enum AngleType {
    CCC,
    XCX,
    XXX,
}


#[derive(Debug)]
pub struct Angle {
    pub atms: [usize;3],
    pub angle_type: AngleType,
    pub angle: f64,
    pub e_bend: f64,  // bending energy
}

impl Angle {
    fn new() -> Angle {
        Angle {
            atms: [0usize;3],
            angle_type: AngleType::XXX,
            angle: 0.0f64,
            e_bend: 0.0f64,
        }
    }

    fn from(atms: [usize;3]) -> Angle {
        Angle {
            atms,
            angle_type: AngleType::XXX,
            angle: 0.0f64,
            e_bend: 0.0f64,
        }
    }

    fn get_angle(&mut self, coord: &Vec<[f64;3]>) {
        self.angle = coord::get_angle(&coord[self.atms[0]-1], &coord[self.atms[1]-1], &coord[self.atms[2]-1]);
    }

    fn get_angle_type(&mut self, atoms_idx: &Vec<usize>) {
        if atoms_idx[self.atms[0]-1] == 6 && atoms_idx[self.atms[1]-1] == 6 && atoms_idx[self.atms[2]-1] == 6 {
            self.angle_type = AngleType::CCC;
        } else {
            self.angle_type = AngleType::XCX;
        }
    }

    fn get_e_bend(&mut self) {
        self.e_bend = match self.angle_type {
            AngleType::CCC => {data::Ka_CCC*(self.angle - data::A0_CCC).to_radians().powi(2)},
            AngleType::XCX => {data::Ka_XCX*(self.angle - data::A0_XCX).to_radians().powi(2)},
            AngleType::XXX => 0.0,
        }
    }

    
}


#[derive(Debug)]
pub struct Dihedral {
    pub atms: [usize;4],
    pub dihedral: f64,
    pub e_tor: f64,  // torsion energy
}

impl Dihedral {
    fn new() -> Dihedral {
        Dihedral {
            atms: [0usize;4],
            dihedral: 0.0f64,
            e_tor: 0.0f64,
        }
    }

    fn from(atms: [usize;4]) -> Dihedral {
        Dihedral {
            atms: atms,
            dihedral: 0.0f64,
            e_tor: 0.0f64,
        }
    }

    fn get_dihedral(&mut self, coord: &Vec<[f64;3]>) {
        self.dihedral = coord::get_dihedral(&coord[self.atms[0]-1], &coord[self.atms[1]-1], &coord[self.atms[2]-1], &coord[self.atms[3]-1]);
    }

    fn get_e_tor(&mut self) {
        self.e_tor = data::A_XCCX*(1.0 + (data::N_XCCX*self.dihedral.to_radians()).cos());
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
    pub etol_str: f64,  // Total Stretching energy
    pub angle_num: usize,  // Number of angles
    pub angles: Vec<Angle>,  // Angles
    pub etol_bend: f64,  // Total Bending energy
    pub dihedral_num: usize,  // Number of dihedral angles
    pub dihedrals: Vec<Dihedral>,  // Dihedral angles
    pub etol_tor: f64,  // Total Torsion energy
    pub c_number: usize,  // Number of carbons
    pub neighbors: Vec<Vec<usize>>,  // Neighboring atoms of each atom
    pub atm_pair: Vec<Bond>,  // Atoms pairs for LJ calculation
    pub etol_vdw: f64,  // Total VDW energy
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
            etol_str: 0.0f64,
            angle_num: 0usize,
            angles: Vec::new(),
            etol_bend: 0.0f64,
            dihedral_num: 0usize,
            dihedrals: Vec::new(),
            etol_tor: 0.0f64,
            c_number: 0usize,
            neighbors: Vec::new(),
            atm_pair: Vec::new(),
            etol_vdw: 0.0f64,
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


        // read bonds
        for _ in 0..geom.bond_num {
            let line = reader.next().unwrap().unwrap();
            let mut iter = line.split_whitespace();
            let atm1 = iter.next().unwrap().parse::<usize>().unwrap();
            let atm2 = iter.next().unwrap().parse::<usize>().unwrap();
            let bond = Bond::from(atm1, atm2);
            geom.bonds.push(bond);
        }

        geom.build();
        geom.logger();

        geom
    }

    fn build(&mut self) {

        self.atoms_idx = data::elems_to_idx(&self.atoms);

        for i in 0..self.atm_num {
            if self.atoms_idx[i] == 6 {
                self.c_number += 1;
            }
        }

        // build bonds
        for bond in &mut self.bonds {
            bond.get_bond_type(&self.atoms_idx);
            bond.get_bond_length(&self.coord);
            bond.get_e_str();
            self.etol_str += bond.e_str;
        }

        self.get_neighbor();

        // build angles
        let angle_table = self.get_angle_table();

        self.angle_num = angle_table.len();
        for atms in angle_table {
            let mut angle = Angle::from(atms);
            angle.get_angle(&self.coord);
            angle.get_angle_type(&self.atoms_idx);
            angle.get_e_bend();
            self.etol_bend += angle.e_bend;
            self.angles.push(angle);
        }


        // build dihedrals
        let dihedral_table = self.get_dihedral_table();

        self.dihedral_num = dihedral_table.len();
        for atms in dihedral_table {
            let mut dihedral = Dihedral::from(atms);
            dihedral.get_dihedral(&self.coord);
            dihedral.get_e_tor();
            self.etol_tor += dihedral.e_tor;
            self.dihedrals.push(dihedral);
        }
        
        let atm_pair_table = self.get_atom_pair_table();
        for atm_pair in atm_pair_table {
            let mut bond = Bond::from(atm_pair[0], atm_pair[1]);
            bond.get_bond_type(&self.atoms_idx);
            bond.get_bond_length(&self.coord);
            bond.get_e_vdw();
            self.etol_vdw += bond.e_str;
            self.atm_pair.push(bond);
        }



    }


    pub fn logger(&self) {
        println!("The input file has {} atoms.", self.atm_num);

        println!("\nAtoms and coordinates (in Angstrom):");
        for i in 0..self.atm_num {
            println!("{}  {:-10.6}  {:-10.6}  {:-10.6}", self.atoms[i], self.coord[i][0], self.coord[i][1], self.coord[i][2]);
        }

        println!("\nNumber of coordinates:");
        println!("Stretching: {}  Bending: {}  Torsion: {}", self.bond_num, self.angle_num, self.dihedral_num);
        println!("Internal: {}  Cartesian: {}", self.bond_num + self.angle_num + self.dihedral_num, 3*self.atm_num);
        println!("\nPotential energy at input structure:");
        println!("{:-10.6} kcal/mol", self.etol_str + self.etol_bend + self.etol_tor + self.etol_vdw);  //potential energy
        println!("   Stretch      Bend       Torsion      VDW");
        println!("{:-10.6}  {:-10.6}  {:-10.6}  {:-10.6}", self.etol_str, self.etol_bend, self.etol_tor, self.etol_vdw);
        println!("\nList of all bonds: (At1 - At2, with labels, and distance in Angstrom, energy contrib in kcal/mol)");
        for bond in self.bonds.iter() {
            println!("{} {} - {} {}:    {:-8.5}    {:-8.5}", self.atoms[bond.atms[0]-1], bond.atms[0], self.atoms[bond.atms[1]-1], bond.atms[1], bond.bond_length, bond.e_str);
        }
        println!("\nList of all bending angles: (At1 - At2 - At3, with labels, angle in radian then degrees, energy contribution");
        for angle in self.angles.iter() {
            println!("{} {} - {} {} - {} {}:    {:-8.5}    {:-8.3}    {:-8.5}", self.atoms[angle.atms[0]-1], angle.atms[0], self.atoms[angle.atms[1]-1], angle.atms[1], self.atoms[angle.atms[2]-1], angle.atms[2], angle.angle.to_radians(), angle.angle, angle.e_bend);
        }
        println!("\nList of all torsional angles: (At1 - At2 - At3 - At4, with labels, angle in radian then degrees,  energy contrib in kcal/mol");
        for dihedral in self.dihedrals.iter() {
            println!("{} {} - {} {} - {} {} - {} {}:    {:-8.5}    {:-8.3}    {:-8.5}", self.atoms[dihedral.atms[0]-1], dihedral.atms[0], self.atoms[dihedral.atms[1]-1], dihedral.atms[1], self.atoms[dihedral.atms[2]-1], dihedral.atms[2], self.atoms[dihedral.atms[3]-1], dihedral.atms[3], dihedral.dihedral.to_radians(), dihedral.dihedral, dihedral.e_tor);
        }
        println!("\nList of LJ atom pairs: (At1 - At2, with labels, distance,  vdW energy contrib in kcal/mol");
        for atm_pair in self.atm_pair.iter() {
            println!("{} {} - {} {}:    {:-8.5}    {:-8.5}", self.atoms[atm_pair.atms[0]-1], atm_pair.atms[0], self.atoms[atm_pair.atms[1]-1], atm_pair.atms[1], atm_pair.bond_length, atm_pair.e_str);
        }


    }

    /// construct neighbor table for saturated hydrocarbons
    /// 
    ///         H3   H8
    ///    H4 - C1 - C2 - H7
    ///         H5   H6
    /// Output:
    /// [[2, 3, 4, 5],
    ///  [1, 6, 7, 8],
    ///  [1],
    ///  ..
    ///  [2]]
    /// 
    fn get_neighbor(&mut self) {

        self.neighbors = vec![vec![]; self.atm_num];

        for bond in &self.bonds {
            match bond.bond_type {
                BondType::CC => { 
                    self.neighbors[bond.atms[0]-1].push(bond.atms[1]); self.neighbors[bond.atms[1]-1].push(bond.atms[0]);
                },
                BondType::CH => {
                    if self.atoms_idx[bond.atms[0]-1] == 6 {
                        self.neighbors[bond.atms[0]-1].push(bond.atms[1]);
                        self.neighbors[bond.atms[1]-1].push(bond.atms[0]);
                    } else {
                        self.neighbors[bond.atms[1]-1].push(bond.atms[0]);
                        self.neighbors[bond.atms[0]-1].push(bond.atms[1]);
                    }
                },
                BondType::XX => (),
                BondType::HH => (),
            }
        }

    }


    pub fn get_atom_pair_table(&self) -> Vec<[usize;2]> {
        let mut atm_pair = Vec::new();
        for atm in 1..self.atm_num+1 {
            atm_pair.extend(self.search_LJ_neighbor(atm));
        }
        atm_pair
    }

    /// Search for neighbors of at lease 3 bonds away from an atom
    /// Level order traversal
    pub fn search_LJ_neighbor(&self, atm: usize) -> Vec<[usize;2]> {
        let mut neighbors = Vec::new();
        let mut depth = 1;
        let mut q = VecDeque::new();
        let mut visited = vec![false; self.atm_num];
        visited[atm-1] = true;
        q.push_back(atm.clone());

        while !q.is_empty() {
            let size = q.len();
            for _ in 0..size {
                let node = q.pop_front().unwrap();
                for i in self.neighbors[node-1].iter() {
                    if !visited[*i-1] {
                        q.push_back(*i);
                        visited[*i-1] = true;
                        if depth >= 3 && atm < *i {
                            neighbors.push([atm, *i]);
                        }
                    }
                }
            }

            depth += 1;

        }
        

        neighbors

    }


    /// construct dihedral angle table for saturated hydrocarbons
    ///         H3   H8
    ///    H4 - C1 - C2 - H7
    ///         H5   H6
    /// Output:
    ///[[3, 1, 2, 8], [3, 1, 2, 7], [3, 1, 2, 6], 
    /// [4, 1, 2, 8], [4, 1, 2, 7], [4, 1, 2, 6], 
    /// [5, 1, 2, 8], [5, 1, 2, 7], [5, 1, 2, 6]]
    /// 
    pub fn get_dihedral_table(&self) -> Vec<[usize;4]> {
        let mut table = Vec::new();
        for b in self.bonds.iter() {
            match b.bond_type {
                BondType::CC => {
                    let n_list1 = self.neighbors[b.atms[0]-1].iter().filter(|&&x| x != b.atms[1]).collect::<Vec<&usize>>();
                    let n_list2 = self.neighbors[b.atms[1]-1].iter().filter(|&&x| x != b.atms[0]).collect::<Vec<&usize>>();
                    for n1 in n_list1.iter() {
                        for n2 in n_list2.iter() {
                            table.push([**n1, b.atms[0], b.atms[1], **n2]);
                        }
                    }
                },
                _ => {
                    continue;
                },
            }
        }
        table
    }

    /// construct angle table for saturated hydrocarbons
    ///         H3   H8
    ///    H4 - C1 - C2 - H7
    ///         H5   H6
    /// Output:
    /// [[3, 1, 2], [3, 1, 4], [3, 1, 5], [4, 1, 2], [4, 1, 5], [5, 1, 2],
    ///  [8, 2, 1], [8, 2, 7], [8, 2, 6], [7, 2, 1], [7, 2, 6], [6, 2, 1]]
    /// 
    pub fn get_angle_table(&self) -> Vec<[usize;3]> {
        let mut table = Vec::new();
        for i in 0..self.atm_num {
            if self.atoms_idx[i] == 6 {
                let n_list = &self.neighbors[i];
                for n in 0..n_list.len() {
                    for m in n+1..n_list.len() {
                        table.push([n_list[n], i+1, n_list[m]]);
                    }
                }
            }
        }
        table
    }


}