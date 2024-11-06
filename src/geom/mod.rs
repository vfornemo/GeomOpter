#![allow(non_snake_case)]


pub mod coord;
use std::fs::File;
use std::io::{BufReader, BufRead};
use coord::CoordVec;

use crate::data::{self};
use crate::io::Input;
use crate::matrix::MatFull;
use std::collections::VecDeque;
use log::{debug, info, trace};

#[derive(Debug, Clone)]
pub enum BondType {
    CH,
    CC,
    HH,
    XX,
}

#[derive(Debug, Clone)]
pub struct Bond {
    pub atms: [usize;2],
    pub bond_type: BondType,
    pub bond_length: f64,
    pub e_str: f64,  // stretching energy (or vdw energy for atom pairs)
    pub is_LJ_pair: bool,
}


impl Bond {
    fn new() -> Bond {
        Bond {
            atms: [0usize;2],
            bond_type: BondType::XX,
            bond_length: 0.0f64,
            e_str: 0.0f64,
            is_LJ_pair: false,
        }
    }
    
    fn from(atm1: usize, atm2: usize) -> Bond {
        Bond {
            atms: [atm1, atm2],
            bond_type: BondType::XX,
            bond_length: 0.0f64,
            e_str: 0.0f64,
            is_LJ_pair: false,
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

    fn get_bond_length(&mut self, coord: &Vec<CoordVec>) {
        self.bond_length = coord[self.atms[0]-1].distance(&coord[self.atms[1]-1]);
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
        if self.is_LJ_pair {
            self.e_str = match self.bond_type {
                BondType::CH => {4.0*epsilon_CH*((sigma_CH/self.bond_length).powi(12) - (sigma_CH/self.bond_length).powi(6))},
                BondType::HH => {4.0*data::EPSILON_H*((sigma_HH/self.bond_length).powi(12) - (sigma_HH/self.bond_length).powi(6))},
                BondType::CC => {4.0*data::EPSILON_C*((sigma_CC/self.bond_length).powi(12) - (sigma_CC/self.bond_length).powi(6))},
                BondType::XX => 0.0,
            }
        }

    }

    
}


#[derive(Debug, Clone)]
pub enum AngleType {
    CCC,
    XCX,
    XXX,
}


#[derive(Debug, Clone)]
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

    fn get_angle(&mut self, coord: &Vec<CoordVec>) {
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


#[derive(Debug, Clone)]
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
            atms,
            dihedral: 0.0f64,
            e_tor: 0.0f64,
        }
    }

    fn get_dihedral(&mut self, coord: &Vec<CoordVec>) {
        self.dihedral = coord::get_dihedral_deg(&coord[self.atms[0]-1], &coord[self.atms[1]-1], &coord[self.atms[2]-1], &coord[self.atms[3]-1]);
    }

    fn get_e_tor(&mut self) {
        self.e_tor = data::A_XCCX*(1.0 + (data::N_XCCX*self.dihedral.to_radians()).cos());
    }
    
}



#[derive(Debug, Clone)]
pub struct Geom {
    pub input: Input,  // Input file
    pub natm: usize,  // Number of atoms
    pub atoms: Vec<String>,  // Atom sequence with element symbols
    pub atoms_idx: Vec<usize>,  // Atom sequence with atom numbers
    pub coord: Vec<CoordVec>,  // Cartesian coordinates
    pub nbond: usize,  // Number of bonds
    pub bonds: Vec<Bond>,  // Bonds
    pub etot_str: f64,  // Total Stretching energy
    pub nangle: usize,  // Number of angles
    pub angles: Vec<Angle>,  // Angles
    pub etot_bend: f64,  // Total Bending energy
    pub ndihedral: usize,  // Number of dihedral angles
    pub dihedrals: Vec<Dihedral>,  // Dihedral angles
    pub etot_tor: f64,  // Total Torsion energy
    pub nintl: usize,  // Number of internal coordinates
    pub c_number: usize,  // Number of carbons
    pub neighbors: Vec<Vec<usize>>,  // Neighboring atoms of each atom
    pub atm_pair: Vec<Bond>,  // Atoms pairs for LJ calculation
    pub etot_vdw: f64,  // Total VDW energy
    pub e_tot: f64,  // Total potential energy
}

impl Geom {
    pub fn new(input: Input) -> Geom {
        Geom {
            input,
            natm: 0usize,
            atoms: Vec::new(),
            atoms_idx: Vec::new(),
            coord: Vec::new(),
            nbond: 0usize,
            bonds: Vec::new(),
            etot_str: 0.0f64,
            nangle: 0usize,
            angles: Vec::new(),
            etot_bend: 0.0f64,
            ndihedral: 0usize,
            dihedrals: Vec::new(),
            etot_tor: 0.0f64,
            nintl: 0usize,
            c_number: 0usize,
            neighbors: Vec::new(),
            atm_pair: Vec::new(),
            etot_vdw: 0.0f64,
            e_tot: 0.0f64,
        }
    }

    
    pub fn from_mol2(input: Input) -> Geom {
        let mut geom = Geom::new(input);
        // open file
        let file = File::open(&geom.input.path).expect("File not found");
        let mut reader = BufReader::new(file).lines();
        let line1 = reader.next().unwrap().unwrap();
        
        // first line of mol2 file:
        // number of atoms, number of bonds, number of substructures, number of features, number of sets
        // we only consider the first two numbers
        let mut iter = line1.split_whitespace();
        geom.natm = iter.next().unwrap().parse::<usize>().unwrap();
        geom.nbond = iter.next().unwrap().parse::<usize>().unwrap();
        geom.c_number = iter.next().unwrap().parse::<usize>().unwrap();


        // read atoms
        for _ in 0..geom.natm {
            let line = reader.next().unwrap().unwrap();
            let mut iter = line.split_whitespace();
            let mut coord = CoordVec::new();
            coord.x = iter.next().unwrap().parse::<f64>().unwrap();
            coord.y = iter.next().unwrap().parse::<f64>().unwrap();
            coord.z = iter.next().unwrap().parse::<f64>().unwrap();
            let atom = iter.next().unwrap();
            geom.atoms.push(atom.to_string());
            geom.coord.push(coord);
        }


        // read bonds
        for _ in 0..geom.nbond {
            let line = reader.next().unwrap().unwrap();
            let mut iter = line.split_whitespace();
            let atm1 = iter.next().unwrap().parse::<usize>().unwrap();
            let atm2 = iter.next().unwrap().parse::<usize>().unwrap();
            let bond = Bond::from(atm1, atm2);
            geom.bonds.push(bond);
        }

        geom
    }

    pub fn build(&mut self) {

        self.atoms_idx = data::elems_to_idx(&self.atoms);

        // for i in 0..self.natm {
        //     if self.atoms_idx[i] == 6 {
        //         self.c_number += 1;
        //     }
        // }

        // build bonds
        for bond in &mut self.bonds {
            bond.get_bond_type(&self.atoms_idx);
            bond.get_bond_length(&self.coord);
            bond.get_e_str();
            self.etot_str += bond.e_str;
        }

        self.get_neighbor();

        // build angles
        let angle_table = self.get_angle_table();

        self.nangle = angle_table.len();
        for atms in angle_table {
            let mut angle = Angle::from(atms);
            angle.get_angle(&self.coord);
            angle.get_angle_type(&self.atoms_idx);
            angle.get_e_bend();
            self.etot_bend += angle.e_bend;
            self.angles.push(angle);
        }


        // build dihedrals
        let dihedral_table = self.get_dihedral_table();

        self.ndihedral = dihedral_table.len();
        for atms in dihedral_table {
            let mut dihedral = Dihedral::from(atms);
            dihedral.get_dihedral(&self.coord);
            dihedral.get_e_tor();
            self.etot_tor += dihedral.e_tor;
            self.dihedrals.push(dihedral);
        }
        
        let atm_pair_table = self.get_atom_pair_table();
        for atm_pair in atm_pair_table {
            let mut bond = Bond::from(atm_pair.0[0], atm_pair.0[1]);
            bond.is_LJ_pair = atm_pair.1;
            bond.get_bond_type(&self.atoms_idx);
            bond.get_bond_length(&self.coord);
            bond.get_e_vdw();
            self.etot_vdw += bond.e_str;
            self.atm_pair.push(bond);
        }

        self.e_tot = self.etot_str + self.etot_bend + self.etot_tor + self.etot_vdw;
        self.nintl = self.nbond + self.nangle + self.ndihedral;


    }

    pub fn update(&mut self) {
        self.etot_str = 0.0;
        self.etot_bend = 0.0;
        self.etot_tor = 0.0;
        self.etot_vdw = 0.0;
        self.e_tot = 0.0;

        // build bonds
        for b in self.bonds.iter_mut() {
            b.get_bond_length(&self.coord);
            b.get_e_str();
            self.etot_str += b.e_str;
        }

        for a in self.angles.iter_mut() {
            a.get_angle(&self.coord);
            a.get_e_bend();
            self.etot_bend += a.e_bend;
        }

        for d in self.dihedrals.iter_mut() {
            d.get_dihedral(&self.coord);
            d.get_e_tor();
            self.etot_tor += d.e_tor;
        }
        
        for ap in self.atm_pair.iter_mut() {
            ap.get_bond_length(&self.coord);
            ap.get_e_vdw();
            self.etot_vdw += ap.e_str;
        }

        self.e_tot = self.etot_str + self.etot_bend + self.etot_tor + self.etot_vdw;

    }

    pub fn calc_intl_coord(&mut self) -> MatFull<f64> {

        // build bonds
        for b in self.bonds.iter_mut() {
            b.get_bond_length(&self.coord);
        }

        for a in self.angles.iter_mut() {
            a.get_angle(&self.coord);
        }

        for d in self.dihedrals.iter_mut() {
            d.get_dihedral(&self.coord);
        }

        self.to_intl_coord()

    }


    /// Print out the information of the molecular geometry
    /// Input:
    /// level: usize, 0 for basic info, 1 for detailed info
    pub fn logger(&self) {
   
        info!("The input file has {} atoms.", self.natm);
        info!("\nAtoms and coordinates (in Angstrom):");
        self.formated_output_coord("info");
        info!("\nNumber of coordinates:");
        info!("Stretching: {}  Bending: {}  Torsion: {}", self.nbond, self.nangle, self.ndihedral);
        info!("Internal: {}  Cartesian: {}", self.nintl, 3*self.natm);
        info!("\nPotential energy at input structure:");
        info!("{:-10.6} kcal/mol", self.e_tot);  //potential energy
        info!("Stretch, Bend, Torsion, VDW components of potential energy:");
        info!("{:-10.6}  {:-10.6}  {:-10.6}  {:-10.6}", self.etot_str, self.etot_bend, self.etot_tor, self.etot_vdw); 

        trace!("\nList of all bonds: (At1 - At2, with labels, and distance in Angstrom, energy contrib in kcal/mol)");
        for bond in self.bonds.iter() {
            trace!("{} {} - {} {}:    {:-8.5}    {:-8.5}", self.atoms[bond.atms[0]-1], bond.atms[0], self.atoms[bond.atms[1]-1], bond.atms[1], bond.bond_length, bond.e_str);
        }
        trace!("\nList of all bending angles: (At1 - At2 - At3, with labels, angle in radian then degrees, energy contribution");
        for angle in self.angles.iter() {
            trace!("{} {} - {} {} - {} {}:    {:-8.5}    {:-8.3}    {:-8.5}", self.atoms[angle.atms[0]-1], angle.atms[0], self.atoms[angle.atms[1]-1], angle.atms[1], self.atoms[angle.atms[2]-1], angle.atms[2], angle.angle.to_radians(), angle.angle, angle.e_bend);
        }
        trace!("\nList of all torsional angles: (At1 - At2 - At3 - At4, with labels, angle in radian then degrees,  energy contrib in kcal/mol");
        for dihedral in self.dihedrals.iter() {
            trace!("{} {} - {} {} - {} {} - {} {}:    {:-8.5}    {:-8.3}    {:-8.5}", self.atoms[dihedral.atms[0]-1], dihedral.atms[0], self.atoms[dihedral.atms[1]-1], dihedral.atms[1], self.atoms[dihedral.atms[2]-1], dihedral.atms[2], self.atoms[dihedral.atms[3]-1], dihedral.atms[3], dihedral.dihedral.to_radians(), dihedral.dihedral, dihedral.e_tor);
        }
        trace!("\nList of LJ atom pairs: (At1 - At2, with labels, distance,  vdW energy contrib in kcal/mol");
        for atm_pair in self.atm_pair.iter() {
            trace!("{} {} - {} {}:    {:-8.5}    {:-8.5}", self.atoms[atm_pair.atms[0]-1], atm_pair.atms[0], self.atoms[atm_pair.atms[1]-1], atm_pair.atms[1], atm_pair.bond_length, atm_pair.e_str);
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

        self.neighbors = vec![vec![]; self.natm];

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


    pub fn get_atom_pair_table(&self) -> Vec<([usize;2], bool)> {
        let mut atm_pair = Vec::new();
        for atm in 1..self.natm+1 {
            atm_pair.extend(self.search_LJ_neighbor(atm));
        }
        atm_pair
    }

    /// Search for neighbors of at lease 3 bonds away from an atom
    /// Level order traversal
    pub fn search_LJ_neighbor(&self, atm: usize) -> Vec<([usize;2], bool)> {
        let mut neighbors = Vec::new();
        let mut depth = 1;
        let mut q = VecDeque::new();
        let mut visited = vec![false; self.natm];
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
                            neighbors.push(([atm, *i], true));
                        } else if depth < 3 && atm < *i {
                            neighbors.push(([atm, *i], false));
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
        for i in 0..self.natm {
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

    pub fn formated_output_coord(&self, print_level: &str) {

        match print_level {
            "info" => {
                for i in 0..self.natm {
                    info!("{}  {:-10.6}  {:-10.6}  {:-10.6}", self.atoms[i], self.coord[i].x, self.coord[i].y, self.coord[i].z);
                }
            },
            "debug" => {
                for i in 0..self.natm {
                    debug!("{}  {:-10.6}  {:-10.6}  {:-10.6}", self.atoms[i], self.coord[i].x, self.coord[i].y, self.coord[i].z);
                }
            },
            "trace" => {
                for i in 0..self.natm {
                    println!("{}  {:-10.6}  {:-10.6}  {:-10.6}", self.atoms[i], self.coord[i].x, self.coord[i].y, self.coord[i].z);
                }
            },
            "off" => (),
            _ => {
                for i in 0..self.natm {
                    println!("{}  {:-10.6}  {:-10.6}  {:-10.6}", self.atoms[i], self.coord[i].x, self.coord[i].y, self.coord[i].z);
                }
            },
        }
    }



    pub fn to_intl_coord(&self) -> MatFull<f64> {
        let mut intl_coord = Vec::new();
        for b in self.bonds.iter() {
            intl_coord.push(b.bond_length);
        }
        for a in self.angles.iter() {
            intl_coord.push(a.angle.to_radians());
        }
        for d in self.dihedrals.iter() {
            intl_coord.push(d.dihedral.to_radians());
        }
        MatFull::from_vec([self.nintl, 1], intl_coord)
    }

    /// Convert Vec\<CoordVec\> to MatFull\<f64\>
    pub fn to_car_coord(&self) -> MatFull<f64> {
        let car_coord = self.coord.iter().map(|x| [x.x, x.y, x.z]).flatten().collect();
        MatFull::from_vec([3*self.natm, 1], car_coord)
    }

}

