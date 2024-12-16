//! Basic structures for molecular geometry and energy calculation
//! 
//! 
//! # Modules:
//! - [`self`] - Molecular geometry structure [`Geom`] and related functions
//! - [`coord`] - Cartesian coordinate structure, including functions for coordinate operations 
//! - [`bond`] - Bond structure, including functions for bond calculation 
//! - [`angle`] - Angle structure, including functions for angle calculation 
//! - [`dihedral`] - Dihedral angle structure, including functions for dihedral angle calculation

#![allow(non_snake_case)]

pub mod coord;
pub mod bond;
pub mod angle;
pub mod dihedral;

use bond::{Bond, BondType};
use angle::Angle;
use dihedral::Dihedral;
use std::fs::File;
use std::io::{BufReader, BufRead};
use coord::CoordVec;
use crate::io::Input;
use crate::matrix::MatFull;
use crate::utils;
use std::collections::VecDeque;
use log::{debug, info, trace};


/// Structure of a molecular geometry

#[derive(Debug, Clone)]
pub struct Geom {
    /// input file
    pub input: Input,  
    /// number of atoms
    pub natm: usize,
    /// atom sequence with element symbols
    pub atoms: Vec<String>,  
    /// atom sequence with atom numbers
    pub atoms_idx: Vec<usize>, 
    /// Cartesian coordinates
    pub coord: Vec<CoordVec>, 
    /// number of bonds
    pub nbond: usize,
    /// bonds
    pub bonds: Vec<Bond>,  
    /// total stretching energy
    pub etot_str: f64, 
    /// number of angles
    pub nangle: usize, 
    /// angles
    pub angles: Vec<Angle>,
    /// total bending energy
    pub etot_bend: f64,  
    /// number of dihedral angles
    pub ndihedral: usize, 
    /// dihedral angles
    pub dihedrals: Vec<Dihedral>, 
    /// total torsion energy
    pub etot_tor: f64, 
    /// number of internal coordinates
    pub nintl: usize,  
    /// number of carbons
    pub c_number: usize, 
    /// neighboring atoms of each atom
    pub neighbors: Vec<Vec<usize>>, 
    /// atoms pairs for vdw calculation
    pub atm_pair: Vec<Bond>, 
    /// total VDW energy
    pub etot_vdw: f64, 
    /// total potential energy
    pub e_tot: f64, 
}

impl Geom {
    /// Create a new molecular geometry
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

    /// Create [`Geom`] structure from .mol2 file
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

    /// Build the [`Geom`] structure, including all bonds, angles, dihedrals, and vdw interactions
    /// as well as their energies. \
    /// **Note**: [`Geom`] only needs to be built once. Use `update()` to update the energy.
    pub fn build(&mut self) {

        self.atoms_idx = utils::elems_to_idx(&self.atoms);

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

    /// Update the energies of [`Geom`] structure \
    /// **Note**: Use this function to update the energy for already built structure
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

    /// Calculate and gather all internal coordinates into a `MatFull<f64>` with size `[nq,1]`
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

    /// construct neighbor table for all atoms
    /// # Example:
    /// ```
    ///         H3   H8
    ///    H4 - C1 - C2 - H7
    ///         H5   H6
    /// ```
    /// Output:
    /// ```
    /// [[2, 3, 4, 5],
    ///  [1, 6, 7, 8],
    ///  [1],
    ///  ...
    ///  [2]]
    /// ```
    fn get_neighbor(&mut self) {

        self.neighbors = vec![vec![]; self.natm];

        for bond in &self.bonds {
            match bond.bond_type {
                BondType::CC => { 
                    self.neighbors[bond.atms[0]-1].push(bond.atms[1]); 
                    self.neighbors[bond.atms[1]-1].push(bond.atms[0]);
                },
                BondType::CH => {
                    self.neighbors[bond.atms[0]-1].push(bond.atms[1]);
                    self.neighbors[bond.atms[1]-1].push(bond.atms[0]);
                },
                BondType::XX => (),
                BondType::HH => (),
            }
        }

    }

    /// Get atom pairs for vdw calculation
    pub fn get_atom_pair_table(&self) -> Vec<([usize;2], bool)> {
        let mut atm_pair = Vec::new();
        for atm in 1..self.natm+1 {
            atm_pair.extend(self.search_LJ_neighbor(atm));
        }
        atm_pair
    }

    /// Search for atom neighbors \
    /// Neighbors at lease 3 bonds away from an atom are considered as LJ pairs 
    // This is a Level order traversal algorithm
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
    /// # Example:
    /// ```
    ///         H3   H8
    ///    H4 - C1 - C2 - H7
    ///         H5   H6
    /// ```
    /// Output:
    /// ```
    ///[[3, 1, 2, 8], [3, 1, 2, 7], [3, 1, 2, 6], 
    /// [4, 1, 2, 8], [4, 1, 2, 7], [4, 1, 2, 6], 
    /// [5, 1, 2, 8], [5, 1, 2, 7], [5, 1, 2, 6]]
    /// ```
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
    /// # Example:
    /// ```
    ///         H3   H8
    ///    H4 - C1 - C2 - H7
    ///         H5   H6
    /// ```
    /// Output:
    /// ```	
    /// [[3, 1, 2], [3, 1, 4], [3, 1, 5], [4, 1, 2], [4, 1, 5], [5, 1, 2],
    ///  [8, 2, 1], [8, 2, 7], [8, 2, 6], [7, 2, 1], [7, 2, 6], [6, 2, 1]]
    /// ```
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

    /// Print out the Cartesian coordinates
    /// # Input:
    /// - `print_level`: print level, `info`, `debug`, `trace`, `off`, or other strings for stdout
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
                    trace!("{}  {:-10.6}  {:-10.6}  {:-10.6}", self.atoms[i], self.coord[i].x, self.coord[i].y, self.coord[i].z);
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


    /// Gather all internal coordinates into a `MatFull<f64>` with size `[nq,1]`
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

    /// Convert `Vec<CoordVec>` to `MatFull<f64>`
    pub fn to_car_coord(&self) -> MatFull<f64> {
        let car_coord = self.coord.iter().map(|x| [x.x, x.y, x.z]).flatten().collect();
        MatFull::from_vec([3*self.natm, 1], car_coord)
    }

}

