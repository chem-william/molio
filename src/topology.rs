use std::{collections::HashMap, ops::Index};

use crate::{
    angle::Angle,
    atom::Atom,
    bond::{Bond, BondOrder},
    connectivity::Connectivity,
    dihedral::Dihedral,
    error::CError,
    improper::Improper,
    residue::Residue,
};

#[derive(Default, Debug)]
pub struct Topology {
    /// Atoms in the system
    pub atoms: Vec<Atom>,

    /// Connectivity of the system
    connect: Connectivity,

    /// List of residues in the system
    residues: Vec<Residue>,

    /// Association between atom indices and residues indices
    residue_mapping: HashMap<usize, usize>,
}

impl Index<usize> for Topology {
    type Output = Atom;

    fn index(&self, index: usize) -> &Self::Output {
        &self.atoms[index]
    }
}

impl Topology {
    pub fn size(&self) -> usize {
        self.atoms.len()
    }

    pub fn bonds(&self) -> Vec<Bond> {
        self.connect.bonds.iter().cloned().collect()
    }

    pub fn angles(&mut self) -> Vec<Angle> {
        self.connect.angles().iter().cloned().collect()
    }

    pub fn dihedrals(&mut self) -> Vec<Dihedral> {
        self.connect.dihedrals().iter().cloned().collect()
    }

    pub fn impropers(&mut self) -> Vec<Improper> {
        self.connect.impropers().iter().cloned().collect()
    }

    // TODO: should this check be moved to Connectivity::add_bond?
    pub fn add_bond(&mut self, i: usize, j: usize, bond_order: BondOrder) -> Result<(), CError> {
        let amount_atoms = self.size();
        if i >= amount_atoms || j >= amount_atoms {
            return Err(CError::GenericError(format!(
                "out of bounds atomic index. We have {amount_atoms}, but the bond indices are {i} and {j}"
            )));
        };

        self.connect.add_bond(i, j, bond_order);
        Ok(())
    }

    // TODO: should this check be moved to Connectivity::remove_bond?
    pub fn remove_bond(&mut self, i: usize, j: usize) -> Result<(), CError> {
        let amount_atoms = self.size();
        if i >= amount_atoms || j >= amount_atoms {
            return Err(CError::GenericError(format!(
                "out of bounds atomic index. We have {amount_atoms}, but the bond indices are {i} and {j}"
            )));
        };

        self.connect.remove_bond(i, j);
        Ok(())
    }

    pub fn bond_order(&self, i: usize, j: usize) -> Result<BondOrder, CError> {
        self.connect.bond_order(i, j)
    }
}
