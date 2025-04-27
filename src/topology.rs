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
    pub residues: Vec<Residue>,

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
    #[must_use]
    pub fn size(&self) -> usize {
        self.atoms.len()
    }

    pub fn resize(&mut self, size: usize) -> Result<(), CError> {
        for bond in &self.connect.bonds {
            if bond[0] >= size || bond[1] >= size {
                return Err(CError::GenericError(format!(
                    "can not resize the topology to contain {size} as there is a bond between atoms {} - {}",
                    bond[0], bond[1]
                )));
            }
        }
        self.atoms.resize(size, Atom::new("X".to_string()));
        Ok(())
    }

    pub fn bonds(&self) -> Vec<Bond> {
        self.connect.bonds.iter().copied().collect()
    }

    pub fn angles(&mut self) -> Vec<Angle> {
        self.connect.angles().iter().copied().collect()
    }

    pub fn dihedrals(&mut self) -> Vec<Dihedral> {
        self.connect.dihedrals().iter().copied().collect()
    }

    pub fn impropers(&mut self) -> Vec<Improper> {
        self.connect.impropers().iter().copied().collect()
    }

    pub fn add_residue(&mut self, residue: Residue) -> Result<(), CError> {
        for &atom_id in &residue {
            if self.residue_mapping.contains_key(&atom_id) {
                return Err(CError::GenericError(format!(
                    "cannot add this residue: atom {atom_id} is already in another residue"
                )));
            }
        }

        let res_index = self.residues.len();
        self.residues.push(residue);

        // Update the mapping
        for &atom_id in &self.residues[res_index] {
            self.residue_mapping.insert(atom_id, res_index);
        }

        Ok(())
    }

    // TODO: should this check be moved to Connectivity::add_bond?
    pub fn add_bond(&mut self, i: usize, j: usize, bond_order: BondOrder) -> Result<(), CError> {
        let amount_atoms = self.size();
        if i >= amount_atoms || j >= amount_atoms {
            return Err(CError::GenericError(format!(
                "out of bounds atomic index. We have {amount_atoms}, but the bond indices are {i} and {j}"
            )));
        }

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
        }

        self.connect.remove_bond(i, j);
        Ok(())
    }

    pub fn bond_order(&self, i: usize, j: usize) -> Result<BondOrder, CError> {
        self.connect.bond_order(i, j)
    }

    #[must_use]
    pub fn residue_for_atom(&self, index: usize) -> Option<Residue> {
        self.residue_mapping
            .get(&index)
            .map(|residue_index| self.residues[*residue_index].clone())
    }

    #[must_use]
    pub fn are_linked(&self, first: &Residue, second: &Residue) -> bool {
        if first == second {
            return true;
        }

        let bonds = self.connect.bonds.clone();
        for &bond_i in first {
            for &bond_j in second {
                let check_bond = Bond::new(bond_i, bond_j);
                if bonds.contains(&check_bond) {
                    return true;
                }
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::{Atom, BondOrder, Residue, Topology};

    #[test]
    fn check_size() {
        let mut topology = Topology::default();
        assert_eq!(topology.size(), 0);

        topology.atoms.push(Atom::new("C".to_string()));
        topology.atoms.push(Atom::new("H".to_string()));
        assert_eq!(topology.size(), 2);
    }

    #[test]
    fn resize_success() {
        let mut topology = Topology::default();
        topology.atoms.push(Atom::new("C".to_string()));

        // Resize to a larger size should succeed
        assert!(topology.resize(3).is_ok());
        assert_eq!(topology.size(), 3);

        // Resize to the same size should succeed
        assert!(topology.resize(3).is_ok());
        assert_eq!(topology.size(), 3);

        // Resize to a smaller size should succeed when no bonds exist
        assert!(topology.resize(2).is_ok());
        assert_eq!(topology.size(), 2);
    }

    #[test]
    fn check_resize_with_bonds() {
        let mut topology = Topology::default();

        // Add three atoms
        topology.atoms.push(Atom::new("C".to_string()));
        topology.atoms.push(Atom::new("H".to_string()));
        topology.atoms.push(Atom::new("O".to_string()));

        // Add a bond
        assert!(topology.add_bond(0, 2, BondOrder::Single).is_ok());

        // Resize to a smaller size that would break the bond should fail
        let result = topology.resize(2);
        assert!(result.is_err());

        // Topology size should remain unchanged
        assert_eq!(topology.size(), 3);
    }

    #[test]
    fn add_bond() {
        let mut topology = Topology::default();

        // Add atoms
        topology.atoms.push(Atom::new("C".to_string()));
        topology.atoms.push(Atom::new("H".to_string()));

        // Adding a valid bond should succeed
        assert!(topology.add_bond(0, 1, BondOrder::Single).is_ok());

        // Adding an out-of-bounds bond should fail
        let result = topology.add_bond(0, 2, BondOrder::Single);
        assert!(result.is_err());
    }

    #[test]
    fn remove_bond() {
        let mut topology = Topology::default();

        // Add atoms
        topology.atoms.push(Atom::new("C".to_string()));
        topology.atoms.push(Atom::new("H".to_string()));

        // Add a bond
        assert!(topology.add_bond(0, 1, BondOrder::Single).is_ok());

        // Removing an existing bond should succeed
        assert!(topology.remove_bond(0, 1).is_ok());

        // Bond should no longer exist
        assert!(topology.bond_order(0, 1).is_err());

        // Removing a non-existent bond should still return Ok
        assert!(topology.remove_bond(0, 1).is_ok());

        // Removing an out-of-bounds bond should fail
        let result = topology.remove_bond(0, 2);
        assert!(result.is_err());
    }

    #[test]
    fn add_residue() {
        let mut topology = Topology::default();

        // Add atoms
        topology.atoms.push(Atom::new("C".to_string()));
        topology.atoms.push(Atom::new("H".to_string()));
        topology.atoms.push(Atom::new("O".to_string()));

        // Create a residue
        let mut residue = Residue::new("ALA".to_string(), 1);
        residue.add_atom(0);
        residue.add_atom(1);

        // Adding a valid residue should succeed
        assert!(topology.add_residue(residue).is_ok());
        assert_eq!(topology.residues.len(), 1);

        // Create another residue with an atom already assigned to the first residue
        let mut residue2 = Residue::new("GLY".to_string(), 2);
        residue2.add_atom(1); // Already in first residue
        residue2.add_atom(2);

        // Adding a residue with overlapping atoms should fail
        let result = topology.add_residue(residue2);
        assert!(result.is_err());
        assert_eq!(topology.residues.len(), 1);
    }

    #[test]
    fn residue_for_atom() {
        let mut topology = Topology::default();

        // Add atoms
        topology.atoms.push(Atom::new("C".to_string()));
        topology.atoms.push(Atom::new("H".to_string()));

        // Create and add a residue
        let mut residue = Residue::new("ALA".to_string(), 1);
        residue.add_atom(0);
        assert!(topology.add_residue(residue).is_ok());

        // Should find residue for atom 0
        let found = topology.residue_for_atom(0);
        assert!(found.is_some());
        assert_eq!(found.unwrap().name, "ALA");

        // Should not find residue for atom 1
        assert!(topology.residue_for_atom(1).is_none());
    }

    #[test]
    fn are_linked() {
        let mut topology = Topology::default();

        // Add atoms
        topology.atoms.push(Atom::new("C".to_string()));
        topology.atoms.push(Atom::new("H".to_string()));
        topology.atoms.push(Atom::new("O".to_string()));
        topology.atoms.push(Atom::new("N".to_string()));

        // Create and add two residues
        let mut residue1 = Residue::new("ALA".to_string(), 1);
        residue1.add_atom(0);
        residue1.add_atom(1);
        assert!(topology.add_residue(residue1.clone()).is_ok());

        let mut residue2 = Residue::new("GLY".to_string(), 2);
        residue2.add_atom(2);
        residue2.add_atom(3);
        assert!(topology.add_residue(residue2.clone()).is_ok());

        // Initially residues are not linked
        assert!(!topology.are_linked(&residue1, &residue2));

        // Add a bond between atoms in different residues
        assert!(topology.add_bond(1, 2, BondOrder::Single).is_ok());

        // Now residues should be linked
        assert!(topology.are_linked(&residue1, &residue2));

        // Same residue is always linked to itself
        assert!(topology.are_linked(&residue1, &residue1));
    }
}
