// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use smallvec::SmallVec;
use std::collections::BTreeSet;

use crate::{
    angle::Angle,
    bond::{Bond, BondOrder},
    dihedral::Dihedral,
    error::CError,
    improper::Improper,
};

#[derive(Default, Debug)]
pub struct Connectivity {
    /// Bonds in this connectivity
    pub bonds: BTreeSet<Bond>,

    /// Angles in this connectivity
    pub angles: BTreeSet<Angle>,

    /// Dihedrals in this connectivity
    pub dihedrals: BTreeSet<Dihedral>,

    /// Impropers in this connectivity
    pub impropers: BTreeSet<Improper>,

    /// Bond orders in this connectivity
    pub bond_orders: Vec<BondOrder>,

    /// Is the cached content up to date?
    up_to_date: bool,

    /// Biggest index within the atoms we know about. Used to pre-allocate
    /// memory when recomputing bonds.
    biggest_atom: usize,
}

/// Add a bond between atoms `i` and `j`
/// Remove any bond between atoms `i` and `j`
/// Update the indices of the bonds after atom removal
///
/// This shifts all indices bigger than `index` in the
/// bonds/angles/dihedrals/impropers lists by -1
/// Get the bond order of the bond between `i` and `j`
impl Connectivity {
    pub fn angles(&mut self) -> &BTreeSet<Angle> {
        if !self.up_to_date {
            self.recalculate();
        }
        &self.angles
    }

    pub fn dihedrals(&mut self) -> &BTreeSet<Dihedral> {
        if !self.up_to_date {
            self.recalculate();
        }
        &self.dihedrals
    }

    pub fn impropers(&mut self) -> &BTreeSet<Improper> {
        if !self.up_to_date {
            self.recalculate();
        }
        &self.impropers
    }

    pub fn add_bond(&mut self, i: usize, j: usize, bond_order: BondOrder) {
        self.up_to_date = false;

        let bond = Bond::new(i, j);
        let was_new = self.bonds.insert(bond);

        if i > self.biggest_atom {
            self.biggest_atom = i;
        }

        if j > self.biggest_atom {
            self.biggest_atom = j;
        }

        if was_new {
            let diff = self
                .bonds
                .iter()
                .position(|b| *b == bond)
                .expect("we just inserted the element");
            self.bond_orders.insert(diff, bond_order);
        }
    }

    pub fn remove_bond(&mut self, i: usize, j: usize) {
        let bond = Bond::new(i, j);
        let pos = self.bonds.iter().position(|b| *b == bond);

        if let Some(found_pos) = pos {
            self.up_to_date = false;
            self.bonds.remove(&bond);

            self.bond_orders.remove(found_pos);

            debug_assert_eq!(self.bond_orders.len(), self.bonds.len());
        }
    }

    pub fn bond_order(&self, i: usize, j: usize) -> Result<BondOrder, CError> {
        let bond = Bond::new(i, j);
        self.bonds
            .iter()
            .position(|b| *b == bond)
            .map(|pos| self.bond_orders[pos])
            .ok_or_else(|| {
                CError::GenericError(format!(
                    "out of bounds atomic index. No bond between {i} and {j} exists"
                ))
            })
    }

    fn recalculate(&mut self) {
        self.angles.clear();
        self.dihedrals.clear();
        self.impropers.clear();

        // Pre-allocate space for bonded_to vectors
        let mut bonded_to = vec![SmallVec::<[usize; 4]>::new(); self.biggest_atom + 1];

        // Generate the list of which atom is bonded to which one
        for bond in &self.bonds {
            debug_assert!(bond[0] < bonded_to.len());
            debug_assert!(bond[1] < bonded_to.len());
            bonded_to[bond[0]].push(bond[1]);
            bonded_to[bond[1]].push(bond[0]);
        }

        // Generate list of angles
        for bond in &self.bonds {
            let i = bond[0];
            let j = bond[1];

            for &k in &bonded_to[i] {
                if k != j {
                    self.angles.insert(Angle::new(k, i, j));
                }
            }

            for &k in &bonded_to[j] {
                if k != i {
                    self.angles.insert(Angle::new(i, j, k));
                }
            }
        }

        // Generate list of dihedrals and impropers
        for angle in &self.angles {
            let i = angle[0];
            let j = angle[1];
            let k = angle[2];

            for &m in &bonded_to[i] {
                if m != j && m != k {
                    self.dihedrals.insert(Dihedral::new(m, i, j, k));
                }
            }

            for &m in &bonded_to[k] {
                if m != i && m != j {
                    self.dihedrals.insert(Dihedral::new(i, j, k, m));
                }
            }

            for &m in &bonded_to[j] {
                if m != i && m != k {
                    self.impropers.insert(Improper::new(i, j, k, m));
                }
            }
        }

        self.up_to_date = true;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_bond() {
        let mut connectivity = Connectivity::default();
        // Add a bond
        connectivity.add_bond(1, 2, BondOrder::Single);

        // Check that the bond was added
        assert_eq!(connectivity.bonds.len(), 1);
        assert!(connectivity.bonds.contains(&Bond::new(1, 2)));

        // Check bond order
        assert_eq!(connectivity.bond_order(1, 2).unwrap(), BondOrder::Single);

        // Adding the same bond again should not increase count
        connectivity.add_bond(1, 2, BondOrder::Double);
        assert_eq!(connectivity.bonds.len(), 1);

        // Check that the bond order was not updated
        assert_eq!(connectivity.bond_order(1, 2).unwrap(), BondOrder::Single);

        // Adding the same bond in reverse order should not increase count
        connectivity.add_bond(2, 1, BondOrder::Single);
        assert_eq!(connectivity.bonds.len(), 1);

        // Check the bond order was updated
        assert_eq!(connectivity.bond_order(1, 2).unwrap(), BondOrder::Single);
    }

    #[test]
    fn test_remove_bond() {
        let mut connectivity = Connectivity::default();
        connectivity.add_bond(1, 2, BondOrder::Single);
        connectivity.add_bond(2, 3, BondOrder::Double);
        connectivity.add_bond(3, 4, BondOrder::Triple);

        assert_eq!(connectivity.bonds.len(), 3);

        // Remove a bond
        connectivity.remove_bond(1, 2);
        assert_eq!(connectivity.bonds.len(), 2);
        assert!(!connectivity.bonds.contains(&Bond::new(1, 2)));

        // Removing a non-existent bond should do nothing
        connectivity.remove_bond(1, 2);
        assert_eq!(connectivity.bonds.len(), 2);

        // Check that bond orders array has same length as bonds set
        assert_eq!(connectivity.bond_orders.len(), connectivity.bonds.len());
    }

    #[test]
    fn test_bond_order() {
        let mut connectivity = Connectivity::default();
        connectivity.add_bond(1, 2, BondOrder::Single);
        connectivity.add_bond(2, 3, BondOrder::Double);
        connectivity.add_bond(3, 4, BondOrder::Triple);

        // Check bond orders
        assert_eq!(connectivity.bond_order(1, 2).unwrap(), BondOrder::Single);
        assert_eq!(connectivity.bond_order(2, 3).unwrap(), BondOrder::Double);
        assert_eq!(connectivity.bond_order(3, 4).unwrap(), BondOrder::Triple);

        // Check bond orders with reversed indices
        assert_eq!(connectivity.bond_order(2, 1).unwrap(), BondOrder::Single);
        assert_eq!(connectivity.bond_order(3, 2).unwrap(), BondOrder::Double);
        assert_eq!(connectivity.bond_order(4, 3).unwrap(), BondOrder::Triple);

        // Bond order for non-existent bond should return error
        assert!(connectivity.bond_order(1, 3).is_err());
    }

    #[test]
    fn test_angles_generation() {
        let mut connectivity = Connectivity::default();

        // Create a simple chain: 1-2-3
        connectivity.add_bond(1, 2, BondOrder::Single);
        connectivity.add_bond(2, 3, BondOrder::Single);

        let angles = connectivity.angles();
        assert_eq!(angles.len(), 1);
        assert!(angles.contains(&Angle::new(1, 2, 3)));

        // Add another atom to create a branch: 1-2-3
        //                                         |
        //                                         4
        connectivity.add_bond(2, 4, BondOrder::Single);

        // Force recalculation by getting angles again
        let angles = connectivity.angles();
        assert_eq!(angles.len(), 3);
        assert!(angles.contains(&Angle::new(1, 2, 3)));
        assert!(angles.contains(&Angle::new(1, 2, 4)));
        assert!(angles.contains(&Angle::new(3, 2, 4)));
    }

    #[test]
    fn test_dihedrals_generation() {
        let mut connectivity = Connectivity::default();

        // Create a simple chain: 1-2-3-4
        connectivity.add_bond(1, 2, BondOrder::Single);
        connectivity.add_bond(2, 3, BondOrder::Single);
        connectivity.add_bond(3, 4, BondOrder::Single);

        let dihedrals = connectivity.dihedrals();
        assert_eq!(dihedrals.len(), 1);
        assert!(dihedrals.contains(&Dihedral::new(1, 2, 3, 4)));

        // Add another atom to create a branch: 1-2-3-4
        //                                          |
        //                                          5
        connectivity.add_bond(3, 5, BondOrder::Single);

        // Force recalculation
        let dihedrals = connectivity.dihedrals().clone();
        let impropers = connectivity.impropers().clone();
        assert_eq!(dihedrals.len(), 2);
        assert_eq!(impropers.len(), 1);
        assert!(dihedrals.contains(&Dihedral::new(1, 2, 3, 4)));
        assert!(dihedrals.contains(&Dihedral::new(1, 2, 3, 5)));
        assert!(impropers.contains(&Improper::new(2, 3, 4, 5)));
    }

    #[test]
    fn test_impropers_generation() {
        let mut connectivity = Connectivity::default();

        // Create a structure with potential improper:
        //       1
        //       |
        //   4---2---3
        connectivity.add_bond(1, 2, BondOrder::Single);
        connectivity.add_bond(2, 3, BondOrder::Single);
        connectivity.add_bond(2, 4, BondOrder::Single);

        let impropers = connectivity.impropers();
        assert_eq!(impropers.len(), 1);
        assert!(impropers.contains(&Improper::new(1, 2, 3, 4)));

        // Add another branch
        //       1
        //       |
        //   4---2---3
        //       |
        //       5
        connectivity.add_bond(2, 5, BondOrder::Single);

        // Force recalculation
        let impropers = connectivity.impropers();
        assert_eq!(impropers.len(), 4);
        // All possible combinations with 2 as the central atom
        assert!(impropers.contains(&Improper::new(1, 2, 3, 4)));
        assert!(impropers.contains(&Improper::new(1, 2, 3, 5)));
        assert!(impropers.contains(&Improper::new(1, 2, 4, 5)));
        assert!(impropers.contains(&Improper::new(3, 2, 4, 5)));
    }
}
