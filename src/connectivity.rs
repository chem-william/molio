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
    pub fn angles(&mut self) -> BTreeSet<Angle> {
        if !self.up_to_date {
            self.recalculate()
        }
        self.angles.clone()
    }

    pub fn dihedrals(&mut self) -> BTreeSet<Dihedral> {
        if !self.up_to_date {
            self.recalculate()
        }
        self.dihedrals.clone()
    }

    pub fn impropers(&mut self) -> BTreeSet<Improper> {
        if !self.up_to_date {
            self.recalculate()
        }
        self.impropers.clone()
    }

    pub fn add_bond(&mut self, i: usize, j: usize, bond_order: BondOrder) {
        self.up_to_date = false;

        let bond = Bond::new(i, j);
        let was_new = self.bonds.insert(bond);

        if i > self.biggest_atom {
            self.biggest_atom = i;
        };
        if j > self.biggest_atom {
            self.biggest_atom = j;
        };

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
        let pos = self.bonds.iter().position(|b| *b == bond);

        match pos {
            Some(found_pos) => Ok(self.bond_orders[found_pos]),
            None => Err(CError::GenericError(format!(
                "out of bounds atomic index. No bond between {i} and {j} exists"
            ))),
        }
    }

    fn recalculate(&mut self) {
        self.angles.clear();
        self.dihedrals.clear();
        self.impropers.clear();

        // Generate the list of which atom is bonded to which one
        let mut bonded_to = vec![Vec::new(); self.biggest_atom + 1];
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
            for k in &bonded_to[i] {
                if *k != j {
                    self.angles.insert(Angle::new(*k, i, j));
                }
            }

            for k in &bonded_to[j] {
                if *k != i {
                    self.angles.insert(Angle::new(i, j, *k));
                }
            }
        }

        // Generate list of dihedrals
        for angle in &self.angles {
            let i = angle[0];
            let j = angle[1];
            let k = angle[2];

            for m in &bonded_to[i] {
                if *m != j && *m != k {
                    self.dihedrals.insert(Dihedral::new(*m, i, j, k));
                }
            }

            for m in &bonded_to[k] {
                if *m != i && *m != j {
                    self.dihedrals.insert(Dihedral::new(i, j, k, *m));
                }
            }

            for m in &bonded_to[j] {
                if *m != i && *m != k {
                    self.impropers.insert(Improper::new(i, j, k, *m));
                }
            }
        }

        self.up_to_date = true;
    }
}
