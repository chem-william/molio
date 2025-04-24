use crate::bond::BondOrder;
use crate::error::CError;
use crate::property::Properties;
use crate::residue::Residue;
use crate::unit_cell::UnitCell;
use crate::{atom::Atom, topology::Topology};
use std::ops::{Index, IndexMut};

#[derive(Debug, Default)]
pub struct Frame {
    pub unit_cell: UnitCell,
    pub properties: Properties,
    topology: Topology,
}

impl Frame {
    pub fn new() -> Self {
        Frame {
            unit_cell: UnitCell::new(),
            properties: Properties::new(),
            topology: Topology::default(),
        }
    }

    /// Get a const reference to the topology of this frame
    ///
    /// It is not possible to get a modifiable reference to the topology,
    /// because it would then be possible to remove/add atoms without changing
    /// the actual positions and velocity storage. Instead, all the mutating
    /// functionalities of the topology are mirrored on the frame (adding and
    /// removing bonds, adding residues, *etc.*)
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn topology_as_mut(&mut self) -> &mut Topology {
        &mut self.topology
    }
    pub fn size(&self) -> usize {
        self.topology.size()
    }

    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.topology
            .atoms
            .iter()
            .map(|a| [a.x, a.y, a.z])
            .collect()
    }

    pub fn add_atom(&mut self, atom: Atom) {
        self.topology.atoms.push(atom)
    }

    pub fn add_residue(&mut self, residue: Residue) -> Result<(), CError> {
        self.topology.add_residue(residue)
    }

    /// Add a bond in the system, between the atoms at index `i` and
    /// `j`.
    pub fn add_bond(&mut self, i: usize, j: usize, bond_order: BondOrder) -> Result<(), CError> {
        self.topology.add_bond(i, j, bond_order)
    }
}

impl Index<usize> for Frame {
    type Output = Atom;

    fn index(&self, index: usize) -> &Self::Output {
        &self.topology.atoms[index]
    }
}

impl IndexMut<usize> for Frame {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.topology.atoms[index]
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    #[test]
    fn test_frame_indexing() {
        let mut frame = Frame::new();
        let atom1 = Atom {
            x: 1.0,
            y: 2.0,
            z: 3.0,
            symbol: "H".to_string(),
            name: "hydrogen".to_string(),
            properties: Properties::new(),
        };
        let atom2 = Atom {
            x: 4.0,
            y: 5.0,
            z: 6.0,
            symbol: "O".to_string(),
            name: "oxygen".to_string(),
            properties: Properties::new(),
        };

        frame.add_atom(atom1);
        frame.add_atom(atom2);

        // Test read access
        assert_eq!(frame[0].symbol, "H");
        assert_eq!(frame[1].symbol, "O");
        assert_approx_eq!(frame[0].x, 1.0);
        assert_approx_eq!(frame[1].x, 4.0);

        // Test write access
        frame[0].x = 10.0;
        assert_approx_eq!(frame[0].x, 10.0);
    }

    #[test]
    #[should_panic]
    fn test_frame_indexing_out_of_bounds() {
        let frame = Frame::new();
        let _ = frame[0];
    }
}
