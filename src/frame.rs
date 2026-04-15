// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::bond::BondOrder;
use crate::error::CError;
use crate::property::Properties;
use crate::residue::Residue;
use crate::unit_cell::UnitCell;
use crate::{atom::Atom, topology::Topology};
use std::ops::{Deref, Index, IndexMut};

#[derive(Debug, Default)]
pub struct Frame {
    /// Unit cell of the system
    pub unit_cell: UnitCell,

    /// Properties stored in this frame
    pub properties: Properties,

    /// Positions of the particles
    positions: Vec<[f64; 3]>,

    /// Velocities of the particles.
    pub(crate) velocities: Option<Vec<[f64; 3]>>,

    /// Topology of the described system
    topology: Topology,
}

impl Frame {
    #[must_use]
    pub fn new() -> Self {
        Frame {
            unit_cell: UnitCell::new(),
            properties: Properties::new(),
            positions: vec![],
            velocities: None,
            topology: Topology::default(),
        }
    }

    #[must_use]
    pub fn from_unitcell(unit_cell: UnitCell) -> Self {
        Frame {
            unit_cell,
            properties: Properties::new(),
            positions: vec![],
            velocities: None,
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
    #[must_use]
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn set_topology(&mut self, topology: Topology) -> Result<(), CError> {
        if topology.size() != self.size() {
            return Err(CError::GenericError(format!(
                "the topology contains {} atoms, but the frame contains {} atoms",
                topology.size(),
                self.size(),
            )));
        }

        self.topology = topology;

        Ok(())
    }

    pub fn topology_mut(&mut self) -> &mut Topology {
        &mut self.topology
    }

    /// Get a reference to the unit cell of this [`Frame`]
    #[must_use]
    pub fn cell(&self) -> &UnitCell {
        &self.unit_cell
    }

    /// Set the unit cell of this frame to `cell`
    pub fn set_unit_cell(&mut self, cell: UnitCell) {
        self.unit_cell = cell;
    }

    #[must_use]
    pub fn positions(&self) -> &Vec<[f64; 3]> {
        &self.positions
    }

    pub fn positions_mut(&mut self) -> &mut Vec<[f64; 3]> {
        &mut self.positions
    }

    #[must_use]
    pub fn velocities(&self) -> Option<&Vec<[f64; 3]>> {
        self.velocities.as_ref()
    }

    pub fn velocities_mut(&mut self) -> Option<&mut Vec<[f64; 3]>> {
        self.velocities.as_mut()
    }

    pub fn add_atom(&mut self, atom: Atom, position: [f64; 3]) {
        self.topology.atoms.push(atom);
        self.positions.push(position);
    }

    pub fn add_atom_with_velocity(&mut self, atom: Atom, position: [f64; 3], velocity: [f64; 3]) {
        self.topology.atoms.push(atom);
        self.positions.push(position);
        if let Some(velocities) = self.velocities.as_mut() {
            velocities.push(velocity);
        }
    }

    /// Remove all connectivity information in the frame's topology
    pub fn clear_bonds(&mut self) {
        self.topology.clear_bonds();
    }

    pub fn add_residue(&mut self, residue: Residue) -> Result<(), CError> {
        self.topology.add_residue(residue)
    }

    pub fn resize(&mut self, size: usize) -> Result<(), CError> {
        self.topology.resize(size)?;
        self.positions.resize(size, [0.0; 3]);
        if let Some(velocities) = self.velocities.as_mut() {
            velocities.resize(size, [0.0; 3]);
        }
        Ok(())
    }

    /// Reserves capacity for at least `size` more elements in this frames
    /// [`Topology`] and [`Self::positions`].
    ///
    /// # Panics
    /// Panics if the new capacity of the underlying `Vec`s exceed `isize::MAX` bytes
    pub fn reserve(&mut self, size: usize) {
        self.topology.reserve(size);
        self.positions.reserve(size);
        // if self.velocities.is_some() {
        //     self.velocities.reserve();
        // }
    }

    /// Add a bond in the system, between the atoms at index `i` and
    /// `j`.
    pub fn add_bond(&mut self, i: usize, j: usize, bond_order: BondOrder) -> Result<(), CError> {
        self.topology.add_bond(i, j, bond_order)
    }

    pub(crate) fn add_velocities(&mut self) {
        if self.velocities.is_none() {
            self.velocities = Some(vec![[0.0; 3]; self.size()]);
        }
    }

    pub(crate) fn size(&self) -> usize {
        debug_assert!(self.positions.len() == self.topology.size());

        if let Some(velocities) = self.velocities.as_ref() {
            debug_assert!(self.positions.len() == velocities.len());
        }

        self.positions.len()
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

impl Deref for Frame {
    type Target = Vec<Atom>;

    fn deref(&self) -> &Self::Target {
        &self.topology.atoms
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
            symbol: "H".to_string(),
            name: "hydrogen".to_string(),
            charge: 0.0,
            mass: 0.0,
            properties: Properties::new(),
        };
        let atom2 = Atom {
            symbol: "O".to_string(),
            name: "oxygen".to_string(),
            charge: 0.0,
            mass: 0.0,
            properties: Properties::new(),
        };

        frame.add_atom(atom1, [1.0, 2.0, 3.0]);
        frame.add_atom(atom2, [4.0, 5.0, 6.0]);

        // Test read access
        assert_eq!(frame[0].symbol, "H");
        assert_eq!(frame[1].symbol, "O");
        let pos = frame.positions();
        assert_approx_eq!(pos[0][0], 1.0);
        assert_approx_eq!(pos[1][0], 4.0);
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn test_frame_indexing_out_of_bounds() {
        let frame = Frame::new();
        let _ = frame[0];
    }
}
