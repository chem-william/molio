// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use std::ops::Index;

/// Stores the type of bond
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondOrder {
    /// Bond order is unknown or unspecified
    Unknown,

    /// Single bond
    Single,

    /// Double bond
    Double,

    /// Triple bond
    Triple,

    /// Quadruplet bond
    Quadruple,

    /// Quintuplet bond
    Quintuplet,

    /// Single bond direction from first atom to second is 'down'. Used for cis-trans isomers
    Down,

    /// Single bond direction from first atom to second is 'up'. Used for cis-trans isomers
    Up,

    /// Dative bond where the electrons are localized to the first atom
    DativeR,

    /// Dative bond where the electrons are localized to the second atom
    DativeL,

    /// Amide bond (C(=O)-NH)
    Amide,

    /// Aromatic bond (for example the ring bonds in benzene)
    Aromatic,
}

/// Ensure a canonical representation of a bond two atoms
#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug)]
pub struct Bond {
    data: [usize; 2],
}

impl Index<usize> for Bond {
    type Output = usize;

    /// Access one of the two atom indices in the bond.
    ///
    /// # Panics
    ///
    /// Panics if `index >= 2`.
    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < 2, "cannot access atom n° {index} in bond");
        &self.data[index]
    }
}

impl Bond {
    /// Create a new bond between `i` and `j`.
    ///
    /// # Panics
    ///
    /// Panics if `i == j`.
    pub fn new(i: usize, j: usize) -> Self {
        assert!(i != j, "can not have a bond between an atom and itself");
        let a = std::cmp::min(i, j);
        let b = std::cmp::max(i, j);
        Bond { data: [a, b] }
    }
}
