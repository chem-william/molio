// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::property::{Properties, Property};
use std::collections::{BTreeSet, btree_set::Iter};

#[derive(Default, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
// the specific order of the struct matters in the following when deriving `[Ord]`
// if resname is not compared last, tests will fail
pub struct FullResidueId {
    /// Chain identifier
    pub chain: char,

    /// Residue id
    pub resid: i64,

    /// Insertion code of the residue
    pub insertion_code: char,

    /// Residue name
    // This needs to come after `insertion_code`
    pub resname: String,
}

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Residue {
    /// Name of the residue
    pub name: String,

    /// Index of the residue in the initial topology file
    pub id: Option<i64>,

    /// Indexes of the atoms in this residue. These indexes refers to the
    /// associated topology.
    pub atoms: BTreeSet<usize>,

    /// Additional properties of this residue
    pub properties: Properties,
}

impl Residue {
    #[must_use]
    pub fn new(name: impl Into<String>, id: i64) -> Self {
        Self {
            name: name.into(),
            id: Some(id),
            ..Default::default()
        }
    }

    #[must_use]
    pub fn new_from_name(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }

    /// Add an atom with index `i` to this residue.
    ///
    /// If the atom is already in the residue, this does nothing.
    pub fn add_atom(&mut self, index: usize) {
        self.atoms.insert(index);
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Checks whether [`Self`] is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Check if the residue contains a given atom with index `i`
    #[must_use]
    pub fn contains(&self, index: usize) -> bool {
        self.atoms.contains(&index)
    }

    #[must_use]
    pub fn get(&self, name: &str) -> Option<&Property> {
        self.properties.get(name)
    }
}

impl<'a> IntoIterator for &'a Residue {
    type Item = &'a usize;
    type IntoIter = Iter<'a, usize>;

    fn into_iter(self) -> Self::IntoIter {
        self.atoms.iter()
    }
}

#[allow(dead_code)]
impl Residue {
    fn iter(&self) -> Iter<'_, usize> {
        <&Self as IntoIterator>::into_iter(self)
    }
}
