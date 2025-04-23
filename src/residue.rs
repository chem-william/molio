use crate::property::Properties;
use std::{
    cmp::Ordering,
    collections::{BTreeSet, btree_set::Iter},
};

#[derive(Default, PartialEq, Eq, Debug)]
pub struct FullResidueId {
    /// Chain identifier
    pub chain: char,
    /// Residue id
    pub resid: i64,
    /// Residue name
    pub resname: String,
    /// Insertion code of the residue
    pub insertion_code: char,
}

impl Ord for FullResidueId {
    fn cmp(&self, other: &Self) -> Ordering {
        self.chain
            .cmp(&other.chain)
            .then_with(|| self.resid.cmp(&other.resid))
            .then_with(|| self.insertion_code.cmp(&other.insertion_code))
            .then_with(|| self.resname.cmp(&other.resname))
    }
}

impl PartialOrd for FullResidueId {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Default, Debug, Clone)]
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
    pub fn add_atom(&mut self, index: usize) {
        self.atoms.insert(index);
    }

    pub fn size(&self) -> usize {
        self.atoms.len()
    }
}

impl<'a> IntoIterator for &'a Residue {
    type Item = &'a usize;
    type IntoIter = Iter<'a, usize>;

    fn into_iter(self) -> Self::IntoIter {
        self.atoms.iter()
    }
}
