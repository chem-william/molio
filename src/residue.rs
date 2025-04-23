use crate::property::{Properties, Property};
use std::collections::{BTreeSet, btree_set::Iter};

#[derive(Default, PartialEq, Eq, Debug, PartialOrd, Ord)]
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

    pub fn contains(&self, index: usize) -> bool {
        self.atoms.contains(&index)
    }

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
