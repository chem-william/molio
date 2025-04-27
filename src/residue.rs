use crate::property::{Properties, Property};
use std::collections::{BTreeSet, btree_set::Iter};

#[derive(Default, Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
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
    pub fn new(name: String, id: i64) -> Self {
        Self {
            name,
            id: Some(id),
            ..Default::default()
        }
    }

    pub fn new_from_name(name: String) -> Self {
        Self {
            name,
            ..Default::default()
        }
    }

    pub fn add_atom(&mut self, index: usize) {
        self.atoms.insert(index);
    }

    #[must_use]
    pub fn size(&self) -> usize {
        self.atoms.len()
    }

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
