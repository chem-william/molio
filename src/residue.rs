use crate::property::Properties;
use std::collections::{BTreeSet, btree_set::Iter};

#[derive(Default, PartialEq, Eq, PartialOrd, Ord)]
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

#[derive(Default, Debug)]
pub struct Residue {
    pub name: String,
    pub id: Option<i64>,
    pub atoms: BTreeSet<usize>,
    pub properties: Properties,
}

impl Residue {
    pub fn add_atom(&mut self, index: usize) {
        self.atoms.insert(index);
    }
}
impl<'a> IntoIterator for &'a Residue {
    type Item = &'a usize;
    type IntoIter = Iter<'a, usize>;

    fn into_iter(self) -> Self::IntoIter {
        self.atoms.iter()
    }
}
