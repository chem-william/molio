use std::{collections::BTreeSet, ops::Index};

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

pub struct Bond {
    pub data: [usize; 2],
}

impl PartialEq for Bond {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}
impl Eq for Bond {}

impl PartialOrd for Bond {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.data.partial_cmp(&other.data)
    }
}

impl Ord for Bond {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.data.cmp(&other.data)
    }
}

impl Index<usize> for Bond {
    type Output = usize;

    /// Access one of the two atom indices in the bond.
    ///
    /// # Panics
    ///
    /// Panics if `index >= 2`.
    fn index(&self, index: usize) -> &Self::Output {
        if index >= 2 {
            panic!("can not access atom n° {} in bond", index);
        }
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
        if i == j {
            panic!("can not have a bond between an atom and itself");
        }
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        Bond { data: [a, b] }
    }
}
