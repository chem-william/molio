use std::ops::Index;

/// Ensures a canonical representation of an improper
/// dihedral angle between four atoms.
///
/// An improper dihedral angle is formed by three bonds around a central atom:
///
/// ```text
///     |  i       k  |
///     |    \   /    |
///     |      j      |
///     |      |      |
///     |      m      |
/// ```
///
///
/// The second atom of the improper is always the central atom.
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy)]
pub struct Improper {
    data: [usize; 4],
}

impl Improper {
    /// Create a new `Improper` containing the atoms `i`, `j`, `k` and `m`. `j`
    /// must be the central atom of the improper.
    ///
    /// # Panics
    ///
    /// Panics if any of `i`, `j`, `k`, `m` has the same value as another
    pub fn new(i: usize, j: usize, k: usize, m: usize) -> Self {
        if j == i || j == k || j == m {
            panic!("cannot have an atom linked to itself in an improper dihedral angle");
        }
        if i == k || i == m || k == m {
            panic!("cannot have an atom twice in an improper dihedral angle");
        }

        // Sort i, k, m
        let mut others = [i, k, m];
        others.sort_unstable();

        Improper {
            data: [others[0], j, others[1], others[2]],
        }
    }
}

impl Index<usize> for Improper {
    type Output = usize;

    /// Get the index of the `i`th atom (`i` can be 0, 1, 2 or 3) in the
    /// improper.
    ///
    /// # Panics
    ///
    /// Panics if `index` is not 0, 1, 2 or 3.
    fn index(&self, index: usize) -> &Self::Output {
        if index >= 4 {
            panic!("can not access atom n° {} in dihedral", index)
        }

        &self.data[index]
    }
}
