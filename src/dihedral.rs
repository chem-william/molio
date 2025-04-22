use std::ops::Index;

/// Ensures a canonical representation of a dihedral angle
/// between four atoms.
///
/// A dihedral angle is formed by three consecutive bonds:
///
/// ```text
///     |  i       k     |
///     |    \   /   \   |
///     |      j      m  |
/// ```
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Dihedral {
    data: [usize; 4],
}

impl Index<usize> for Dihedral {
    type Output = usize;

    /// Get the index of `index` atom (`index` can be 0, 1, 2 or 3) in the
    /// dihedral.
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

impl Dihedral {
    /// Create a new `Dihedral` containing the atoms `i`, `j`, `k` and `m`.
    pub fn new(i: usize, j: usize, k: usize, m: usize) -> Self {
        if i == j || j == k || k == m {
            panic!("cannot have an atom linked to itself in a dihedral angle");
        }
        if i == k || j == m || i == m {
            panic!("cannot have an atom twice in a dihedral angle");
        }

        let data = if std::cmp::max(i, j) < std::cmp::max(k, m) {
            [i, j, k, m]
        } else {
            [m, k, j, i]
        };

        Dihedral { data }
    }
}
