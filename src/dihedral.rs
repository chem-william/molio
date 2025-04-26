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
#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug)]
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
            panic!("can not access atom n° {index} in dihedral")
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dihedral_creation() {
        let dihedral = Dihedral::new(2, 1, 3, 4);
        // Check canonical ordering
        assert_eq!(dihedral[0], 2);
        assert_eq!(dihedral[1], 1);
        assert_eq!(dihedral[2], 3);
        assert_eq!(dihedral[3], 4);

        // Test with reverse order - should flip the atoms since max(4, 3) > max(1, 2)
        let dihedral_reversed = Dihedral::new(4, 3, 1, 2);
        assert_eq!(dihedral_reversed[0], 2);
        assert_eq!(dihedral_reversed[1], 1);
        assert_eq!(dihedral_reversed[2], 3);
        assert_eq!(dihedral_reversed[3], 4);
    }

    #[test]
    fn test_dihedral_equality() {
        // These should be equal due to canonical representation
        let dihedral1 = Dihedral::new(1, 2, 3, 4);
        let dihedral2 = Dihedral::new(4, 3, 2, 1);

        assert_eq!(dihedral1, dihedral2);

        // Different dihedral
        let dihedral3 = Dihedral::new(1, 2, 3, 5);
        assert_ne!(dihedral1, dihedral3);
    }

    #[test]
    #[should_panic(expected = "cannot have an atom linked to itself in a dihedral angle")]
    fn test_dihedral_with_linked_duplicate_i_j() {
        Dihedral::new(1, 1, 3, 4);
    }

    #[test]
    #[should_panic(expected = "cannot have an atom twice in a dihedral angle")]
    fn test_dihedral_with_distant_duplicate_i_k() {
        Dihedral::new(1, 2, 1, 4);
    }

    #[test]
    #[should_panic(expected = "can not access atom n° 4 in dihedral")]
    fn test_dihedral_index_out_of_bounds() {
        let dihedral = Dihedral::new(1, 2, 3, 4);
        let _ = dihedral[4]; // This should panic
    }
}
