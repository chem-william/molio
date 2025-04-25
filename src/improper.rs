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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_improper_creation() {
        let improper = Improper::new(5, 1, 8, 3);
        // Check canonical ordering (all atoms except central are sorted)
        assert_eq!(improper[0], 3); // smallest of i,k,m
        assert_eq!(improper[1], 1); // central atom j
        assert_eq!(improper[2], 5); // middle of i,k,m
        assert_eq!(improper[3], 8); // largest of i,k,m
    }

    #[test]
    fn test_improper_equality() {
        // These should be equal due to canonical representation
        let improper1 = Improper::new(5, 1, 8, 3);
        let improper2 = Improper::new(3, 1, 5, 8);
        let improper3 = Improper::new(8, 1, 3, 5);

        assert_eq!(improper1, improper2);
        assert_eq!(improper1, improper3);
        assert_eq!(improper2, improper3);

        // Different improper
        let improper4 = Improper::new(5, 2, 8, 3);
        assert_ne!(improper1, improper4);
    }

    #[test]
    #[should_panic(expected = "cannot have an atom linked to itself in an improper dihedral angle")]
    fn test_improper_with_duplicate_central() {
        Improper::new(1, 2, 2, 3);
    }

    #[test]
    #[should_panic(expected = "cannot have an atom twice in an improper dihedral angle")]
    fn test_improper_with_duplicate_outer() {
        Improper::new(1, 2, 1, 3);
    }

    #[test]
    #[should_panic(expected = "can not access atom n° 4 in dihedral")]
    fn test_improper_index_out_of_bounds() {
        let improper = Improper::new(5, 1, 8, 3);
        let _ = improper[4]; // This should panic
    }
}
