// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use std::ops::Index;

/// Ensures a canonical representation of an angle between
/// three atoms.
///
/// An angle is formed by two consecutive bonds:
///
/// ```text
///     |  i       k  |
///     |    \   /    |
///     |      j      |
/// ```
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Clone, Copy)]
pub struct Angle {
    data: [usize; 3],
}

impl Index<usize> for Angle {
    type Output = usize;

    /// Get the index of the `i`th atom (`i == 0`, `i == 1` or `i == 2`) in the
    /// angle.
    ///
    /// # Panics
    ///
    /// Panics if `index` is not 0, 1, or 2
    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < 3, "can not access atom n° {index} in angle");
        &self.data[index]
    }
}

impl Angle {
    pub fn new(i: usize, j: usize, k: usize) -> Self {
        assert!(
            !(i == j || i == k || j == k),
            "can not have the same atom twice in an angle"
        );
        Angle {
            data: [std::cmp::min(i, k), j, std::cmp::max(i, k)],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_angle_creation() {
        let angle = Angle::new(5, 1, 8);
        // Check canonical ordering (first and last atoms are sorted)
        assert_eq!(angle[0], 5); // min of i,k
        assert_eq!(angle[1], 1); // central atom j
        assert_eq!(angle[2], 8); // max of i,k

        // Test with reverse order of outer atoms
        let angle_reversed = Angle::new(8, 1, 5);
        assert_eq!(angle_reversed[0], 5); // min of i,k
        assert_eq!(angle_reversed[1], 1); // central atom j
        assert_eq!(angle_reversed[2], 8); // max of i,k
    }

    #[test]
    fn test_angle_equality() {
        // These should be equal due to canonical representation
        let angle1 = Angle::new(5, 1, 8);
        let angle2 = Angle::new(8, 1, 5);

        assert_eq!(angle1, angle2);

        // Different angle
        let angle3 = Angle::new(5, 2, 8);
        assert_ne!(angle1, angle3);
    }

    #[test]
    #[should_panic(expected = "can not have the same atom twice in an angle")]
    fn test_angle_with_duplicate_atoms_i_j() {
        Angle::new(1, 1, 3);
    }

    #[test]
    #[should_panic(expected = "can not access atom n° 3 in angle")]
    fn test_angle_index_out_of_bounds() {
        let angle = Angle::new(5, 1, 8);
        let _ = angle[3]; // This should panic
    }
}
