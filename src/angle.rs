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
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Angle {
    pub data: [usize; 3],
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
        if index >= 3 {
            panic!("can not access atom n° {} in angle", index)
        }
        &self.data[index]
    }
}

impl Angle {
    pub fn new(i: usize, j: usize, k: usize) -> Self {
        if i == j || i == k || j == k {
            panic!("can not have the same atom twice in an angle")
        }
        Angle {
            data: [std::cmp::min(i, k), j, std::cmp::max(i, k)],
        }
    }
}
