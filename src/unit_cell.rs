use crate::error::CError;
use core::f64;
use nalgebra::Matrix3;

type Vec3D = [f64; 3];

#[derive(Debug, PartialEq, Eq, Default)]
pub enum CellShape {
    #[default]
    Orthorhombic,
    Triclinic,
    Infinite,
}

mod utils {
    use super::Vec3D;
    use nalgebra::Matrix3;

    /// Check if a value is approximately zero (within 1e-5)
    pub fn is_roughly_zero(x: f64) -> bool {
        x.abs() < 1e-5
    }

    /// Check if an angle is approximately 90 degrees (within 1e-3)
    pub fn is_roughly_90(value: f64) -> bool {
        (value - 90.0).abs() < 1e-3
    }

    /// Calculate the lengths of the cell vectors from the matrix
    pub fn calc_lengths_from_matrix(matrix: Matrix3<f64>) -> Vec3D {
        let v1 = matrix.row(0);
        let v2 = matrix.row(1);
        let v3 = matrix.row(2);

        [v1.norm(), v2.norm(), v3.norm()]
    }

    /// Calculate the angles between cell vectors from the matrix
    pub fn calc_angles_from_matrix(matrix: Matrix3<f64>) -> Vec3D {
        let v1 = matrix.row(0);
        let v2 = matrix.row(1);
        let v3 = matrix.row(2);

        [
            (v2.dot(&v3) / (v2.norm() * v3.norm())).acos().to_degrees(),
            (v1.dot(&v3) / (v1.norm() * v3.norm())).acos().to_degrees(),
            (v1.dot(&v2) / (v1.norm() * v2.norm())).acos().to_degrees(),
        ]
    }

    /// Calculate cosine of an angle in degrees
    pub fn cos_degree(theta: f64) -> f64 {
        theta.to_radians().cos()
    }

    /// Calculate sine of an angle in degrees
    pub fn sin_degree(theta: f64) -> f64 {
        theta.to_radians().sin()
    }

    /// Check if a cell is orthorhombic based on its lengths and angles
    pub fn is_orthorhombic(lengths: Vec3D, angles: Vec3D) -> bool {
        if is_infinite(lengths) {
            return false;
        }

        // we support cells with one or two lengths of 0 which results in NaN angles
        angles
            .iter()
            .all(|&angle| is_roughly_90(angle) || angle.is_nan())
    }

    /// Check if a cell is infinite (all lengths are zero)
    pub fn is_infinite(lengths: Vec3D) -> bool {
        lengths.iter().all(|&x| is_roughly_zero(x))
    }

    /// Check if a matrix is diagonal (all off-diagonal elements are zero)
    pub fn is_diagonal(matrix: Matrix3<f64>) -> bool {
        matrix
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 4 != 0) // Skip diagonal elements
            .all(|(_, &x)| is_roughly_zero(x))
    }
}

mod validation {
    use super::{utils, CError, Vec3D};

    /// Validate cell angles
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Any length is negative
    pub fn check_lengths(lengths: &Vec3D) -> Result<(), CError> {
        if let Some(&length) = lengths.iter().find(|&&x| x < 0.0) {
            return Err(CError::GenericError(format!(
                "negative length found: {length}"
            )));
        }

        Ok(())
    }

    /// Validate cell angles
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Any angle is negative
    /// - Any angle is zero
    /// - Any angle is 180 degrees or greater
    pub fn check_angles(angles: &Vec3D) -> Result<(), CError> {
        if let Some(&angle) = angles.iter().find(|&&x| x < 0.0) {
            return Err(CError::GenericError(format!(
                "negative angle found: {angle}"
            )));
        }

        if let Some(&angle) = angles.iter().find(|&&x| utils::is_roughly_zero(x)) {
            return Err(CError::GenericError(format!("zero angle found: {angle}")));
        }

        if let Some(&angle) = angles.iter().find(|&&x| x >= 180.0) {
            return Err(CError::GenericError(format!("angle too large: {angle}")));
        }

        Ok(())
    }
}

mod matrix {
    use super::{utils, validation, CError, Matrix3, Vec3D};

    /// Create a cell matrix from lengths and angles
    ///
    /// # Errors
    ///
    /// Returns an error if the lengths or angles are invalid
    pub fn cell_matrix_from_lengths_angles(
        lengths: Vec3D,
        angles: &mut Vec3D,
    ) -> Result<Matrix3<f64>, CError> {
        validation::check_lengths(&lengths)?;
        validation::check_angles(angles)?;

        // Normalize angles to 90 degrees if they're close enough
        if angles.iter().all(|&x| utils::is_roughly_90(x)) {
            angles.iter_mut().for_each(|x| *x = 90.0);
        }

        let mut cell_matrix = Matrix3::zeros();
        cell_matrix[(0, 0)] = lengths[0];

        let cos_gamma = utils::cos_degree(angles[2]);
        let sin_gamma = utils::sin_degree(angles[2]);
        cell_matrix[(1, 0)] = cos_gamma * lengths[1];
        cell_matrix[(1, 1)] = sin_gamma * lengths[1];

        let cos_beta = utils::cos_degree(angles[1]);
        let cos_alpha = utils::cos_degree(angles[0]);

        cell_matrix[(2, 0)] = cos_beta;
        cell_matrix[(2, 1)] = (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
        cell_matrix[(2, 2)] =
            (1.0 - cell_matrix[(2, 0)].powi(2) - cell_matrix[(2, 1)].powi(2)).sqrt();
        cell_matrix[(2, 0)] *= lengths[2];
        cell_matrix[(2, 1)] *= lengths[2];
        cell_matrix[(2, 2)] *= lengths[2];

        Ok(cell_matrix)
    }
}

#[derive(Debug, Default)]
pub struct UnitCell {
    pub matrix: Matrix3<f64>,
    pub shape: CellShape,
    pub matrix_inv_transposed: Option<Matrix3<f64>>,
}

impl PartialEq for UnitCell {
    fn eq(&self, other: &Self) -> bool {
        self.matrix
            .iter()
            .zip(other.matrix.iter())
            .all(|(a, b)| (a - b).abs() < f64::EPSILON)
    }
}

impl UnitCell {
    #[must_use]
    pub fn new() -> Self {
        UnitCell::new_from_lengths([0.0, 0.0, 0.0])
    }

    #[must_use]
    pub fn new_from_lengths(lengths: Vec3D) -> Self {
        UnitCell::new_from_lengths_angles(lengths, &mut [90.0, 90.0, 90.0]).unwrap()
    }

    /// # Errors
    ///
    /// Returns an error if:
    /// - Any length is negative
    /// - Any angle is negative
    /// - Any angle is zero
    /// - Any angle is 180 degrees or greater
    pub fn new_from_lengths_angles(lengths: Vec3D, angles: &mut Vec3D) -> Result<Self, CError> {
        let matrix = matrix::cell_matrix_from_lengths_angles(lengths, angles)?;
        Self::new_from_matrix(matrix)
    }

    /// # Errors
    ///
    /// Will return `Err` if the determinant of the input matrix < 0.0
    ///
    /// Will return `Err` if the input matrix [`CellShape::Orthorhombic`], but
    /// not diagonal
    pub fn new_from_matrix(matrix: Matrix3<f64>) -> Result<Self, CError> {
        if matrix.determinant() < 0.0 {
            return Err(CError::GenericError(
                "invalid unit cell matrix with negative determinant".to_string(),
            ));
        }

        let lengths = utils::calc_lengths_from_matrix(matrix);
        let angles = utils::calc_angles_from_matrix(matrix);

        let is_diagonal_matrix = utils::is_diagonal(matrix);
        if !is_diagonal_matrix && utils::is_orthorhombic(lengths, angles) {
            return Err(CError::GenericError("orthorhombic cell must have their a vector along x axis, b vector along y axis and c vector along z axis".to_string()));
        }

        let shape = if is_diagonal_matrix {
            if matrix.diagonal().iter().all(|&x| utils::is_roughly_zero(x)) {
                CellShape::Infinite
            } else {
                CellShape::Orthorhombic
            }
        } else {
            CellShape::Triclinic
        };

        let matrix_inv_transposed = matrix.try_inverse().map(|m| m.transpose());

        Ok(UnitCell {
            matrix,
            shape,
            matrix_inv_transposed,
        })
    }

    #[must_use]
    pub fn lengths(&self) -> Vec3D {
        match self.shape {
            CellShape::Infinite => [0.0, 0.0, 0.0],
            CellShape::Orthorhombic => [
                self.matrix[(0, 0)],
                self.matrix[(1, 1)],
                self.matrix[(2, 2)],
            ],
            CellShape::Triclinic => utils::calc_lengths_from_matrix(self.matrix),
        }
    }

    #[must_use]
    pub fn angles(&self) -> Vec3D {
        match self.shape {
            CellShape::Infinite | CellShape::Orthorhombic => [90.0, 90.0, 90.0],
            CellShape::Triclinic => utils::calc_angles_from_matrix(self.matrix),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use nalgebra::Matrix3;

    #[test]
    fn test_unit_cell_new() {
        let cell = UnitCell::new();
        assert_eq!(cell.matrix, Matrix3::zeros());
        assert_eq!(cell.shape, CellShape::Infinite);
        assert!(cell.matrix_inv_transposed.is_none());
    }

    #[test]
    fn test_unit_cell_from_lengths() {
        let lengths = [10.0, 20.0, 30.0];
        let cell = UnitCell::new_from_lengths(lengths);
        assert_approx_eq!(cell.matrix[(0, 0)], 10.0);
        assert_approx_eq!(cell.matrix[(1, 1)], 20.0);
        assert_approx_eq!(cell.matrix[(2, 2)], 30.0);
        assert_eq!(cell.shape, CellShape::Orthorhombic);
    }

    #[test]
    fn test_unit_cell_from_lengths_angles() {
        let lengths = [10.0, 20.0, 30.0];
        let mut angles = [90.0, 90.0, 90.0];
        let cell = UnitCell::new_from_lengths_angles(lengths, &mut angles).unwrap();
        assert_approx_eq!(cell.matrix[(0, 0)], 10.0);
        assert_approx_eq!(cell.matrix[(1, 1)], 20.0);
        assert_approx_eq!(cell.matrix[(2, 2)], 30.0);
        assert_eq!(cell.shape, CellShape::Orthorhombic);
    }

    #[test]
    fn test_unit_cell_from_matrix() {
        let mut matrix = Matrix3::zeros();
        matrix[(0, 0)] = 10.0;
        matrix[(1, 1)] = 20.0;
        matrix[(2, 2)] = 30.0;
        let cell = UnitCell::new_from_matrix(matrix).unwrap();
        assert_eq!(cell.matrix, matrix);
        assert_eq!(cell.shape, CellShape::Orthorhombic);
    }

    #[test]
    fn test_unit_cell_lengths() {
        let cell = UnitCell::new_from_lengths([10.0, 20.0, 30.0]);
        let expected = [10.0, 20.0, 30.0];
        for (t, e) in expected.iter().zip(cell.lengths()) {
            assert_approx_eq!(t, e);
        }
    }

    #[test]
    fn test_unit_cell_angles() {
        let cell = UnitCell::new_from_lengths([10.0, 20.0, 30.0]);
        for i in 0..3 {
            assert_approx_eq!(cell.angles()[i], 90.0);
        }
    }

    #[test]
    fn test_check_lengths_valid() {
        let lengths = [1.0, 2.0, 3.0];
        assert!(validation::check_lengths(&lengths).is_ok());
    }

    #[test]
    #[should_panic(expected = "negative length found: -1")]
    fn test_check_lengths_invalid_first_negative() {
        let lengths = [-1.0, 2.0, 3.0];
        validation::check_lengths(&lengths).unwrap();
    }

    #[test]
    fn test_check_angles_valid() {
        let angles = [90.0, 90.0, 90.0];
        assert!(validation::check_angles(&angles).is_ok());
    }

    #[test]
    fn test_check_angles_valid_non_orthogonal() {
        let angles = [60.0, 60.0, 60.0];
        assert!(validation::check_angles(&angles).is_ok());
    }

    #[test]
    #[should_panic(expected = "negative angle found: -90")]
    fn test_check_angles_invalid_negative() {
        let angles = [-90.0, 90.0, 90.0];
        validation::check_angles(&angles).unwrap();
    }

    #[test]
    #[should_panic(expected = "zero angle found: 0")]
    fn test_check_angles_invalid_zero() {
        let angles = [0.0, 90.0, 90.0];
        validation::check_angles(&angles).unwrap();
    }

    #[test]
    #[should_panic(expected = "angle too large: 180")]
    fn test_check_angles_invalid_180_and_greater() {
        let angles = [180.0, 90.0, 90.0];
        validation::check_angles(&angles).unwrap();
        let angles = [200.0, 90.0, 90.0];
        validation::check_angles(&angles).unwrap();
    }

    #[test]
    fn test_cell_matrix_from_lengths_angles_valid() {
        let lengths = [10.0, 20.0, 30.0];
        let mut angles = [90.0, 90.0, 90.0];
        let matrix = matrix::cell_matrix_from_lengths_angles(lengths, &mut angles).unwrap();
        assert_approx_eq!(matrix[(0, 0)], 10.0);
        assert_approx_eq!(matrix[(1, 1)], 20.0);
        assert_approx_eq!(matrix[(2, 2)], 30.0);
    }

    #[test]
    fn test_cell_matrix_from_lengths_angles_triclinic() {
        let lengths = [10.0, 20.0, 30.0];
        let mut angles = [60.0, 60.0, 60.0];
        let matrix = matrix::cell_matrix_from_lengths_angles(lengths, &mut angles).unwrap();
        assert!(!utils::is_diagonal(matrix));
    }

    #[test]
    fn test_is_orthorhombic() {
        let lengths = [10.0, 20.0, 30.0];
        let angles = [90.0, 90.0, 90.0];
        assert!(utils::is_orthorhombic(lengths, angles));
    }

    #[test]
    fn test_is_orthorhombic_with_nan() {
        let lengths = [10.0, 20.0, 0.0];
        let angles = [90.0, 90.0, f64::NAN];
        assert!(utils::is_orthorhombic(lengths, angles));
    }

    #[test]
    fn test_is_infinite() {
        let lengths = [0.0, 0.0, 0.0];
        assert!(utils::is_infinite(lengths));
    }

    #[test]
    fn test_is_diagonal() {
        let mut matrix = Matrix3::zeros();
        matrix[(0, 0)] = 1.0;
        matrix[(1, 1)] = 2.0;
        matrix[(2, 2)] = 3.0;
        assert!(utils::is_diagonal(matrix));
    }

    #[test]
    fn test_is_not_diagonal() {
        let mut matrix = Matrix3::zeros();
        matrix[(0, 0)] = 1.0;
        matrix[(1, 1)] = 2.0;
        matrix[(2, 2)] = 3.0;
        matrix[(0, 1)] = 0.1;
        assert!(!utils::is_diagonal(matrix));
    }

    #[test]
    fn from_lengths_angles() {
        let expected = UnitCell::new_from_lengths_angles(
            [8.431_160_35, 14.505_106_13, 15.609_114_68],
            &mut [73.316_992_12, 85.702_005_82, 89.375_015_29],
        )
        .unwrap();
        let mut true_cell = UnitCell::new();
        true_cell.matrix[(0, 0)] = 8.431_160_35;
        true_cell.matrix[(1, 0)] = 0.158_219_155_128;
        true_cell.matrix[(1, 1)] = 14.504_243_186_3;
        true_cell.matrix[(2, 0)] = 1.169_806_636_24;
        true_cell.matrix[(2, 1)] = 4.468_514_985_5;
        true_cell.matrix[(2, 2)] = 14.910_009_640_5;
        let diff = expected.matrix - true_cell.matrix;
        assert!((diff).iter().all(|&x| x.abs() < 1e-6), "diff: {diff}");
    }
}
