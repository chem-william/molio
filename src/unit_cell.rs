use crate::error::CError;
use core::f64;
use nalgebra::Matrix3;

type Vec3D = [f64; 3];

#[derive(Default, Debug)]
pub struct UnitCell {
    pub cell_matrix: Matrix3<f64>,
}

impl PartialEq for UnitCell {
    fn eq(&self, other: &Self) -> bool {
        self.cell_matrix
            .iter()
            .zip(other.cell_matrix.iter())
            .all(|(a, b)| (a - b).abs() < f64::EPSILON)
    }
}

impl UnitCell {
    const EPSILON: f64 = 1e-5;

    fn deg2rad(x: f64) -> f64 {
        x * f64::consts::PI / 180.0
    }

    fn rad2deg(x: f64) -> f64 {
        x * 180.0 / f64::consts::PI
    }

    fn cos_degree(theta: f64) -> f64 {
        Self::deg2rad(theta).cos()
    }

    fn sin_degree(theta: f64) -> f64 {
        Self::deg2rad(theta).sin()
    }

    pub fn new() -> Self {
        UnitCell {
            cell_matrix: Matrix3::zeros(),
        }
    }

    pub fn parse(lattice: &str) -> Self {
        let mut cell_matrix = Matrix3::zeros();

        cell_matrix
            .iter_mut()
            .zip(lattice.split_whitespace())
            .for_each(|(matrix_entry, lattice_item)| {
                *matrix_entry = lattice_item.parse::<f64>().expect("expected float");
            });
        UnitCell { cell_matrix }
    }

    fn check_lengths(lengths: &Vec3D) -> Result<(), CError> {
        if lengths.iter().any(|&x| x < 0.0) {
            return Err(CError::GenericError(
                "lengths cannot be negative".to_string(),
            ));
        };

        Ok(())
    }

    fn check_angles(angles: &Vec3D) -> Result<(), CError> {
        if angles.iter().any(|&x| x < 0.0) {
            return Err(CError::GenericError(
                "angles cannot be negative".to_string(),
            ));
        };

        if angles.iter().any(|&x| x.abs() < Self::EPSILON) {
            return Err(CError::GenericError(
                "angles cannot be (roughly) zero".to_string(),
            ));
        }

        if angles.iter().any(|&x| x >= 180.0) {
            return Err(CError::GenericError(
                "angles cannot be larger than or equal to 180 degrees".to_string(),
            ));
        }

        Ok(())
    }
    pub fn cell_matrix_from_lengths_angles(
        lengths: Vec3D,
        angles: &mut Vec3D,
    ) -> Result<Self, CError> {
        Self::check_lengths(&lengths)?;
        Self::check_angles(angles)?;

        if angles.iter().all(|&x| (x - 90.0).abs() < 1e-3) {
            angles.iter_mut().for_each(|x| *x = 90.0);
        }
        let mut cell_matrix: Matrix3<f64> = Matrix3::zeros();
        cell_matrix[(0, 0)] = lengths[0];

        cell_matrix[(1, 0)] = Self::cos_degree(angles[2]) * lengths[1];
        cell_matrix[(1, 1)] = Self::sin_degree(angles[2]) * lengths[1];

        cell_matrix[(2, 0)] = Self::cos_degree(angles[1]);
        cell_matrix[(2, 1)] = (Self::cos_degree(angles[0])
            - Self::cos_degree(angles[1]) * Self::cos_degree(angles[2]))
            / Self::sin_degree(angles[2]);
        cell_matrix[(2, 2)] = (1.0
            - cell_matrix[(2, 0)] * cell_matrix[(2, 0)]
            - cell_matrix[(2, 1)] * cell_matrix[(2, 1)])
            .sqrt();
        cell_matrix[(2, 0)] *= lengths[2];
        cell_matrix[(2, 1)] *= lengths[2];
        cell_matrix[(2, 2)] *= lengths[2];

        Ok(UnitCell { cell_matrix })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Matrix3;

    #[test]
    fn test_unit_cell_new() {
        let cell = UnitCell::new();
        assert_eq!(cell.cell_matrix, Matrix3::zeros());
    }

    #[test]
    fn test_unit_cell_parse() {
        let lattice = "1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0";
        let cell = UnitCell::parse(lattice);

        let expected = Matrix3::new(1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0);
        assert_eq!(cell.cell_matrix, expected);
    }

    #[test]
    fn test_check_lengths_valid() {
        let lengths = [1.0, 2.0, 3.0];
        assert!(UnitCell::check_lengths(&lengths).is_ok());
    }

    #[test]
    #[should_panic(expected = "lengths cannot be negative")]
    fn test_check_lengths_invalid() {
        let lengths = [-1.0, -2.0, -3.0];
        UnitCell::check_lengths(&lengths).unwrap();
    }

    #[test]
    fn test_check_angles_valid() {
        let angles = [90.0, 90.0, 90.0];
        assert!(UnitCell::check_angles(&angles).is_ok());
    }

    #[test]
    fn from_lengths_angles() {
        let expected = UnitCell::cell_matrix_from_lengths_angles(
            [8.43116035, 14.50510613, 15.60911468],
            &mut [73.31699212, 85.70200582, 89.37501529],
        )
        .unwrap();
        let mut true_cell = UnitCell::new();
        true_cell.cell_matrix[(0, 0)] = 8.43116035;
        true_cell.cell_matrix[(1, 0)] = 0.158219155128;
        true_cell.cell_matrix[(1, 1)] = 14.5042431863;
        true_cell.cell_matrix[(2, 0)] = 1.16980663624;
        true_cell.cell_matrix[(2, 1)] = 4.4685149855;
        true_cell.cell_matrix[(2, 2)] = 14.9100096405;
        let diff = expected.cell_matrix - true_cell.cell_matrix;
        assert!((diff).iter().all(|&x| x.abs() < 1e-6), "diff: {diff}");
    }

    #[test]
    #[should_panic(expected = "angles cannot be negative")]
    fn test_check_angles_invalid_negative() {
        let angles = [-90.0, -90.0, -90.0];
        UnitCell::check_angles(&angles).unwrap();
    }

    #[test]
    #[should_panic(expected = "angles cannot be (roughly) zero")]
    fn test_check_angles_invalid_zero() {
        let angles = [0.0, 0.0, 0.0];
        UnitCell::check_angles(&angles).unwrap();
    }

    #[test]
    #[should_panic(expected = "angles cannot be larger than or equal to 180 degrees")]
    fn test_check_angles_invalid_180() {
        let angles = [180.0, 180.0, 180.0];
        UnitCell::check_angles(&angles).unwrap();
    }

    #[test]
    #[should_panic(expected = "lengths cannot be negative")]
    fn test_cell_matrix_from_lengths_angles_invalid() {
        let lengths = [-10.0, 20.0, 30.0];
        let mut angles = [90.0, 90.0, 90.0];
        UnitCell::cell_matrix_from_lengths_angles(lengths, &mut angles).unwrap();
    }
}
