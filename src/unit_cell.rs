use nalgebra::Matrix3;

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
}
