// pub struct PDBFormat;
// impl FileFormat for PDBFormat {
//     fn read_next(&self, path: &Path) -> Result<Frame, CError> {
//         println!("Reading {:?} as PDB format", path);
//         // Replace with real parsing logic.
//         Ok(Frame { atoms: vec![] })
//     }

//     fn write(&self, path: &Path, frame: &Frame) -> Result<(), CError> {
//         println!(
//             "Writing {:?} as PDB format with {} atoms",
//             path,
//             frame.size()
//         );
//         Ok(())
//     }
//     fn read(&self) -> Result<Frame, CError> {
//         println!("Reading as PDB format");
//         // Replace with real parsing logic.
//         Ok(Frame { atoms: vec![] })
//     }
// }

// fn main() -> Result<(), CError> {
//     let path = Path::new("structure.xyz");
//     let mut trajectory = Trajectory::new(path)?;

//     let frame = trajectory.read_at(0)?;
//     println!(
//         "There are {} atoms in the first frame using read_at(0)",
//         frame.size()
//     );
//     let frame = trajectory.read_at(1)?;
//     println!(
//         "There are {} atoms in the second frame using read_at(1)",
//         frame.size()
//     );
//     let frame = trajectory.read_at(0)?;
//     println!(
//         "There are {} atoms in the first frame using read_at(0)",
//         frame.size()
//     );
//     let frame = trajectory.read_at(2)?;
//     println!(
//         "There are {} atoms in the first frame using read_at(2)",
//         frame.size()
//     );

//     Ok(())
// }

fn main() {
    println!("hellow world");
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use assert_approx_eq::assert_approx_eq;
    use molio::{trajectory::Trajectory, unit_cell::UnitCell};

    #[test]
    fn trajectory() {
        let path = Path::new("./src/tests-data/xyz/extended.xyz");
        let mut trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 3);

        let frame = trajectory.read_at(0).unwrap();
        assert_eq!(frame.size(), 192);
        let mut unit_cell = UnitCell::new();
        unit_cell.cell_matrix[(0, 0)] = 8.43116035;
        unit_cell.cell_matrix[(0, 1)] = 0.158219155128;
        unit_cell.cell_matrix[(1, 1)] = 14.5042431863;
        unit_cell.cell_matrix[(0, 2)] = 1.16980663624;
        unit_cell.cell_matrix[(1, 2)] = 4.4685149855;
        unit_cell.cell_matrix[(2, 2)] = 14.9100096405;
        assert_eq!(frame.unit_cell, unit_cell);

        let frame = trajectory.read_at(1).unwrap();
        assert_eq!(frame.size(), 62);

        let frame = trajectory.read_at(0).unwrap();
        assert_eq!(frame.size(), 192);

        // Atom level properties
        let positions = frame.positions()[0];
        assert_approx_eq!(positions[0], 2.33827271799, 1e-9);
        assert_approx_eq!(positions[1], 4.55315540425, 1e-9);
        assert_approx_eq!(positions[2], 11.5841360926, 1e-9);
        assert_approx_eq!(frame.atoms[0].properties["CS_0"].expect_double(), 24.10);
        assert_approx_eq!(frame.atoms[0].properties["CS_1"].expect_double(), 31.34);

        // Frame level properties
        assert_eq!(frame.properties["ENERGY"].expect_double(), -2069.84934116);
        assert_eq!(frame.properties["Natoms"].expect_double(), 192.0);
        assert_eq!(frame.properties["NAME"].expect_string(), "COBHUW");
        assert!(frame.properties["IsStrange"].expect_bool());

        let frame = trajectory.read_at(2).unwrap();
        assert_eq!(frame.size(), 8);

        let mut unit_cell = UnitCell::new();
        unit_cell.cell_matrix[(0, 0)] = 4.0;
        unit_cell.cell_matrix[(1, 1)] = 7.0;
        unit_cell.cell_matrix[(2, 2)] = 3.0;

        assert_eq!(frame.unit_cell, unit_cell);
    }
}
