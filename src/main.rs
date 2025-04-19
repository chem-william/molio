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
