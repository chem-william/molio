use molio::{
    extendedxyzparser::{ExtendedXyzParser, Property},
    unit_cell::UnitCell,
};
use std::{
    fs::File,
    io::{BufRead, BufReader, Seek},
    path::Path,
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum CError {
    #[error("Unsupported file format: `{0}`")]
    UnsupportedFileFormat(String),
    #[error("{0}")]
    IoError(#[from] std::io::Error),
    #[error("generic file error")]
    GenericFileError,
    #[error("{format} format: not enough lines at step {step} (expected {expected}, got {got})")]
    UnexpectedEof {
        format: String,
        step: usize,
        expected: usize,
        got: usize,
    },
    #[error("unknown format: {0}")]
    UnknownFormat(String),
    #[error("")]
    UnexpectedSymbol,
    #[error("Failed to parse float: {0}")]
    ParseFloatError(#[from] std::num::ParseFloatError),
    #[error("Missing token")]
    MissingToken,
}

#[derive(Debug)]
pub struct Atom {
    pub id: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Default)]
pub struct Frame {
    pub atoms: Vec<Atom>,
    pub unit_cell: UnitCell,
}

impl Frame {
    pub fn new() -> Self {
        Frame {
            atoms: vec![],
            unit_cell: UnitCell::new(),
        }
    }
    pub fn size(&self) -> usize {
        self.atoms.len()
    }

    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.atoms.iter().map(|a| [a.x, a.y, a.z]).collect()
    }
    pub fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom)
    }
}

pub trait FileFormat {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError>;
    // fn read(&self) -> Result<Frame, CError>;
    // fn read_at(&mut self, index: usize) -> Result<Frame, CError>;
    fn write(&self, path: &Path, frame: &Frame) -> Result<(), CError>;
    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError>;
}
pub struct XYZFormat;
enum PropertyKind {
    Bool,
    Double,
    String,
    Vector3D,
}
struct ExtendedProperty {
    name: String,
    kind: PropertyKind,
}

type PropertiesList = Vec<ExtendedProperty>;
impl XYZFormat {
    fn read_extended_comment_line(line: &str, frame: &mut Frame) -> Result<PropertiesList, CError> {
        println!("what is the line: {line}");
        let contains_properties = line.contains("species:S:1:pos:R:3");
        let contains_lattice = line.contains("Lattice");

        if !(contains_properties || contains_lattice) {
            return Ok(Vec::new());
        }

        let mut extxyz_parser = ExtendedXyzParser::new(line);
        let properties = extxyz_parser.parse();
        println!("are there properties: {:?}", properties);

        for (name, value) in &properties {
            if name == "Lattice" || name == "Properties" {
                continue;
            }
            // frame.set(name.clone(), value.clone());
        }

        if let Some(Property::Str(lattice)) = properties.get("Lattice") {
            let cell = UnitCell::parse(lattice);
            frame.unit_cell = cell;
        }

        // if let Some(prop_string) = properties.get("Properties") {
        //     return Ok(parse_property_list(prop_string.as_string()?));
        // }

        Ok(Vec::new())
    }
}
impl FileFormat for XYZFormat {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        println!("Reading as XYZ format");
        let mut line = String::new();
        let _ = reader.read_line(&mut line)?;

        let n_atoms = line
            .trim()
            .parse::<usize>()
            .expect("expected number of atoms");

        line.clear();
        let _ = reader.read_line(&mut line)?;
        let mut frame = Frame::new();
        let properties = XYZFormat::read_extended_comment_line(&line, &mut frame);

        for _ in 0..n_atoms {
            line.clear();
            let _ = reader.read_line(&mut line)?;
            let mut tokens = line.split_whitespace();

            let name = tokens.next().ok_or(CError::UnexpectedSymbol)?.to_string();

            let x: f64 = tokens.next().ok_or(CError::MissingToken)?.parse()?;
            let y: f64 = tokens.next().ok_or(CError::MissingToken)?.parse()?;
            let z: f64 = tokens.next().ok_or(CError::MissingToken)?.parse()?;

            assert_eq!(tokens.next(), None, "extxyz not yet implemented");

            let atom = Atom { x, y, z, id: name };
            // let remaining_line = tokens.collect::<Vec<_>>().join(" ");
            // read_atomic_properties(&properties, &remaining_line, &mut atom)?;

            frame.add_atom(atom);
        }

        Ok(frame)
    }

    // fn read(&self) -> Result<Frame, CError> {
    //     println!("Reading as XYZ format");
    //     Ok(Frame { atoms: vec![] })
    // }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        println!("calling forward from XYZ");
        let mut line_position = 0;

        let mut line = String::new();

        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 || line.trim().is_empty() {
            return Ok(None);
        }
        let n_atoms: usize = line.trim().parse().unwrap();

        for i in 0..=n_atoms {
            line.clear();
            let bytes = reader.read_line(&mut line)?;
            if bytes == 0 {
                return Err(CError::UnexpectedEof {
                    format: "XYZ".to_string(),
                    step: line_position,
                    expected: n_atoms + 2, // first count line + n_atoms atom lines + blank/comment?
                    got: i + 1,            // how many we actually read
                });
            }
        }
        let position = reader.stream_position().unwrap();
        line_position += 1;
        Ok(Some(position))
    }

    fn write(&self, path: &Path, frame: &Frame) -> Result<(), CError> {
        println!(
            "Writing {:?} as XYZ format with {} atoms",
            path,
            frame.size()
        );
        Ok(())
    }
}

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

pub enum TextFormat {
    XYZ,
    // PDB,
    Guess,
}

pub enum Format {
    XYZ(XYZFormat),
    // PDB(PDBFormat),
}

impl Format {
    pub fn new(path: &Path) -> Result<Self, CError> {
        let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");

        match ext.to_lowercase().as_str() {
            "xyz" => Ok(Format::XYZ(XYZFormat)),
            _ => Err(CError::GenericFileError),
        }
    }
    pub fn new_from_format(fmt: TextFormat, path: &Path) -> Result<Self, CError> {
        match fmt {
            TextFormat::XYZ => Ok(Format::XYZ(XYZFormat)),
            TextFormat::Guess => Self::new(path),
        }
    }
}

impl FileFormat for Format {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        match self {
            Format::XYZ(format) => format.read_next(reader),
        }
    }
    // fn read(&self) -> Result<Frame, CError> {
    //     match self {
    //         Format::XYZ(format) => format.read(),
    //     }
    // }

    fn write(&self, path: &Path, frame: &Frame) -> Result<(), CError> {
        match self {
            Format::XYZ(format) => format.write(path, frame),
        }
    }
    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        match self {
            Format::XYZ(format) => format.forward(reader),
        }
    }
}

pub struct Trajectory<'a> {
    pub frames: Vec<Frame>,

    path: &'a Path,
    strategy: Format,
    reader: BufReader<File>,
    size: usize,
    frame_positions: Vec<u64>,
}

impl<'a> Trajectory<'a> {
    /// Constructs a `Trajectory` that will guess the format.
    pub fn new(path: &'a Path) -> Result<Self, CError> {
        Self::with_format(path, TextFormat::Guess)
    }

    /// Constructs a `Trajectory` with an explicitly provided format.
    pub fn with_format(path: &'a Path, fmt: TextFormat) -> Result<Self, CError> {
        let strategy = Format::new_from_format(fmt, path)?;

        let file = File::open(path).map_err(CError::IoError)?;
        let mut reader = BufReader::new(file);
        let frame_positions = Trajectory::scan_all(&mut reader, &strategy);
        let size = frame_positions.len() - 1;

        Ok(Trajectory {
            frames: vec![],
            path,
            strategy,
            frame_positions,
            reader,
            size,
        })
    }

    // pub fn read(&mut self) -> Result<Frame, CError> {
    //     self.strategy.read()
    // }

    pub fn read_at(&mut self, index: usize) -> Result<Frame, CError> {
        if index >= self.frame_positions.len() {
            if self.frame_positions.is_empty() {
                return Err(CError::GenericFileError);
            }
            return Err(CError::GenericFileError);
        }
        self.reader
            .seek(std::io::SeekFrom::Start(self.frame_positions[index]))
            .unwrap();
        self.strategy.read_next(&mut self.reader)
    }

    pub fn write(&self, frame: &Frame) -> Result<(), CError> {
        self.strategy.write(self.path, frame)
    }

    fn scan_all(reader: &mut BufReader<File>, strategy: &Format) -> Vec<u64> {
        let mut frame_positions = vec![0];
        while let Some(pos) = strategy.forward(reader).unwrap() {
            frame_positions.push(pos)
        }
        reader.rewind().unwrap();
        frame_positions
    }
}

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
    use super::*;

    #[test]
    fn trajectory() {
        let path = Path::new("./src/tests-data/xyz/extended.xyz");
        let mut trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 3);

        let frame = trajectory.read_at(0).unwrap();
        assert_eq!(frame.size(), 192);

        let frame = trajectory.read_at(1).unwrap();
        assert_eq!(frame.size(), 62);

        let frame = trajectory.read_at(0).unwrap();
        assert_eq!(frame.size(), 192);

        let frame = trajectory.read_at(2).unwrap();
        assert_eq!(frame.size(), 8);

        let mut unit_cell = UnitCell::new();
        unit_cell.cell_matrix[(0, 0)] = 4.0;
        unit_cell.cell_matrix[(1, 1)] = 7.0;
        unit_cell.cell_matrix[(2, 2)] = 3.0;

        assert_eq!(frame.unit_cell, unit_cell);
    }
}
