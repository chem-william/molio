use crate::atom::Atom;
use crate::error::CError;
use crate::extendedxyzparser::ExtendedXyzParser;
use crate::format::FileFormat;
use crate::frame::Frame;
use crate::property::Properties;
use crate::property::{Property, PropertyKind};
use crate::unit_cell::UnitCell;
use nalgebra::Matrix3;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek};
use std::path::Path;
use std::str::SplitWhitespace;

pub struct XYZFormat;
type PropertiesList = BTreeMap<String, PropertyKind>;

// Helper trait for parsing values
trait ValueParser {
    fn parse(value: &str) -> Result<Property, CError>;
}

// Implement parsers for each property type
struct StringParser;
struct BoolParser;
struct DoubleParser;
struct Vector3DParser;
struct Matrix3x3Parser;
struct VectorXDParser;

impl ValueParser for StringParser {
    fn parse(value: &str) -> Result<Property, CError> {
        Ok(Property::String(value.to_string()))
    }
}

impl ValueParser for BoolParser {
    fn parse(value: &str) -> Result<Property, CError> {
        match value.to_lowercase().as_str() {
            "t" | "true" => Ok(Property::Bool(true)),
            "f" | "false" => Ok(Property::Bool(false)),
            _ => Err(CError::GenericError(format!(
                "Invalid boolean value: {}",
                value
            ))),
        }
    }
}

impl ValueParser for DoubleParser {
    fn parse(value: &str) -> Result<Property, CError> {
        value
            .parse::<f64>()
            .map(Property::Double)
            .map_err(|e| CError::GenericError(format!("Failed to parse number: {}", e)))
    }
}

impl ValueParser for Vector3DParser {
    fn parse(value: &str) -> Result<Property, CError> {
        let parts: Vec<&str> = value.split_whitespace().collect();
        if parts.len() != 3 {
            return Err(CError::GenericError(format!(
                "Vector3D requires exactly 3 components, got {}: {:?}",
                parts.len(),
                parts
            )));
        }
        let x = parts[0]
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("Failed to parse x component: {}", e)))?;
        let y = parts[1]
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("Failed to parse y component: {}", e)))?;
        let z = parts[2]
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("Failed to parse z component: {}", e)))?;
        Ok(Property::Vector3D([x, y, z]))
    }
}

impl ValueParser for Matrix3x3Parser {
    fn parse(value: &str) -> Result<Property, CError> {
        let parts: Vec<&str> = value.split_whitespace().collect();
        if parts.len() != 9 {
            return Err(CError::GenericError(format!(
                "Matrix3x3 requires exactly 9 components, got {}: {:?}",
                parts.len(),
                parts
            )));
        }
        let nums: Result<Vec<f64>, _> = parts.iter().map(|p| p.parse::<f64>()).collect();
        nums.map(|n| Property::Matrix3x3(Matrix3::from_iterator(n)))
            .map_err(|e| CError::GenericError(format!("Failed to parse matrix components: {}", e)))
    }
}

impl ValueParser for VectorXDParser {
    fn parse(value: &str) -> Result<Property, CError> {
        let parts: Vec<&str> = value.split_whitespace().collect();
        let nums: Result<Vec<f64>, _> = parts.iter().map(|p| p.parse::<f64>()).collect();
        nums.map(Property::VectorXD)
            .map_err(|e| CError::GenericError(format!("Failed to parse vector components: {}", e)))
    }
}

impl XYZFormat {
    fn parse_value(value: &str, kind: PropertyKind) -> Result<Property, CError> {
        match kind {
            PropertyKind::String => StringParser::parse(value),
            PropertyKind::Bool => BoolParser::parse(value),
            PropertyKind::Double => DoubleParser::parse(value),
            PropertyKind::Vector3D => Vector3DParser::parse(value),
            PropertyKind::Matrix3x3 => Matrix3x3Parser::parse(value),
            PropertyKind::VectorXD => VectorXDParser::parse(value),
        }
    }

    fn read_atomic_properties(
        properties: &PropertiesList,
        tokens: &mut SplitWhitespace,
        atom: &mut Atom,
    ) -> Result<(), CError> {
        for (name, kind) in properties {
            match kind {
                PropertyKind::Vector3D => {
                    // For Vector3D, collect all three components first
                    let x = tokens.next().ok_or(CError::MissingToken)?;
                    let y = tokens.next().ok_or(CError::MissingToken)?;
                    let z = tokens.next().ok_or(CError::MissingToken)?;
                    let value = format!("{} {} {}", x, y, z);
                    let property = Self::parse_value(&value, kind.clone())?;
                    atom.properties.insert(name.clone(), property);
                }
                _ => {
                    let value = tokens.next().ok_or(CError::MissingToken)?;
                    let property = Self::parse_value(value, kind.clone())?;
                    atom.properties.insert(name.clone(), property);
                }
            }
        }
        Ok(())
    }

    fn parse_property_list(line: &str) -> Result<PropertiesList, CError> {
        const PREFIX: &str = "species:S:1:pos:R:3:";

        let initial_input = &line;
        let rest = line.strip_prefix(PREFIX).ok_or_else(|| {
            CError::GenericError(
                "Invalid property list format: missing expected prefix".to_string(),
            )
        })?;

        let fields: Vec<&str> = rest.split(':').collect();

        if fields.len() % 3 != 0 {
            return Err(CError::GenericError(
                "Invalid property list format: property definitions must be in groups of 3 (name:type:count)".to_string(),
            ));
        }

        let mut properties = PropertiesList::new();

        for chunk in fields.chunks_exact(3) {
            let name = chunk[0];
            let kind = match chunk[1] {
                "R" | "I" => PropertyKind::Double,
                "S" => PropertyKind::String,
                "L" => PropertyKind::Bool,
                unknown => {
                    return Err(CError::GenericError(format!(
                        "Unknown property type: {}",
                        unknown
                    )));
                }
            };

            let count = chunk[2].parse::<usize>().map_err(|e| {
                CError::GenericError(format!("Invalid property count '{}': {}", chunk[2], e))
            })?;

            if count == 0 {
                return Err(CError::GenericError(format!(
                    "Invalid count of 0 for property '{}'",
                    name
                )));
            }

            match (count, &kind) {
                (3, PropertyKind::Double) => {
                    properties.insert(name.to_string(), PropertyKind::Vector3D);
                }
                (1, k) => {
                    properties.insert(name.to_string(), k.clone());
                }
                (n, k) => {
                    for i in 0..n {
                        properties.insert(format!("{name}_{i}"), k.clone());
                    }
                }
            }
        }

        Ok(properties)
    }

    fn parse_frame_property(value: &str) -> Property {
        // Try parsing as bool first
        let lowercased = value.to_lowercase();
        if ["t", "true", "f", "false"].contains(&lowercased.as_str()) {
            return Self::parse_value(&lowercased, PropertyKind::Bool)
                .unwrap_or_else(|_| Property::String(value.to_string()));
        }

        // Try parsing as vector/matrix
        let parts: Vec<&str> = value.split_whitespace().collect();
        match parts.len() {
            1 => Self::parse_value(value, PropertyKind::Double)
                .unwrap_or_else(|_| Property::String(value.to_string())),
            3 => Self::parse_value(value, PropertyKind::Vector3D)
                .unwrap_or_else(|_| Property::String(value.to_string())),
            9 => Self::parse_value(value, PropertyKind::Matrix3x3)
                .unwrap_or_else(|_| Property::String(value.to_string())),
            _ => Self::parse_value(value, PropertyKind::VectorXD)
                .unwrap_or_else(|_| Property::String(value.to_string())),
        }
    }

    fn read_extended_comment_line(line: &str, frame: &mut Frame) -> Result<PropertiesList, CError> {
        if !(line.contains("species:S:1:pos:R:3") || line.contains("Lattice")) {
            return Ok(PropertiesList::new());
        }

        let extxyz_parser = ExtendedXyzParser::new(line);
        let properties = extxyz_parser.parse();

        for (k, v) in properties
            .iter()
            .filter(|(k, _)| k.as_str() != "Lattice" && k.as_str() != "Properties")
        {
            frame
                .properties
                .insert(k.clone(), Self::parse_frame_property(v));
        }

        if let Some(lattice) = properties.get("Lattice") {
            frame.unit_cell = UnitCell::parse(lattice);
        }

        if let Some(prop_string) = properties.get("Properties") {
            return Self::parse_property_list(prop_string);
        }

        Ok(PropertiesList::new())
    }
}

impl FileFormat for XYZFormat {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        let mut line = String::new();
        let _ = reader.read_line(&mut line)?;

        let n_atoms = line
            .trim()
            .parse::<usize>()
            .expect("expected number of atoms");

        line.clear();
        let _ = reader.read_line(&mut line)?;
        let mut frame = Frame::new();
        let properties = XYZFormat::read_extended_comment_line(&line, &mut frame)?;

        for _ in 0..n_atoms {
            line.clear();
            let _ = reader.read_line(&mut line)?;
            let mut tokens = line.split_whitespace();

            let name = tokens.next().ok_or(CError::UnexpectedSymbol)?.to_string();

            let x: f64 = tokens.next().ok_or(CError::MissingToken)?.parse()?;
            let y: f64 = tokens.next().ok_or(CError::MissingToken)?.parse()?;
            let z: f64 = tokens.next().ok_or(CError::MissingToken)?.parse()?;

            let mut atom = Atom {
                x,
                y,
                z,
                symbol: name,
                properties: Properties::new(),
            };
            XYZFormat::read_atomic_properties(&properties, &mut tokens, &mut atom)?;

            frame.add_atom(atom);
        }

        Ok(frame)
    }

    fn read(&self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        // TODO: replace with has_data_left when stabilized
        if reader.fill_buf().map(|b| !b.is_empty()).unwrap() {
            Ok(Some(self.read_next(reader).unwrap()))
        } else {
            Ok(None)
        }
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
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

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::{frame::Frame, trajectory::Trajectory, unit_cell::UnitCell};
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/xyz/trajectory.xyz");
        let trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 2);

        let path = Path::new("./src/tests-data/xyz/helium.xyz");
        let trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 397);

        let path = Path::new("./src/tests-data/xyz/topology.xyz");
        let trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 1);
    }

    #[test]
    fn extended_xyz() {
        let path = Path::new("./src/tests-data/xyz/extended.xyz");
        let mut trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 3);

        let frame = trajectory.read_at(0).unwrap();
        assert_eq!(frame.size(), 192);

        // Reading the unit cell
        let mut unit_cell = UnitCell::new();
        unit_cell.cell_matrix[(0, 0)] = 8.43116035;
        unit_cell.cell_matrix[(0, 1)] = 0.158219155128;
        unit_cell.cell_matrix[(1, 1)] = 14.5042431863;
        unit_cell.cell_matrix[(0, 2)] = 1.16980663624;
        unit_cell.cell_matrix[(1, 2)] = 4.4685149855;
        unit_cell.cell_matrix[(2, 2)] = 14.9100096405;
        assert_eq!(frame.unit_cell, unit_cell);

        // Frame level properties
        assert_eq!(frame.properties["ENERGY"].expect_double(), -2069.84934116);
        assert_eq!(frame.properties["Natoms"].expect_double(), 192.0);
        assert_eq!(frame.properties["NAME"].expect_string(), "COBHUW");
        let virial = frame.properties["virial"].expect_matrix3x3();
        let true_virial = [
            222.64807906606316,
            21.561432674983322,
            57.08570699995725,
            21.561432674983322,
            247.82551996948678,
            -13.553549356442273,
            57.08570699995725,
            -13.553549356442273,
            156.87636103904026,
        ];
        for i in 0..9 {
            assert_approx_eq!(virial[i], true_virial[i], 1e-12);
        }
        let stress = frame.properties["stress"].expect_matrix3x3();
        let true_stress = [
            -0.0002564107895985581,
            0.0,
            0.0,
            0.0,
            -1.221_130_654_549_481e-5,
            0.0,
            0.0,
            0.0,
            -0.00026623200192425463,
        ];
        for i in 0..9 {
            assert_approx_eq!(stress[i], true_stress[i], 1e-12);
        }
        assert!(frame.properties["IsStrange"].expect_bool());

        // Atom level properties
        let positions = frame.positions()[0];
        assert_approx_eq!(positions[0], 2.33827271799, 1e-9);
        assert_approx_eq!(positions[1], 4.55315540425, 1e-9);
        assert_approx_eq!(positions[2], 11.5841360926, 1e-9);
        assert_approx_eq!(frame[0].properties["CS_0"].expect_double(), 24.10);
        assert_approx_eq!(frame[0].properties["CS_1"].expect_double(), 31.34);
        assert_approx_eq!(frame[51].properties["CS_0"].expect_double(), -73.98);
        assert_approx_eq!(frame[51].properties["CS_1"].expect_double(), -81.85);

        // different types
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 62);
        let vector3d = frame[0].properties["CS"].expect_vector3d();
        assert_approx_eq!(vector3d[0], 198.20, 1e-12);
        assert_approx_eq!(vector3d[1], 202.27, 1e-12);
        assert_approx_eq!(vector3d[2], 202.27, 1e-12);

        // Different syntaxes for bool values
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 8);
        assert!(frame[0].properties["bool"].expect_bool());
        assert!(frame[1].properties["bool"].expect_bool());
        assert!(frame[2].properties["bool"].expect_bool());
        assert!(frame[3].properties["bool"].expect_bool());
        assert!(!frame[4].properties["bool"].expect_bool());
        assert!(!frame[5].properties["bool"].expect_bool());
        assert!(!frame[6].properties["bool"].expect_bool());
        assert!(!frame[7].properties["bool"].expect_bool());

        assert_eq!(frame[0].properties["int"].expect_double(), 33.0);
        assert_eq!(frame[0].properties["strings_0"].expect_string(), "bar");
        assert_eq!(frame[0].properties["strings_1"].expect_string(), "\"test\"");
    }

    #[test]
    fn read_whole_file() {
        let path = Path::new("./src/tests-data/xyz/helium.xyz");
        let mut trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 397);

        let mut frame = Frame::new();
        while let Some(next_frame) = trajectory.read().unwrap() {
            frame = next_frame;
        }
        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], -1.186037, 1e-12);
        assert_approx_eq!(positions[0][1], 11.439334, 1e-12);
        assert_approx_eq!(positions[0][2], 0.529939, 1e-12);

        assert_approx_eq!(positions[124][0], 5.208778, 1e-12);
        assert_approx_eq!(positions[124][1], 12.707273, 1e-12);
        assert_approx_eq!(positions[124][2], 10.940157, 1e-12);
    }
    #[test]
    fn various_files_formatting() {
        let path = Path::new("./src/tests-data/xyz/spaces.xyz");
        let mut trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 1);
        let frame = trajectory.read().unwrap().unwrap();
        let positions = frame.positions();

        assert_approx_eq!(positions[10][0], 0.8336);
        assert_approx_eq!(positions[10][1], 0.3006);
        assert_approx_eq!(positions[10][2], 0.4968);
    }

    macro_rules! assert_read_at_fails {
        ($index:expr) => {
            let path = Path::new(BAD_EXTENDED);
            let mut trajectory = Trajectory::new(path).unwrap();
            trajectory.read_at($index).unwrap();
        };
    }

    const BAD_EXTENDED: &str = "./src/tests-data/xyz/bad_extended.xyz";

    #[test]
    #[should_panic(expected = "MissingToken")]
    fn read_bad_files0() {
        assert_read_at_fails!(0);
    }
    #[test]
    #[should_panic(expected = "Failed to parse number")]
    fn read_bad_files1() {
        assert_read_at_fails!(1);
    }
    #[test]
    #[should_panic(expected = "MissingToken")]
    fn read_bad_files2() {
        assert_read_at_fails!(2);
    }
    #[test]
    #[should_panic(expected = "Failed to parse number")]
    fn read_bad_files3() {
        assert_read_at_fails!(3);
    }
    #[test]
    #[should_panic(expected = "MissingToken")]
    fn read_bad_files4() {
        assert_read_at_fails!(4);
    }
    #[test]
    #[should_panic(expected = "Failed to parse z component")]
    fn read_bad_files5() {
        assert_read_at_fails!(5);
    }
    #[test]
    #[should_panic(expected = "MissingToken")]
    fn read_bad_files6() {
        assert_read_at_fails!(6);
    }
    #[test]
    #[should_panic(expected = "Invalid boolean value: ok")]
    fn read_bad_files7() {
        assert_read_at_fails!(7);
    }
    #[test]
    #[should_panic(expected = "MissingToken")]
    fn read_bad_files8() {
        assert_read_at_fails!(8);
    }
}
