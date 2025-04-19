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

    // fn read(&self) -> Result<Frame, CError> {
    //     println!("Reading as XYZ format");
    //     Ok(Frame { atoms: vec![] })
    // }

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

    use crate::{trajectory::Trajectory, unit_cell::UnitCell};
    use assert_approx_eq::assert_approx_eq;

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

        assert_eq!(frame.atoms[0].symbol, "O");
        assert_eq!(frame.atoms[1].symbol, "O");

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
