// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use nalgebra::Matrix3;

use crate::atom::Atom;
use crate::error::CError;
use crate::extendedxyzparser::ExtendedXyzParser;
use crate::format::FileFormat;
use crate::frame::Frame;
use crate::property::{Properties, Property, PropertyKind};
use crate::unit_cell::{self, UnitCell};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, Write};
use std::str::SplitWhitespace;

pub struct XYZFormat;
type PropertiesList = BTreeMap<String, PropertyKind>;

impl XYZFormat {
    pub fn parse_unitcell(lattice: &str) -> UnitCell {
        let mut matrix = Matrix3::zeros();

        matrix.iter_mut().zip(lattice.split_whitespace()).for_each(
            |(matrix_entry, lattice_item)| {
                *matrix_entry = fast_float::parse(lattice_item).expect("expected float");
            },
        );

        UnitCell::new_from_matrix(matrix).expect("failed to convert lattice to unit cell")
    }

    fn read_atomic_properties(
        properties: &PropertiesList,
        tokens: &mut SplitWhitespace,
        atom: &mut Atom,
    ) -> Result<(), CError> {
        for (name, kind) in properties {
            if kind == &PropertyKind::Vector3D {
                // For Vector3D, collect all three components first
                let x = tokens.next().ok_or(CError::MissingToken)?;
                let y = tokens.next().ok_or(CError::MissingToken)?;
                let z = tokens.next().ok_or(CError::MissingToken)?;
                let value = format!("{x} {y} {z}");
                let property = Property::parse_value(&value, kind)?;
                atom.properties.insert(name.clone(), property);
            } else {
                let value = tokens.next().ok_or(CError::MissingToken)?;
                let property = Property::parse_value(value, kind)?;
                atom.properties.insert(name.clone(), property);
            }
        }

        Ok(())
    }

    fn parse_property_list(line: &str) -> Result<PropertiesList, CError> {
        const PREFIX: &str = "species:S:1:pos:R:3:";

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
                        "Unknown property type: {unknown}"
                    )));
                }
            };

            let count = chunk[2].parse::<usize>().map_err(|e| {
                CError::GenericError(format!("Invalid property count '{}': {e}", chunk[2]))
            })?;

            if count == 0 {
                return Err(CError::GenericError(format!(
                    "Invalid count of 0 for property '{name}'"
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
            return Property::parse_value(&lowercased, &PropertyKind::Bool)
                .unwrap_or_else(|_| Property::String(value.to_string()));
        }

        // Try parsing as vector/matrix
        let parts: Vec<&str> = value.split_whitespace().collect();
        match parts.len() {
            1 => Property::parse_value(value, &PropertyKind::Double)
                .unwrap_or_else(|_| Property::String(value.to_string())),
            3 => Property::parse_value(value, &PropertyKind::Vector3D)
                .unwrap_or_else(|_| Property::String(value.to_string())),
            9 => Property::parse_value(value, &PropertyKind::Matrix3x3)
                .unwrap_or_else(|_| Property::String(value.to_string())),
            _ => Property::parse_value(value, &PropertyKind::VectorXD)
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
            frame.unit_cell = XYZFormat::parse_unitcell(lattice);
        }

        if let Some(prop_string) = properties.get("Properties") {
            return Self::parse_property_list(prop_string);
        }

        Ok(PropertiesList::new())
    }

    fn is_valid_property_name(name: &str) -> bool {
        let mut chars = name.chars();
        // Must have at least one character and start with an ASCII letter
        match chars.next() {
            Some(first) if first.is_ascii_alphabetic() => {
                chars.all(|c| c.is_ascii_alphanumeric() || c == '_')
            }
            _ => false,
        }
    }

    fn should_be_quoted(s: &str) -> bool {
        // TODO: ASE also allow [] {} and () to function as quotes. This should
        // be updated when a specification is agreed on.
        s.chars()
            .any(|c| c.is_ascii_whitespace() || c == '=' || c == '\'' || c == '"')
    }

    fn get_atom_properties(frame: &Frame) -> PropertiesList {
        if frame.size() == 0 {
            return PropertiesList::new();
        }

        let mut all_properties = HashMap::new();
        let mut partially_defined_already_warned = HashSet::new();

        // Process first atom to establish baseline properties
        let first_atom = &frame[0];
        if !first_atom.properties.is_empty() {
            for (name, property) in &first_atom.properties {
                if !XYZFormat::is_valid_property_name(name) {
                    eprintln!(
                        "warning: '{name}' is not a valid property name for extended XYZ. it will not be saved"
                    );
                    partially_defined_already_warned.insert(name.clone());
                    continue;
                }

                if let Some(string_prop) = property.as_string() {
                    if XYZFormat::should_be_quoted(string_prop) {
                        eprint!(
                            "warning: string value for property '{name}' on atom 0 cannot be saved as an atomic property"
                        );
                        continue;
                    }
                }
                all_properties.insert(name.clone(), property.kind());
            }
        }

        // Check remaining atoms
        for atom in frame.iter().skip(1) {
            // Use retain to filter properties in-place
            all_properties.retain(|prop_name, prop_kind| {
                match atom.properties.get(prop_name) {
                    None => {
                        // Property missing on this atom
                        eprintln!(
                            "warning: property '{prop_name}' is only defined for a subset of atoms. it will not be saved"
                        );
                        false // Remove this property
                    }
                    Some(prop) if prop.kind() != *prop_kind => {
                        // Property has different type on this atom
                        eprintln!(
                            "warning: property '{prop_name}' is defined with different types on different atoms. it will not be saved",
                        );
                        partially_defined_already_warned.insert(prop_name.clone());
                        false // Remove this property
                    }
                    _ => true, // Keep this property
                }
            });

            // Check for properties on this atom that weren't on the first atom
            for (prop_name, _) in &atom.properties {
                if !all_properties.contains_key(prop_name)
                    && !partially_defined_already_warned.contains(prop_name)
                {
                    eprintln!(
                        "warning: property '{prop_name}' is only defined for a subset of atoms. it will not be saved",
                    );
                    partially_defined_already_warned.insert(prop_name.clone());
                }
            }
        }

        // Convert to BTreeMap for result
        let mut results = PropertiesList::new();
        for (name, kind) in all_properties {
            results.insert(name, kind);
        }
        results
    }

    fn write_extended_comment_line(&self, frame: &Frame, properties: &PropertiesList) -> String {
        let mut result = "Properties=species:S:1:pos:R:3".to_string();

        for property in properties {
            let (prop_type, count) = match *property.1 {
                PropertyKind::String => ('S', 1),
                PropertyKind::Bool => ('L', 1),
                PropertyKind::Double => ('R', 1),
                PropertyKind::Vector3D => ('R', 3),
                PropertyKind::Matrix3x3 => ('R', 9),
                PropertyKind::VectorXD => todo!(),
            };
            write!(result, ":{}:{prop_type}:{count}", property.0).unwrap();
        }

        // support for lattice
        if frame.unit_cell.shape != unit_cell::CellShape::Infinite {
            // set small elements to 0
            let m = frame
                .unit_cell
                .matrix
                .map(|m| if m.abs() < 1e-12 { 0.0 } else { m });
            write!(
                result,
                " Lattice=\"{:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}\"",
                m[(0, 0)],
                m[(0, 1)],
                m[(0, 2)],
                m[(1, 0)],
                m[(1, 1)],
                m[(1, 2)],
                m[(2, 0)],
                m[(2, 1)],
                m[(2, 2)]
            )
            .unwrap();
        }

        // sort properties to have reproducible output
        let sorted_properties: BTreeMap<_, _> = frame.properties.iter().collect();

        // support for generic frame properties
        for item in sorted_properties {
            if XYZFormat::should_be_quoted(item.0) {
                // quote the string
                if !item.0.contains('\"') {
                    write!(result, " \"{}\"=", item.0).unwrap();
                } else if !item.0.contains('\'') {
                    write!(result, " '{}'=", item.0).unwrap();
                } else {
                    eprintln!(
                        "warning: frame property '{}' contains both single and double quote. it will not be saved",
                        item.0
                    );
                    continue;
                }
            } else {
                write!(result, " {}=", item.0).unwrap();
            }

            match item.1.kind() {
                PropertyKind::String => write!(result, "\"{}\"", item.1.expect_string()).unwrap(),
                PropertyKind::Bool => result.push_str(if item.1.expect_bool() { "T" } else { "F" }),
                PropertyKind::Double => write!(result, "{:?}", item.1.expect_double()).unwrap(),
                PropertyKind::Vector3D => {
                    let v = item.1.expect_vector3d();
                    write!(result, "\"{:?} {:?} {:?}\"", v[0], v[1], v[2]).unwrap();
                }
                PropertyKind::Matrix3x3 => {
                    let v = item.1.expect_matrix3x3();
                    write!(
                        result,
                        "\"{:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}\"",
                        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]
                    )
                    .unwrap();
                }
                PropertyKind::VectorXD => todo!(),
            }
        }

        result
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

            let symbol = tokens.next().ok_or(CError::UnexpectedSymbol)?.to_string();

            let x: f64 = tokens.next().ok_or(CError::MissingToken)?.parse()?;
            let y: f64 = tokens.next().ok_or(CError::MissingToken)?.parse()?;
            let z: f64 = tokens.next().ok_or(CError::MissingToken)?.parse()?;

            let mut atom = Atom {
                symbol,
                name: String::new(),
                properties: Properties::new(),
            };
            XYZFormat::read_atomic_properties(&properties, &mut tokens, &mut atom)?;

            frame.add_atom(atom, [x, y, z]);
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
                    expected: n_atoms + 2, // first count line + n_atoms atom lines + blank/comment?
                    got: i + 1,            // how many we actually read
                });
            }
        }
        let position = reader.stream_position().unwrap();
        Ok(Some(position))
    }

    fn write_next(&self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
        let positions = frame.positions();
        let properties = XYZFormat::get_atom_properties(frame);

        writeln!(writer, "{}", frame.size())?;
        writeln!(
            writer,
            "{}",
            self.write_extended_comment_line(frame, &properties)
        )?;

        for (atom, pos) in frame.iter().zip(positions) {
            let mut name = atom.name.clone();

            if name.is_empty() {
                name = "X".to_string();
            }
            write!(writer, "{name} {:?} {:?} {:?}", pos[0], pos[1], pos[2])?;

            for property in &properties {
                let val = atom.properties.get(property.0).unwrap();

                if *property.1 == PropertyKind::String {
                    write!(writer, " {}", val.expect_string())?;
                } else if *property.1 == PropertyKind::Bool {
                    if val.expect_bool() {
                        write!(writer, " T")?;
                    } else {
                        write!(writer, " F")?;
                    }
                } else if *property.1 == PropertyKind::Double {
                    write!(writer, " {:?}", val.expect_double())?;
                } else if *property.1 == PropertyKind::Vector3D {
                    let v = val.expect_vector3d();
                    write!(writer, " {:?} {:?} {:?}", v[0], v[1], v[2])?;
                }
            }
            writeln!(writer)?;
        }

        Ok(())
    }

    fn finalize(&self, _writer: &mut BufWriter<File>) -> Result<(), CError> {
        // XYZ format doesn't need any special finalization
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::{
        atom::Atom, frame::Frame, property::Property, trajectory::Trajectory, unit_cell::UnitCell,
    };
    use assert_approx_eq::assert_approx_eq;
    use nalgebra::Matrix3;
    use tempfile::Builder;

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/xyz/single_struct.xyz");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 1);

        let path = Path::new("./src/tests-data/xyz/trajectory.xyz");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 2);

        let path = Path::new("./src/tests-data/xyz/helium.xyz");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 397);

        let path = Path::new("./src/tests-data/xyz/topology.xyz");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 1);
    }

    #[test]
    fn extended_xyz() {
        let path = Path::new("./src/tests-data/xyz/extended.xyz");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 3);

        let frame = trajectory.read_at(0).unwrap().unwrap();
        assert_eq!(frame.size(), 192);

        // Reading the unit cell
        let mut unit_cell = UnitCell::new();
        unit_cell.matrix[(0, 0)] = 8.431_160_35;
        unit_cell.matrix[(0, 1)] = 0.158_219_155_128;
        unit_cell.matrix[(1, 1)] = 14.504_243_186_3;
        unit_cell.matrix[(0, 2)] = 1.169_806_636_24;
        unit_cell.matrix[(1, 2)] = 4.468_514_985_5;
        unit_cell.matrix[(2, 2)] = 14.910_009_640_5;
        assert_eq!(frame.unit_cell, unit_cell);

        // Frame level properties
        assert_approx_eq!(
            frame.properties["ENERGY"].expect_double(),
            -2_069.849_341_16
        );
        assert_approx_eq!(frame.properties["Natoms"].expect_double(), 192.0);
        assert_eq!(frame.properties["NAME"].expect_string(), "COBHUW");
        let virial = frame.properties["virial"].expect_matrix3x3();
        let true_virial = [
            222.648_079_066_063_16,
            21.561_432_674_983_322,
            57.085_706_999_957_25,
            21.561_432_674_983_322,
            247.825_519_969_486_78,
            -13.553_549_356_442_273,
            57.085_706_999_957_25,
            -13.553_549_356_442_273,
            156.876_361_039_040_26,
        ];
        for i in 0..9 {
            assert_approx_eq!(virial[i], true_virial[i], 1e-12);
        }
        let stress = frame.properties["stress"].expect_matrix3x3();
        let true_stress = [
            -0.000_256_410_789_598_558_1,
            0.0,
            0.0,
            0.0,
            -1.221_130_654_549_481e-5,
            0.0,
            0.0,
            0.0,
            -0.000_266_232_001_924_254_63,
        ];
        for i in 0..9 {
            assert_approx_eq!(stress[i], true_stress[i], 1e-12);
        }
        assert!(frame.properties["IsStrange"].expect_bool());

        // Atom level properties
        let positions = frame.positions()[0];
        assert_approx_eq!(positions[0], 2.338_272_717_99, 1e-9);
        assert_approx_eq!(positions[1], 4.553_155_404_25, 1e-9);
        assert_approx_eq!(positions[2], 11.584_136_092_6, 1e-9);
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

        assert_approx_eq!(frame[0].properties["int"].expect_double(), 33.0);
        assert_eq!(frame[0].properties["strings_0"].expect_string(), "bar");
        assert_eq!(frame[0].properties["strings_1"].expect_string(), "\"test\"");
    }

    #[test]
    fn read_whole_file() {
        let path = Path::new("./src/tests-data/xyz/helium.xyz");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 397);

        let mut frame = Frame::new();
        while let Some(next_frame) = trajectory.read().unwrap() {
            frame = next_frame;
        }
        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], -1.186_037, 1e-12);
        assert_approx_eq!(positions[0][1], 11.439_334, 1e-12);
        assert_approx_eq!(positions[0][2], 0.529_939, 1e-12);

        assert_approx_eq!(positions[124][0], 5.208_778, 1e-12);
        assert_approx_eq!(positions[124][1], 12.707_273, 1e-12);
        assert_approx_eq!(positions[124][2], 10.940_157, 1e-12);
    }
    #[test]
    fn various_files_formatting() {
        let path = Path::new("./src/tests-data/xyz/spaces.xyz");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 1);
        let frame = trajectory.read().unwrap().unwrap();
        let positions = frame.positions();

        assert_approx_eq!(positions[10][0], 0.833_6);
        assert_approx_eq!(positions[10][1], 0.300_6);
        assert_approx_eq!(positions[10][2], 0.496_8);
    }

    macro_rules! assert_read_at_fails {
        ($index:expr) => {
            let path = Path::new(BAD_EXTENDED);
            let mut trajectory = Trajectory::open(&path).unwrap();
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

    const EXPECTED_CONTENT: &str = r#"4
Properties=species:S:1:pos:R:3:bool:L:1:double:R:1:string:S:1:vector:R:3 name="Test" stress="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"
A 1.0 2.0 3.0 T 10.0 atom_0 10.0 20.0 30.0
B 1.0 2.0 3.0 F 11.0 atom_1 11.0 21.0 31.0
C 1.0 2.0 3.0 T 12.0 atom_2 12.0 22.0 32.0
D 1.0 2.0 3.0 T 13.0 atom_2 13.0 23.0 33.0
6
Properties=species:S:1:pos:R:3 Lattice="12.0 0.0 0.0 0.0 13.0 0.0 0.0 0.0 14.0" direction="1.0 0.0 2.0" is_open=F name="Test" 'quotes"'=T "quotes'"=T speed=33.4 stress="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" "with space"=T
A 1.0 2.0 3.0
B 1.0 2.0 3.0
C 1.0 2.0 3.0
D 1.0 2.0 3.0
E 4.0 5.0 6.0
F 4.0 5.0 6.0
"#;
    #[test]
    fn test_xyz_file_contents() {
        let named_tmpfile = Builder::new()
            .prefix("temporary-xyz")
            .suffix(".xyz")
            .tempfile()
            .unwrap();

        // Write the expected content into the temp file
        let mut trajectory = Trajectory::create(named_tmpfile.path()).unwrap();
        let mut frame = Frame::new();
        frame
            .properties
            .insert("name".to_string(), Property::String("Test".to_string()));
        frame
            .properties
            .insert("stress".to_string(), Property::Matrix3x3(Matrix3::zeros()));
        let atom = Atom::with_symbol("A".to_string(), "O".to_string());
        frame.add_atom(atom, [1.0, 2.0, 3.0]);
        let atom = Atom::new("B".to_string());
        frame.add_atom(atom, [1.0, 2.0, 3.0]);
        let atom = Atom::new("C".to_string());
        frame.add_atom(atom, [1.0, 2.0, 3.0]);
        let atom = Atom::new("D".to_string());
        frame.add_atom(atom, [1.0, 2.0, 3.0]);

        // atomic properties
        frame[0]
            .properties
            .insert("string".to_string(), Property::String("atom_0".to_string()));
        frame[1]
            .properties
            .insert("string".to_string(), Property::String("atom_1".to_string()));
        frame[2]
            .properties
            .insert("string".to_string(), Property::String("atom_2".to_string()));
        frame[3]
            .properties
            .insert("string".to_string(), Property::String("atom_2".to_string()));

        frame[0]
            .properties
            .insert("bool".to_string(), Property::Bool(true));
        frame[1]
            .properties
            .insert("bool".to_string(), Property::Bool(false));
        frame[2]
            .properties
            .insert("bool".to_string(), Property::Bool(true));
        frame[3]
            .properties
            .insert("bool".to_string(), Property::Bool(true));

        frame[0]
            .properties
            .insert("double".to_string(), Property::Double(10.0));
        frame[1]
            .properties
            .insert("double".to_string(), Property::Double(11.0));
        frame[2]
            .properties
            .insert("double".to_string(), Property::Double(12.0));
        frame[3]
            .properties
            .insert("double".to_string(), Property::Double(13.0));

        frame[0]
            .properties
            .insert("vector".to_string(), Property::Vector3D([10.0, 20.0, 30.0]));
        frame[1]
            .properties
            .insert("vector".to_string(), Property::Vector3D([11.0, 21.0, 31.0]));
        frame[2]
            .properties
            .insert("vector".to_string(), Property::Vector3D([12.0, 22.0, 32.0]));
        frame[3]
            .properties
            .insert("vector".to_string(), Property::Vector3D([13.0, 23.0, 33.0]));

        // not saved, bad property name
        frame[0]
            .properties
            .insert("value with spaces".to_string(), Property::Double(0.0));
        frame[1]
            .properties
            .insert("value with spaces".to_string(), Property::Double(0.0));
        frame[2]
            .properties
            .insert("value with spaces".to_string(), Property::Double(0.0));
        frame[3]
            .properties
            .insert("value with spaces".to_string(), Property::Double(0.0));

        // not saved, different types
        frame[0]
            .properties
            .insert("value".to_string(), Property::Double(0.0));
        frame[1]
            .properties
            .insert("value".to_string(), Property::String("0.0".to_string()));
        frame[2]
            .properties
            .insert("value".to_string(), Property::Bool(false));
        frame[3]
            .properties
            .insert("value".to_string(), Property::Double(0.0));

        trajectory.write(&frame).unwrap();

        frame.unit_cell = UnitCell::new_from_lengths([12.0, 13.0, 14.0]);
        frame
            .properties
            .insert("is_open".to_string(), Property::Bool(false));
        frame
            .properties
            .insert("speed".to_string(), Property::Double(33.4));
        frame
            .properties
            .insert("direction".to_string(), Property::Vector3D([1.0, 0.0, 2.0]));
        frame
            .properties
            .insert("with space".to_string(), Property::Bool(true));
        frame
            .properties
            .insert("quotes'".to_string(), Property::Bool(true));
        frame
            .properties
            .insert("quotes\"".to_string(), Property::Bool(true));

        // properties with two types of quotes are skipped
        frame
            .properties
            .insert("all_quotes'\"".to_string(), Property::Bool(true));

        let atom = Atom::new("E".to_string());
        frame.add_atom(atom, [4.0, 5.0, 6.0]);

        let atom = Atom::new("F".to_string());
        frame.add_atom(atom, [4.0, 5.0, 6.0]);

        trajectory.write(&frame).unwrap();

        let contents = std::fs::read_to_string(named_tmpfile.path()).unwrap();
        assert_eq!(contents, EXPECTED_CONTENT);
    }

    #[test]
    fn triclinic_cell_with_negative_values() {
        let mut matrix = Matrix3::zeros();
        matrix[(0, 0)] = 6.923_95;
        matrix[(1, 0)] = -3.224_55;
        matrix[(2, 0)] = 0.000_0;
        matrix[(0, 1)] = 0.000_00;
        matrix[(1, 1)] = 5.453_55;
        matrix[(2, 1)] = 0.000_0;
        matrix[(0, 2)] = 0.100_667;
        matrix[(1, 2)] = -3.320_57;
        matrix[(2, 2)] = 7.283_6;
        let mut frame = Frame::new();
        frame.unit_cell = UnitCell::new_from_matrix(matrix).unwrap();
        frame.resize(1).unwrap();

        let named_tmpfile = Builder::new()
            .prefix("triclinic-xyz")
            .suffix(".xyz")
            .tempfile()
            .unwrap();

        // Write the expected content into the temp file
        let mut trajectory = Trajectory::create(named_tmpfile.path()).unwrap();
        trajectory.write(&frame).unwrap();

        let expected = r#"1
Properties=species:S:1:pos:R:3 Lattice="6.92395 0.0 0.100667 -3.22455 5.45355 -3.32057 0.0 0.0 7.2836"
X 0.0 0.0 0.0
"#;
        let contents = std::fs::read_to_string(named_tmpfile.path()).unwrap();
        assert_eq!(contents, expected);
    }
}
