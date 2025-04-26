use nalgebra::Matrix3;

use crate::atom::Atom;
use crate::error::CError;
use crate::extendedxyzparser::ExtendedXyzParser;
use crate::format::FileFormat;
use crate::frame::Frame;
use crate::property::{self, Properties};
use crate::property::{Property, PropertyKind};
use crate::unit_cell::{self, UnitCell};
use std::borrow::BorrowMut;
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
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
            match kind {
                PropertyKind::Vector3D => {
                    // For Vector3D, collect all three components first
                    let x = tokens.next().ok_or(CError::MissingToken)?;
                    let y = tokens.next().ok_or(CError::MissingToken)?;
                    let z = tokens.next().ok_or(CError::MissingToken)?;
                    let value = format!("{} {} {}", x, y, z);
                    let property = Property::parse_value(&value, kind)?;
                    atom.properties.insert(name.clone(), property);
                }
                _ => {
                    let value = tokens.next().ok_or(CError::MissingToken)?;
                    let property = Property::parse_value(value, kind)?;
                    atom.properties.insert(name.clone(), property);
                }
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
        if name.is_empty() {
            return false;
        }

        return name.chars().all(|c| !c.is_ascii_alphanumeric() && c != '_');
    }

    fn should_be_quoted(s: &str) -> bool {
        // TODO: ASE also allow [] {} and () to function as quotes. This should
        // be updated when a specification is agreed on.
        return s
            .chars()
            .any(|c| c.is_ascii_whitespace() || c == '=' || c == '\\' || c == '"');
    }

    fn get_atom_properties(&self, frame: &Frame) -> PropertiesList {
        if frame.size() == 0 {
            return PropertiesList::new();
        };

        let mut all_properties = RefCell::new(HashMap::new());
        let mut partially_defined_already_warned = HashSet::new();
        let first_atom = frame[0].clone();
        if !first_atom.properties.is_empty() {
            for property in first_atom.properties {
                if !XYZFormat::is_valid_property_name(&property.0) {
                    eprintln!(
                        "warning: '{}' is not a valid property name for extended XYZ. it will not be saved",
                        property.0
                    );
                    partially_defined_already_warned.insert(property.0.clone());
                    continue;
                }

                if let Some(string_prop) = property.1.as_string() {
                    if XYZFormat::should_be_quoted(&string_prop) {
                        eprint!(
                            "warning: string value for property '{}' on atom 0 cannot be saved as an atomic property",
                            property.0
                        );
                        continue;
                    }
                }
                all_properties
                    .borrow_mut()
                    .insert(property.0, property.1.kind());
            }
        }

        for atom in frame.iter() {
            let mut to_remove = vec![];
            let properties = all_properties.borrow();
            for property in properties.iter() {
                let current_property = atom.properties.get(property.0);
                if current_property.is_none() {
                    // this property was present for all atoms until now, but not on the current one
                    eprintln!(
                        "warning: property '{}' wis only defined for a subset of atoms. it will not be saved",
                        property.0
                    );

                    to_remove.push(property.0);
                    continue;
                }

                if current_property.unwrap().kind() != *property.1 {
                    eprintln!(
                        "warning: property '{}' is defined with different types on different atoms. it will not be saved",
                        property.0
                    );

                    partially_defined_already_warned.insert(property.0.clone());
                    to_remove.push(property.0);
                }
            }

            let atom_properties = &atom.properties;
            if !atom_properties.is_empty()
                && (atom_properties.len() > all_properties.borrow().len())
            {
                // warn for properties defined on this atom but not on others
                for property in &atom.properties {
                    if !all_properties.borrow().contains_key(property.0)
                        && !partially_defined_already_warned.contains(property.0)
                    {
                        eprintln!(
                            "warning: property '{}' is only defined for a subset of atoms. it will not be saved",
                            property.0
                        );
                        partially_defined_already_warned.insert(property.0.clone());
                    }
                }
            }

            for name in to_remove {
                all_properties.borrow_mut().remove(name);
            }
        }

        let mut results = PropertiesList::new();
        for property in all_properties.borrow().iter() {
            results.insert(property.0.clone(), property.1.clone());
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
            result.push_str(&format!("{}:{prop_type}:{count}", property.0));
        }

        // support for lattice

        if frame.unit_cell.shape != unit_cell::CellShape::Infinite {
            // set small elements to 0
            let m = frame
                .unit_cell
                .matrix
                .map(|m| if m.abs() < 1e-12 { 0.0 } else { m });
            result.push_str(&format!(
                " Lattice=\"{:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e} {:e}\"",
                m[(0, 0)],
                m[(0, 1)],
                m[(0, 2)],
                m[(1, 0)],
                m[(1, 1)],
                m[(1, 2)],
                m[(2, 0)],
                m[(2, 1)],
                m[(2, 2)]
            ));
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
                name: "".to_string(),
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

        let mut line_position = 0;
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

    fn write_next(&self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
        let positions = frame.positions();
        let properties = self.get_atom_properties(frame);

        write!(writer, "{}", frame.size())?;
        write!(
            writer,
            "{}",
            self.write_extended_comment_line(frame, &properties)
        )?;
        dbg!(properties);

        Ok(())
        //         const auto& positions = frame.positions();
        // auto properties = get_atom_properties(frame);

        // file_.print("{}\n", frame.size());
        // file_.print("{}\n", write_extended_comment_line(frame, properties));

        // for (size_t i = 0; i < frame.size(); i++) {
        //     const auto& atom = frame[i];

        //     auto name = atom.name();
        //     if (name.empty()) {
        //         name = "X";
        //     }

        //     file_.print("{} {:g} {:g} {:g}",
        //         name, positions[i][0], positions[i][1], positions[i][2]
        //     );

        //     for (const auto& property: properties) {
        //         const auto& value = atom.get(property.name).value();

        //         if (property.type == Property::STRING) {
        //             file_.print(" {}", value.as_string());
        //         } else if (property.type == Property::BOOL) {
        //             if (value.as_bool()) {
        //                 file_.print(" T");
        //             } else {
        //                 file_.print(" F");
        //             }
        //         } else if (property.type == Property::DOUBLE) {
        //             file_.print(" {:g}", value.as_double());
        //         } else if (property.type == Property::VECTOR3D) {
        //             const auto& vector = value.as_vector3d();
        //             file_.print(" {:g} {:g} {:g}", vector[0], vector[1], vector[2]);
        //         }
        //     }

        //     file_.print("\n");
        // }
    }

    // fn write(&self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
    //     self.write_next(writer, frame);
    //     self.frame
    //     println!("Writingmut  as XYZ format with {} atoms", frame.size());
    //     Ok(())
    // }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::{
        frame::Frame,
        property::Property,
        trajectory::{FileMode, Trajectory},
        unit_cell::UnitCell,
    };
    use assert_approx_eq::assert_approx_eq;
    use tempfile::NamedTempFile;

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/xyz/single_struct.xyz");
        let trajectory = Trajectory::new(path, FileMode::Read).unwrap();
        assert_eq!(trajectory.size, 1);

        let path = Path::new("./src/tests-data/xyz/trajectory.xyz");
        let trajectory = Trajectory::new(path, FileMode::Read).unwrap();
        assert_eq!(trajectory.size, 2);

        let path = Path::new("./src/tests-data/xyz/helium.xyz");
        let trajectory = Trajectory::new(path, FileMode::Read).unwrap();
        assert_eq!(trajectory.size, 397);

        let path = Path::new("./src/tests-data/xyz/topology.xyz");
        let trajectory = Trajectory::new(path, FileMode::Read).unwrap();
        assert_eq!(trajectory.size, 1);
    }

    #[test]
    fn extended_xyz() {
        let path = Path::new("./src/tests-data/xyz/extended.xyz");
        let mut trajectory = Trajectory::new(path, FileMode::Read).unwrap();
        assert_eq!(trajectory.size, 3);

        let frame = trajectory.read_at(0).unwrap();
        assert_eq!(frame.size(), 192);

        // Reading the unit cell
        let mut unit_cell = UnitCell::new();
        unit_cell.matrix[(0, 0)] = 8.43116035;
        unit_cell.matrix[(0, 1)] = 0.158219155128;
        unit_cell.matrix[(1, 1)] = 14.5042431863;
        unit_cell.matrix[(0, 2)] = 1.16980663624;
        unit_cell.matrix[(1, 2)] = 4.4685149855;
        unit_cell.matrix[(2, 2)] = 14.9100096405;
        assert_eq!(frame.unit_cell, unit_cell);

        // Frame level properties
        assert_approx_eq!(frame.properties["ENERGY"].expect_double(), -2069.84934116);
        assert_approx_eq!(frame.properties["Natoms"].expect_double(), 192.0);
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

        assert_approx_eq!(frame[0].properties["int"].expect_double(), 33.0);
        assert_eq!(frame[0].properties["strings_0"].expect_string(), "bar");
        assert_eq!(frame[0].properties["strings_1"].expect_string(), "\"test\"");
    }

    #[test]
    fn read_whole_file() {
        let path = Path::new("./src/tests-data/xyz/helium.xyz");
        let mut trajectory = Trajectory::new(path, FileMode::Read).unwrap();
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
        let mut trajectory = Trajectory::new(path, FileMode::Read).unwrap();
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
            let mut trajectory = Trajectory::new(path, FileMode::Read).unwrap();
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

    #[test]
    fn test_xyz_file_contents() -> std::io::Result<()> {
        let mut tmpfile = NamedTempFile::new()?;
        const EXPECTED_CONTENT: &str = r#"4
Properties=species:S:1:pos:R:3:bool:L:1:double:R:1:string:S:1:vector:R:3 name="Test"
A 1 2 3 T 10 atom_0 10 20 30
B 1 2 3 F 11 atom_1 11 21 31
C 1 2 3 T 12 atom_2 12 22 32
D 1 2 3 T 13 atom_2 13 23 33
6
Properties=species:S:1:pos:R:3 Lattice="12 0 0 0 13 0 0 0 14" direction="1 0 2" is_open=F name="Test" 'quotes"'=T "quotes'"=T speed=33.4 "with space"=T
A 1 2 3
B 1 2 3
C 1 2 3
D 1 2 3
E 4 5 6
F 4 5 6
"#;

        // Write the expected content into the temp file
        let mut trajectory = Trajectory::new(tmpfile.path(), FileMode::Write).unwrap();
        let mut frame = Frame::new();
        frame
            .properties
            .insert("name".to_string(), Property::String("Test".to_string()));
        // frame.add_atom(Atom("A","O"), {1, 2, 3});
        // frame.add_atom(Atom("B"), {1, 2, 3});
        // frame.add_atom(Atom("C"), {1, 2, 3});
        // frame.add_atom(Atom("D"), {1, 2, 3});

        trajectory.write(&frame);
        // tmpfile.write_all(EXPECTED_CONTENT.as_bytes())?;

        // …now you can reopen/read the file, feed it into your parser, assert on the results, etc.

        Ok(())
    }

    // TEST_CASE("Write files in XYZ format") {
    //     auto tmpfile = NamedTempPath(".xyz");
    //     const auto* EXPECTED_CONTENT =
    // R"(4
    // Properties=species:S:1:pos:R:3:bool:L:1:double:R:1:string:S:1:vector:R:3 name="Test"
    // A 1 2 3 T 10 atom_0 10 20 30
    // B 1 2 3 F 11 atom_1 11 21 31
    // C 1 2 3 T 12 atom_2 12 22 32
    // D 1 2 3 T 13 atom_2 13 23 33
    // 6
    // Properties=species:S:1:pos:R:3 Lattice="12 0 0 0 13 0 0 0 14" direction="1 0 2" is_open=F name="Test" 'quotes"'=T "quotes'"=T speed=33.4 "with space"=T
    // A 1 2 3
    // B 1 2 3
    // C 1 2 3
    // D 1 2 3
    // E 4 5 6
    // F 4 5 6
    // )";

    //     auto frame = Frame();
    //     frame.set("name", "Test");
    //     frame.add_atom(Atom("A","O"), {1, 2, 3});
    //     frame.add_atom(Atom("B"), {1, 2, 3});
    //     frame.add_atom(Atom("C"), {1, 2, 3});
    //     frame.add_atom(Atom("D"), {1, 2, 3});

    //     // atomic properties
    //     frame[0].set("string", "atom_0");
    //     frame[1].set("string", "atom_1");
    //     frame[2].set("string", "atom_2");
    //     frame[3].set("string", "atom_2");

    //     frame[0].set("bool", true);
    //     frame[1].set("bool", false);
    //     frame[2].set("bool", true);
    //     frame[3].set("bool", true);

    //     frame[0].set("double", 10);
    //     frame[1].set("double", 11);
    //     frame[2].set("double", 12);
    //     frame[3].set("double", 13);

    //     frame[0].set("vector", Vector3D{10, 20, 30});
    //     frame[1].set("vector", Vector3D{11, 21, 31});
    //     frame[2].set("vector", Vector3D{12, 22, 32});
    //     frame[3].set("vector", Vector3D{13, 23, 33});

    //     // not saved, bad property name
    //     frame[0].set("value with spaces", 0);
    //     frame[1].set("value with spaces", 0);
    //     frame[2].set("value with spaces", 0);
    //     frame[3].set("value with spaces", 0);

    //     // not saved, different types
    //     frame[0].set("value", 0);
    //     frame[1].set("value", "0");
    //     frame[2].set("value", false);
    //     frame[3].set("value", 0);

    //     auto file = Trajectory(tmpfile, 'w');
    //     file.write(frame);

    //     frame.set_cell(UnitCell({12, 13, 14}));
    //     frame.set("is_open", false);
    //     frame.set("speed", 33.4);
    //     frame.set("direction", Vector3D(1, 0, 2));
    //     frame.set("with space", true);
    //     frame.set("quotes'", true);
    //     frame.set("quotes\"", true);

    //     // properties with two type of quotes are skipped
    //     frame.set("all_quotes'\"", true);

    //     frame.add_atom(Atom("E"), {4, 5, 6});
    //     frame.add_atom(Atom("F"), {4, 5, 6});

    //     file.write(frame);
    //     file.close();

    //     auto content = read_text_file(tmpfile);
    //     CHECK(content == EXPECTED_CONTENT);
    // }
}
