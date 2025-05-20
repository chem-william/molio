// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.
//

use std::fmt::Write;
use std::io::Write as _;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Seek},
};

use crate::property::PropertyKind;
use crate::{atom::Atom, bond::BondOrder, property::Property};
use crate::{error::CError, format::FileFormat, frame::Frame};
use log::warn;

/// Properties are indiscriminately encoded as strings. This means that you should always
/// [`Property::expect_string()`] when reading properties from Properties
pub struct SDFFormat;

impl FileFormat for SDFFormat {
    fn read_next(&mut self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        let mut line = String::new();
        let mut frame = Frame::new();
        let _ = reader.read_line(&mut line)?;

        if !line.is_empty() {
            frame
                .properties
                .insert("name".to_string(), Property::String(line.to_string()));
        }
        reader.read_line(&mut line)?; // Program line - skip it
        reader.read_line(&mut line)?; // Comment line - skip it

        line.clear();
        reader.read_line(&mut line)?;
        let natoms = line[..3].trim().parse::<usize>().map_err(|_| {
            CError::GenericError(format!(
                "could not parse atom count in SDF file: '{}'",
                &line[..3]
            ))
        })?;

        let nbonds = line[3..6].trim().parse::<usize>().map_err(|_| {
            CError::GenericError(format!(
                "could not parse bond count in SDF file: '{}'",
                &line[3..6]
            ))
        })?;

        frame.reserve(natoms);
        for _ in 0..natoms {
            line.clear();
            let _ = reader.read_line(&mut line)?;

            if line.len() < 34 {
                return Err(CError::GenericError(format!(
                    "atom line is too small for SDF: '{line}'"
                )));
            }

            let x = line[..10].trim().parse::<f64>()?;
            let y = line[10..20].trim().parse::<f64>()?;
            let z = line[20..30].trim().parse::<f64>()?;
            let name = line[31..34].trim();
            let mut atom = Atom::new(name.to_string());

            if line.len() >= 40 {
                let charge_code = line[36..39].trim().parse::<i64>().unwrap_or_else(|e| {
                    warn!("charge code is not numeric '{}': {e}", &line[36..39]);
                    0
                });
                match charge_code {
                    0 => {}
                    1 => atom.charge = 3.0,
                    2 => atom.charge = 2.0,
                    3 => atom.charge = 1.0,
                    5 => atom.charge = -1.0,
                    6 => atom.charge = -2.0,
                    7 => atom.charge = -3.0,
                    _ => warn!("unknown charge code: {charge_code}"),
                }
            }
            frame.add_atom(atom, [x, y, z]);
        }

        for _ in 0..nbonds {
            line.clear();
            let _ = reader.read_line(&mut line)?;
            let atom1 = line[..3].trim().parse::<usize>()?;
            let atom2 = line[3..6].trim().parse::<usize>()?;
            let order = line[6..9].trim().parse::<usize>()?;

            let bond_order = match order {
                1 => BondOrder::Single,
                2 => BondOrder::Double,
                3 => BondOrder::Triple,
                4 => BondOrder::Aromatic,

                // 8 specifically means unspecified
                8 => BondOrder::Unknown,
                _ => BondOrder::Unknown,
            };
            frame.add_bond(atom1 - 1, atom2 - 1, bond_order)?;
        }
        // Parsing the file is more or less complete now, but atom properties can
        // still be read (until 'M  END' is reached).
        // This loop breaks when the property block ends or returns on an error

        line.clear();
        while reader.read_line(&mut line)? > 0 {
            if line.is_empty() {
                continue;
            } else if line.starts_with("$$$$") {
                // Ending block, technically wrong - but we can exit safely
                return Ok(frame);
            } else if line.starts_with("M  END") {
                // Proper end of block
                break;
            } // TODO: add actual ATOM property parsing here

            line.clear();
        }
        // This portion of the file is for molecule wide properties.
        // We're done parsing, so just quit if any errors occur
        let mut property_name = String::new();
        let mut property_value = String::new();
        line.clear();
        while reader.read_line(&mut line)? > 0 {
            if line.trim().is_empty() {
                // This breaks a property group - so store now
                if property_name.is_empty() {
                    warn!("missing property name");
                }
                frame.properties.insert(
                    property_name.clone(),
                    Property::String(property_value.clone()),
                );
                line.clear();
            } else if line.starts_with("$$$$") {
                // Molecule ending block
                return Ok(frame);
            } else if line.starts_with("> <") {
                // Get the property name which is formatted like:
                // > <NAMEGOESHERE>
                if let Some(npos) = line.rfind('>') {
                    property_name = line[3..npos].to_string();
                    line.clear();
                    reader.read_line(&mut line)?;
                    property_value = line.trim().to_string();
                } else {
                    warn!("no closing delimiter found for property name: {line}");
                }
            } else {
                // Continuation of a property value
                writeln!(property_value).expect("write to string we control(?)");
                write!(property_value, "{line}").expect("write to string we control(?)");
            }

            line.clear();
        }

        Ok(frame)
    }

    fn read(&mut self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        // TODO: replace with has_data_left when stabilized
        if reader.fill_buf().map(|b| !b.is_empty()).unwrap() {
            Ok(Some(self.read_next(reader).unwrap()))
        } else {
            Ok(None)
        }
    }

    fn write_next(&mut self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
        let topology = frame.topology();
        debug_assert!(frame.size() == topology.size());

        let mut frame_name = frame
            .properties
            .get("name")
            .map(|p| p.expect_string())
            .unwrap_or("");
        if frame_name.len() > 80 {
            warn!("the frame 'name' property is too long for the SDF format. It has been truncated to 80 characters");
            frame_name = &frame_name[..80];
        }
        writeln!(writer, "{frame_name}")?;
        // // TODO: this line can contain more data (file creation time and energy in particular)
        writeln!(writer)?;
        writeln!(writer, "created by molio")?;
        writeln!(
            writer,
            "{:>3}{:>3}  0     0  0  0  0  0  0999 V2000",
            frame.size(),
            topology.bonds().len()
        )?;

        for (atom, pos) in topology.atoms.iter().zip(frame.positions().iter()) {
            let mut symbol = atom.symbol.clone();

            if symbol.is_empty() || symbol.len() > 3 {
                symbol = "Xxx".to_string();
            }

            let charge_int = atom.charge.trunc() as i64;
            let charge_code = match charge_int {
                0 => 0, // we don't want to warn so handled explicitly
                1 => 3,
                2 => 2,
                3 => 1,
                -1 => 5,
                -2 => 6,
                -3 => 7,
                _ => {
                    warn!("charge code not available for '{charge_int}'");
                    0
                }
            };

            writeln!(
                writer,
                "{:>10.4}{:>10.4}{:>10.4} {:3} 0{:3}  0  0  0  0  0  0  0  0  0  0",
                pos[0],
                pos[1],
                pos[2],
                symbol.to_string(),
                charge_code
            )?;
        }

        for bond in topology.bonds() {
            let bo = topology
                .bond_order(bond[0], bond[1])
                .unwrap_or_else(|_| panic!("a bond between {} and {}", bond[0], bond[1]));

            let bond_order = match bo {
                BondOrder::Single => "  1".to_string(),
                BondOrder::Double => "  2".to_string(),
                BondOrder::Triple => "  3".to_string(),
                BondOrder::Aromatic => "  4".to_string(),
                BondOrder::Unknown => "  8".to_string(),
                _ => "  8".to_string(),
            };

            writeln!(
                writer,
                "{:>3}{:>3}{}  0  0  0  0",
                bond[0] + 1,
                bond[1] + 1,
                bond_order
            )?;
        }

        writeln!(writer, "M  END")?;

        for prop in &frame.properties {
            if prop.0 == "name" {
                continue;
            }

            writeln!(writer, "> <{}>", prop.0)?;

            match prop.1.kind() {
                PropertyKind::String => writeln!(writer, "{}\n", prop.1.expect_string())?,
                PropertyKind::Double => writeln!(writer, "{:?}\n", prop.1.expect_double())?,
                PropertyKind::Bool => writeln!(writer, "{}\n", prop.1.expect_bool())?,
                PropertyKind::Vector3D => {
                    let v = prop.1.expect_vector3d();
                    writeln!(writer, "{:?} {:?} {:?}\n", v[0], v[1], v[2])?;
                }
                _ => unreachable!("unexpected property"),
            }
        }

        writeln!(writer, "$$$$")?;

        Ok(())
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        let mut buf = Vec::with_capacity(256); // Pre-allocate a reasonable size

        // Ignore header lines: molecule name, metadata, and general comment line
        for _ in 0..3 {
            buf.clear();
            let bytes_read = reader.read_until(b'\n', &mut buf)?;
            if bytes_read == 0 {
                // EOF reached
                return Ok(None);
            }
        }

        // Read counts line
        buf.clear();
        let bytes_read = reader.read_until(b'\n', &mut buf)?;
        if bytes_read == 0 {
            return Ok(None);
        }
        if bytes_read < 10 {
            return Err(CError::GenericError(format!(
                "counts line must have at least 10 characters in SDF file. It has {bytes_read} bytes: '{}'", String::from_utf8_lossy(&buf)
            )));
        }

        // Parse atom and bond counts
        let counts_line = std::str::from_utf8(&buf)
            .map_err(|_| CError::GenericError("invalid UTF-8 in counts line".to_string()))?;

        let natoms = counts_line[..3].trim().parse::<usize>().map_err(|_| {
            CError::GenericError(format!(
                "could not parse atom count in SDF file: '{}'",
                &counts_line[..3]
            ))
        })?;
        let nbonds = counts_line[3..6].trim().parse::<usize>().map_err(|_| {
            CError::GenericError(format!(
                "could not parse bond count in SDF file: '{}'",
                &counts_line[3..6]
            ))
        })?;

        for _ in 0..(natoms + nbonds) {
            buf.clear();
            if reader.read_until(b'\n', &mut buf)? == 0 {
                return Err(CError::GenericError(
                    "unexpected EOF in SDF format".to_string(),
                ));
            }
        }

        // Search for ending character, updating the cursor in the file for the next
        // call to forward
        loop {
            buf.clear();
            let bytes_read = reader.read_until(b'\n', &mut buf)?;
            if bytes_read == 0 {
                break; // EOF reached
            }

            // Check if the line starts with "$$$$"
            if buf.starts_with(b"$$$$") {
                break;
            }
        }

        // We have enough data to parse an entire molecule.
        // So, even if the file may not have an ending string,
        // return the start of this step
        let position = reader.stream_position()?;
        Ok(Some(position))
    }

    fn finalize(&self, _writer: &mut BufWriter<File>) -> Result<(), CError> {
        // SDF format doesn't need any special finalization
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use assert_approx_eq::assert_approx_eq;
    use tempfile::Builder;

    use crate::{
        atom::Atom, bond::BondOrder, frame::Frame, property::Property, trajectory::Trajectory,
    };

    #[test]
    fn check_nsteps_aspirin() {
        let path = Path::new("./src/tests-data/sdf/aspirin.sdf");
        let trajectory = Trajectory::open(path).unwrap();

        assert_eq!(trajectory.size, 1);
    }

    #[test]
    fn check_nsteps_kinases() {
        let path = Path::new("./src/tests-data/sdf/kinases.sdf");
        let trajectory = Trajectory::open(path).unwrap();

        assert_eq!(trajectory.size, 6);
    }

    #[test]
    fn read_next_step() {
        let path = Path::new("./src/tests-data/sdf/kinases.sdf");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 47);

        // Check positions
        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 4.9955, 1e-3);
        assert_approx_eq!(positions[0][1], -2.6277, 1e-3);
        assert_approx_eq!(positions[0][2], 0.2047, 1e-3);

        assert_approx_eq!(positions[46][0], -8.5180, 1e-3);
        assert_approx_eq!(positions[46][1], 0.2962, 1e-3);
        assert_approx_eq!(positions[46][2], 2.1406, 1e-3);

        // Check topology
        let topology = frame.topology();
        assert_eq!(topology.size(), 47);
        assert_eq!(topology[0], Atom::new("O".to_string()));
    }

    #[test]
    fn read_specific_step() {
        let path = Path::new("./src/tests-data/sdf/kinases.sdf");
        let mut trajectory = Trajectory::open(path).unwrap();
        let mut frame = trajectory.read_at(3).unwrap().unwrap();

        let mut positions = frame.positions();
        assert_approx_eq!(positions[0][0], -0.8276, 1e-3);
        assert_approx_eq!(positions[0][1], 0.2486, 1e-3);
        assert_approx_eq!(positions[0][2], -1.0418, 1e-3);

        assert_approx_eq!(positions[67][0], -1.1356, 1e-3);
        assert_approx_eq!(positions[67][1], 5.2260, 1e-3);
        assert_approx_eq!(positions[67][2], 1.3726, 1e-3);

        let topology = frame.topology();
        assert_eq!(topology.size(), 68);
        assert_eq!(topology[0], Atom::new("O".to_string()));

        frame = trajectory.read_at(0).unwrap().unwrap();
        positions = frame.positions();
        assert_approx_eq!(positions[0][0], 4.9955);
        assert_approx_eq!(positions[0][1], -2.6277);
        assert_approx_eq!(positions[0][2], 0.2047);

        assert_approx_eq!(positions[46][0], -8.5180);
        assert_approx_eq!(positions[46][1], 0.2962);
        assert_approx_eq!(positions[46][2], 2.1406);
    }
    #[test]
    fn read_whole_file() {
        let path = Path::new("./src/tests-data/sdf/kinases.sdf");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 6);

        let mut frame = Frame::new();
        while let Some(next_frame) = trajectory.read().unwrap() {
            frame = next_frame;
        }

        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 3.1149);
        assert_approx_eq!(positions[0][1], -1.1207);
        assert_approx_eq!(positions[0][2], 3.0606);

        assert_approx_eq!(positions[49][0], -7.4890);
        assert_approx_eq!(positions[49][1], -0.0147);
        assert_approx_eq!(positions[49][2], -2.1114);
    }

    #[test]
    fn read_various_file_properties() {
        let path = Path::new("./src/tests-data/sdf/aspirin.sdf");
        let mut trajectory = Trajectory::open(path).unwrap();

        let frame = trajectory.read().unwrap().unwrap();
        let prop = frame.properties.get("PUBCHEM_COMPOUND_CID").unwrap();
        assert_eq!(prop.expect_string(), "2244".to_string());

        let prop2 = frame.properties.get("PUBCHEM_MOLECULAR_FORMULA").unwrap();
        assert_eq!(prop2.expect_string(), "C9H8O4");
    }

    #[test]
    fn read_charges() {
        let path = Path::new("./src/tests-data/sdf/aspirin_charged.sdf");
        let mut trajectory = Trajectory::open(path).unwrap();

        let frame = trajectory.read().unwrap().unwrap();

        assert_approx_eq!(frame[0].charge, 0.0);
        assert_approx_eq!(frame[1].charge, 3.0);
        assert_approx_eq!(frame[2].charge, 2.0);
        assert_approx_eq!(frame[3].charge, 1.0);
        assert_approx_eq!(frame[4].charge, 0.0);
        assert_approx_eq!(frame[5].charge, -1.0);
        assert_approx_eq!(frame[6].charge, -2.0);
        assert_approx_eq!(frame[7].charge, -3.0);
        assert_approx_eq!(frame[8].charge, 0.0);
        assert_approx_eq!(frame[9].charge, 0.0);
        assert_approx_eq!(frame[10].charge, 0.0);
    }

    #[test]
    #[should_panic(expected = "atom line is too small for SDF")]
    fn too_small_atom_line() {
        let path = Path::new("./src/tests-data/sdf/bad-atom-line.sdf");
        let mut trajectory = Trajectory::open(path).unwrap();
        trajectory.read().unwrap().unwrap();
    }

    #[test]
    #[should_panic(expected = "could not parse bond count")]
    fn bad_counts_line() {
        let path = Path::new("./src/tests-data/sdf/count-line-not-numbers.sdf");
        let _trajectory = Trajectory::open(path).unwrap();
    }

    #[test]
    #[should_panic(expected = "counts line must have at least 10 characters")]
    fn count_line_too_short() {
        let path = Path::new("./src/tests-data/sdf/count-line-too-short.sdf");
        let _trajectory = Trajectory::open(path).unwrap();
    }

    #[test]
    fn write_files() {
        const EXPECTED: &str = r#"

created by molio
  4  3  0     0  0  0  0  0  0999 V2000
    1.0000    2.0000    3.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    2.0000    3.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    2.0000    3.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    2.0000    3.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  3  0  0  0  0
M  END
> <string-property>
prop1

$$$$
TEST

created by molio
 11  5  0     0  0  0  0  0  0999 V2000
    1.0000    2.0000    3.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    2.0000    3.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
    1.0000    2.0000    3.0000 C   0  2  0  0  0  0  0  0  0  0  0  0
    1.0000    2.0000    3.0000 F   0  1  0  0  0  0  0  0  0  0  0  0
    4.0000    5.0000    6.0000 E   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    5.0000    6.0000 D   0  5  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 G   0  6  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 H   0  7  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 J   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Xxx 0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  3  0  0  0  0
  9 10  8  0  0  0  0
 10 11  4  0  0  0  0
M  END
> <string-property>
prop1

$$$$
TEST

created by molio
  1  0  0     0  0  0  0  0  0999 V2000
    1.0000    2.0000    3.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <string-property>
prop1

$$$$
TEST

created by molio
  1  0  0     0  0  0  0  0  0999 V2000
    1.0000    2.0000    3.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <string-property>
prop1

$$$$
abc dfe ghi jkl mno pqr stu vwx yz 123 456 789 ABC DFE GHI JKL MNO PQR STU VWX Y

created by molio
  0  0  0     0  0  0  0  0  0999 V2000
M  END
$$$$
"#;
        let named_tmpfile = Builder::new()
            .prefix("test-sdf")
            .suffix(".sdf")
            .tempfile()
            .unwrap();

        // Write the expected content into the temp file
        let mut trajectory = Trajectory::create(named_tmpfile.path()).unwrap();

        let mut frame = Frame::new();
        frame.add_atom(
            Atom::with_symbol("A".to_string(), "O".to_string()),
            [1.0, 2.0, 3.0],
        );
        frame.add_atom(
            Atom::with_symbol("B".to_string(), "N".to_string()),
            [1.0, 2.0, 3.0],
        );
        frame.add_atom(Atom::new("C".to_string()), [1.0, 2.0, 3.0]);
        frame.add_atom(Atom::new("F".to_string()), [1.0, 2.0, 3.0]);
        frame.add_bond(0, 2, BondOrder::Single).unwrap();
        frame.add_bond(1, 2, BondOrder::Double).unwrap();
        frame.add_bond(2, 3, BondOrder::Triple).unwrap();
        frame.properties.insert(
            "string-property".to_string(),
            Property::String("prop1".to_string()),
        );
        trajectory.write(&frame).unwrap();

        frame.add_atom(Atom::new("E".to_string()), [4.0, 5.0, 6.0]);
        frame.add_atom(Atom::new("D".to_string()), [4.0, 5.0, 6.0]);
        frame.add_atom(Atom::new("G".to_string()), [0.0, 0.0, 0.0]);
        frame.add_atom(Atom::new("H".to_string()), [0.0, 0.0, 0.0]);
        frame.add_atom(Atom::new("I".to_string()), [0.0, 0.0, 0.0]);
        frame.add_atom(Atom::new("J".to_string()), [0.0, 0.0, 0.0]);
        frame.add_atom(Atom::new("".to_string()), [0.0, 0.0, 0.0]);
        frame[0].charge = 0.05;
        frame[1].charge = 1.0;
        frame[2].charge = 2.0;
        frame[3].charge = 3.0;
        frame[4].charge = 4.0;
        frame[5].charge = -1.0;
        frame[6].charge = -2.0;
        frame[7].charge = -3.0;

        frame.add_bond(9, 10, BondOrder::Aromatic).unwrap();
        frame.add_bond(8, 9, BondOrder::Unknown).unwrap();
        frame
            .properties
            .insert("name".to_string(), Property::String("TEST".to_string()));
        trajectory.write(&frame).unwrap();

        frame.clear_bonds();
        frame.resize(1).unwrap();

        trajectory.write(&frame).unwrap();

        trajectory.write(&frame).unwrap();

        // name is too long for SDF specification
        let mut frame = Frame::new();
        frame.properties.insert("name".to_string(), Property::String("abc dfe ghi jkl mno pqr stu vwx yz 123 456 789 ABC DFE GHI JKL MNO PQR STU VWX YZ 123 456 789".to_string()));
        trajectory.write(&frame).unwrap();

        let contents = std::fs::read_to_string(named_tmpfile.path()).unwrap();
        assert_eq!(contents, EXPECTED);

        frame
            .properties
            .insert("float property".to_string(), Property::Double(1.23));
        frame
            .properties
            .insert("bool property".to_string(), Property::Bool(false));
        frame.properties.insert(
            "vector property".to_string(),
            Property::Vector3D([1.0, 2.0, 3.0]),
        );
        trajectory.write(&frame).unwrap();

        // we have to check the properties separately as Properties is a HashMap thus
        // unordered
        let mut read_trajectory = Trajectory::open(named_tmpfile.path()).unwrap();
        let last_frame = read_trajectory.read_at(5).unwrap().unwrap();
        // properties are indiscriminately encoded as strings
        assert_eq!(
            last_frame
                .properties
                .get("float property")
                .unwrap()
                .expect_string(),
            "1.23"
        );
        assert_eq!(
            last_frame
                .properties
                .get("bool property")
                .unwrap()
                .expect_string(),
            "false"
        );
        let v = last_frame
            .properties
            .get("vector property")
            .unwrap()
            .expect_string();
        assert_eq!(v, "1.0 2.0 3.0");

        trajectory.finish().unwrap();
    }
}
