use crate::atom::Atom;
use crate::bond::BondOrder;
use crate::error::CError;
use crate::format::FileFormat;
use crate::frame::Frame;
use crate::property::Property;
use crate::residue::{FullResidueId, Residue};
use crate::unit_cell::UnitCell;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek};
use std::path::Path;

#[derive(Debug, PartialEq, Eq)]
pub enum Record {
    // Records containing summary data
    HEADER,
    TITLE,
    // Records containing useful data
    CRYST1,
    ATOM,
    HETATM,
    CONECT,
    // Beginning of model
    MODEL,
    // End of model
    ENDMDL,
    // End of chain. May increase atom count
    TER,
    // End of file
    END,
    // Secondary structure
    HELIX,
    SHEET,
    TURN,
    // Ignored records
    IGNORED_,
    // Unknown record type
    UNKNOWN_,
}

/// See http://www.wwpdb.org/documentation/file-format-content/format23/sect5.html
/// for definitions of helix types
pub enum HelixType {
    RightHandedAlphaHelix,
    RightHandedOmegaHelix,
    RightHandedPiHelix,
    RightHandedGammaHelix,
    RightHanded3_10Helix,
    LeftHandedAlphaHelix,
    LeftHandedOmegaHelix,
    LeftHandedGammaHelix,
    Two7RibbonHelix,
    Polyproline,
}

impl HelixType {
    pub fn nth(n: usize) -> Option<&'static str> {
        match n {
            0 => Some("right handed alpha helix"),
            1 => Some("right handed omega helix"),
            2 => Some("right handed pi helix"),
            3 => Some("right handed gamma helix"),
            4 => Some("right handed 3-10 helix"),
            5 => Some("left handed alpha helix"),
            6 => Some("left handed omega helix"),
            7 => Some("left handed gamma helix"),
            8 => Some("two 7 ribbon helix"),
            9 => Some("polyproline"),
            _ => None,
        }
    }
}

pub fn get_record(line: &str) -> Record {
    let rec = if line.len() >= 6 { &line[..6] } else { line };

    match rec {
        "ENDMDL" => Record::ENDMDL,
        _ if rec.starts_with("END") => Record::END,
        "CRYST1" => Record::CRYST1,
        "ATOM  " => Record::ATOM,
        "HETATM" => Record::HETATM,
        "CONECT" => Record::CONECT,
        _ if rec.starts_with("MODEL") => Record::MODEL,
        _ if rec.starts_with("TER") => Record::TER,
        "HELIX " => Record::HELIX,
        "SHEET " => Record::SHEET,
        "TURN  " => Record::TURN,
        "HEADER" => Record::HEADER,
        "TITLE " => Record::TITLE,
        "REMARK" | "MASTER" | "AUTHOR" | "CAVEAT" | "COMPND" | "EXPDTA" | "KEYWDS" | "OBSLTE"
        | "SOURCE" | "SPLIT " | "SPRSDE" | "JRNL  " | "SEQRES" | "HET   " | "REVDAT" | "SCALE1"
        | "SCALE2" | "SCALE3" | "ORIGX1" | "ORIGX2" | "ORIGX3" | "ANISOU" | "SITE  " | "FORMUL"
        | "DBREF " | "HETNAM" | "HETSYN" | "SSBOND" | "LINK  " | "SEQADV" | "MODRES" | "CISPEP" => {
            Record::IGNORED_
        }
        _ if line.trim().is_empty() => Record::IGNORED_,
        _ => Record::UNKNOWN_,
    }
}

const DIGITS_UPPER: &str = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const DIGITS_LOWER: &str = "0123456789abcdefghijklmnopqrstuvwxyz";

fn decode_pure(s: &str) -> i64 {
    let mut result: i64 = 0;
    let base = 36; // Supports 0-9, A-Z, a-z, but lowercase maps like uppercase

    for c in s.chars() {
        result *= base;
        result += digit_to_value(c) as i64;
    }

    result
}

fn encode_pure(digits: &str, value: i64) -> String {
    if value == 0 {
        return digits.chars().nth(1).unwrap().to_string();
    };

    let n = i64::try_from(digits.len()).unwrap();
    let mut result = String::new();
    let mut val = value;

    while val != 0 {
        let rest = val / n;
        let rem = usize::try_from(val - rest * n).unwrap();
        let ch = digits.chars().nth(rem).expect("digit index out of bounds");
        result.push(ch);
        val = rest;
    }

    result.chars().rev().collect()
}

fn digit_to_value(c: char) -> i32 {
    match c {
        '0'..='9' => (c as u8 - b'0') as i32,
        'A'..='Z' => (c as u8 - b'A') as i32 + 10,
        'a'..='z' => (c as u8 - b'a') as i32 + 10,
        _ => panic!("Invalid character: {}", c),
    }
}

fn pow_int(base: usize, power: usize) -> i64 {
    (base.pow(power as u32)) as i64
}

// TODO: maybe this should be changed to u64? can it even return negative numbers?
pub(crate) fn decode_hybrid36(width: usize, line: &str) -> Result<i64, CError> {
    if line.len() > width {
        return Err(CError::GenericError(format!(
            "length of '{line}' is greater than the width '{width}'. this is a bug"
        )));
    }

    let f = line.chars().next();
    if let Some(f) = f {
        if f == ' ' {
            // Space-prefixed numbers are not encoded => return 0
            return Ok(0);
        } else if f == '-' || f.is_ascii_digit() {
            // Try to parse as a number
            return line.parse::<i64>().map_err(|_e| {
                CError::GenericError(format!("expected negative number. got '{f}'"))
            });
        }
    }

    if line.is_empty() {
        return Ok(0);
    }

    if DIGITS_UPPER.contains(f.unwrap()) {
        let is_valid = line
            .chars()
            .all(|c| c.is_ascii_digit() || c.is_ascii_uppercase());

        if !is_valid {
            return Err(CError::GenericError(format!(
                "the value '{}' is not a valid hybrid 36 number",
                line
            )));
        }

        return Ok(decode_pure(line) - 10 * pow_int(36, width - 1) + pow_int(10, width));
    }

    if DIGITS_LOWER.contains(f.unwrap()) {
        let is_valid = line
            .chars()
            .all(|c| c.is_ascii_digit() || c.is_ascii_lowercase());

        if !is_valid {
            return Err(CError::GenericError(format!(
                "the value '{}' is not a valid hybrid 36 number",
                line
            )));
        }

        return Ok(decode_pure(line) + 16 * pow_int(36, width - 1) + pow_int(10, width));
    }

    Err(CError::GenericError(format!(
        "the value '{line}' is not a valid hybrid 36 number"
    )))
}

pub(crate) fn encode_hybrid36(width: usize, value: i64) -> String {
    // the number is too negative to be encoded
    if value < (1 - pow_int(10, width - 1)) {
        return "*".repeat(width);
    };

    // no need to encode
    if value < pow_int(10, width) {
        return value.to_string();
    };

    // use upper case set
    let mut val = value;
    val -= pow_int(10, width);
    if val < 26 * pow_int(36, width - 1) {
        val += 10 * pow_int(36, width - 1);
        return encode_pure(DIGITS_UPPER, val);
    };

    // use lower case set
    val -= 26 * pow_int(36, width - 1);
    if val < 26 * pow_int(36, width - 1) {
        val += 10 * pow_int(36, width - 1);
        return encode_pure(DIGITS_LOWER, val);
    }

    // too large to be encoded
    "*".repeat(width)
}

pub struct PDBFormat {
    /// Residue information in the current step
    pub residues: RefCell<BTreeMap<FullResidueId, Residue>>,

    /// List of all atom offsets. This maybe pushed in read_ATOM or if a TER
    /// record is found. It is reset every time a frame is read.
    pub atom_offsets: RefCell<Vec<usize>>,

    /// This will be None when no secondary structure information should be
    /// read. Else it is set to the final residue of a secondary structure and
    /// the text description which should be set.
    pub current_secinfo: RefCell<Option<(FullResidueId, String)>>,

    /// Store secondary structure information. Keys are the
    /// starting residue of the secondary structure, and values are pairs
    /// containing the ending residue and a string which is a written
    /// description of the secondary structure
    pub secinfo: RefCell<BTreeMap<FullResidueId, (FullResidueId, String)>>,

    /// Number of models read/written to the file
    models: RefCell<usize>,
}

impl PDBFormat {
    fn parse_atom(&self, frame: &mut Frame, line: &str, is_hetatm: bool) -> Result<(), CError> {
        debug_assert!(line[..6] == *"ATOM  " || line[..6] == *"HETATM");

        if line.len() < 54 {
            return Err(CError::InvalidRecord {
                expected_record_type: "ATOM".to_string(),
                actual_record_type: line.to_string(),
                reason: "line too short".to_string(),
            });
        }

        if self.atom_offsets.borrow().is_empty() {
            let initial_offset = decode_hybrid36(5, line[6..11].trim()).unwrap();

            if initial_offset <= 0 {
                println!(
                    "warning: '{}' is too small, assuming id is '1'",
                    initial_offset
                );
                self.atom_offsets.borrow_mut().push(0);
            } else {
                self.atom_offsets.borrow_mut().push(
                    usize::try_from(initial_offset)
                        .expect("decode_hybrid36 returned a negative number")
                        - 1,
                );
            };
        }

        let mut atom = Atom::default();
        let name = &line[12..16];
        if line.len() >= 78 {
            let atom_type = &line[76..78];
            atom.name = name.trim().to_string();
            atom.symbol = atom_type.trim().to_string();
        } else {
            // Read just the atom name and hope for the best
            atom.name = name.to_string();
        }

        let altloc = &line[16..17];
        if altloc != " " {
            atom.properties
                .insert("altloc".to_string(), Property::String(altloc.to_string()));
        }

        let x = Property::parse_value(&line[30..38], crate::property::PropertyKind::Double)?;
        let y = Property::parse_value(&line[38..46], crate::property::PropertyKind::Double)?;
        let z = Property::parse_value(&line[46..54], crate::property::PropertyKind::Double)?;
        atom.x = x.expect_double();
        atom.y = y.expect_double();
        atom.z = z.expect_double();
        frame.add_atom(atom);

        let atom_id = frame.size() - 1;
        let insertion_code = line.chars().nth(26).unwrap();
        let resid = match decode_hybrid36(4, &line[22..26]) {
            Ok(resid) => resid,
            // No residue information so return early
            Err(_) => return Ok(()),
        };

        let chain = &line.chars().nth(21).unwrap();
        let resname = line[17..20].trim().to_string();
        let full_residue_id = FullResidueId {
            chain: *chain,
            resid,
            resname: resname.to_string(),
            ..Default::default()
        };

        if self.residues.borrow().len() == 0 {
            let mut residue = Residue {
                name: resname,
                id: Some(resid),
                ..Default::default()
            };
            residue.add_atom(atom_id);

            if insertion_code != ' ' {
                residue.properties.insert(
                    "insertion_code".to_string(),
                    Property::String(insertion_code.to_string()),
                );
            }

            // Set whether or not the residue is standardized
            residue
                .properties
                .insert("is_standard_pdb".to_string(), Property::Bool(!is_hetatm));

            // This is saved as a string (instead of a number) on purpose
            // to match MMTF format
            residue.properties.insert(
                "chainid".to_string(),
                Property::String(line.chars().nth(21).unwrap().to_string()),
            );

            // PDB format makes no distinction between chainid and chainname
            residue.properties.insert(
                "chainname".to_string(),
                Property::String(line.chars().nth(21).unwrap().to_string()),
            );

            // segment name is not part of the standard, but something added by
            // CHARM/NAMD in the un-used character range 73-76
            if line.len() > 72 {
                let segname = line[72..76].trim();
                if !segname.is_empty() {
                    residue
                        .properties
                        .insert("segname".to_string(), Property::String(segname.to_string()));
                }
            }

            // Are we within a secondary information sequence?
            if let Some(secinfo) = self.current_secinfo.borrow().as_ref() {
                residue.properties.insert(
                    "secondary_structure".to_string(),
                    Property::String(secinfo.1.clone()),
                );

                if secinfo.0 == full_residue_id {
                    self.current_secinfo.replace(None);
                }
            }

            // Are we the start of a secondary information sequence?
            if let Some(secinfo_for_residue) = self.secinfo.borrow().get(&full_residue_id) {
                self.current_secinfo.replace(None);
                residue.properties.insert(
                    "secondary_structure".to_string(),
                    Property::String(secinfo_for_residue.1.clone()),
                );
            }

            self.residues.borrow_mut().insert(full_residue_id, residue);
        } else {
            // Just add this atom to the residue
            self.residues
                .borrow_mut()
                .get_mut(&full_residue_id)
                .unwrap()
                .add_atom(atom_id);
        }

        Ok(())
    }

    fn parse_cryst1(frame: &mut Frame, line: &str) -> Result<(), CError> {
        assert_eq!(&line[..6], "CRYST1");
        if line.len() < 54 {
            return Err(CError::InvalidRecord {
                expected_record_type: "CRYST1".to_string(),
                actual_record_type: line.to_string(),
                reason: "line too short".to_string(),
            });
        }

        let a = line[6..15]
            .trim()
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let b = line[15..24]
            .trim()
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let c = line[24..33]
            .trim()
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let alpha = line[33..40]
            .trim()
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let beta = line[40..47]
            .trim()
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let gamma = line[47..54]
            .trim()
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;

        let unit_cell = UnitCell::new_from_lengths_angles([a, b, c], &mut [alpha, beta, gamma])?;
        frame.unit_cell = unit_cell;

        if line.len() >= 55 {
            let space_group = &line[55..65].trim();
            // TODO: handle this as a warning (somehow)?
            if space_group != &"P 1" && space_group != &"P1" {
                println!("ignoring custom space group ({space_group}), using P1 instead");
            }
        }

        Ok(())
    }
    fn _read_index(&self, line: &str, initial: usize) -> Result<usize, CError> {
        let trimmed_line = &line[initial..initial + 5].trim();
        let mut pdb_atom_id = decode_hybrid36(5, trimmed_line)?;

        // Find the lower bound index
        let lower = self
            .atom_offsets
            .borrow()
            .binary_search(&(pdb_atom_id as usize))
            .unwrap_or_else(|insert_pos| insert_pos)
            - 1;

        pdb_atom_id -= lower as i64 - self.atom_offsets.borrow().first().copied().unwrap() as i64;

        // TODO: is this correct?
        pdb_atom_id -= 1;
        Ok(usize::try_from(pdb_atom_id).unwrap())
    }

    fn _add_bond(&self, frame: &mut Frame, line: &str, i: usize, j: usize) {
        if i >= frame.size() || j >= frame.size() {
            println!(
                "warning: PDB reader: ignoring CONECT ('{}') with atomic indexes bigger than frame size ({})",
                line.trim(),
                frame.size()
            );
            return;
        }
        frame
            .add_bond(i, j, BondOrder::Unknown)
            .expect("could not add bond on frame when parsing CONECT");
    }

    fn parse_conect(&self, frame: &mut Frame, line: &str) -> Result<(), CError> {
        debug_assert_eq!(&line[..6], "CONECT");

        let line_length = line.trim().len();

        let index_i = self
            ._read_index(line, 6)
            .expect("could not read index when parsing CONECT");

        if line_length > 11 {
            let index_j = self._read_index(line, 11).unwrap();
            self._add_bond(frame, line, index_i, index_j);
        };

        if line_length > 16 {
            let index_j = self._read_index(line, 16).unwrap();
            self._add_bond(frame, line, index_i, index_j);
        };

        if line_length > 21 {
            let index_j = self._read_index(line, 21).unwrap();
            self._add_bond(frame, line, index_i, index_j);
        };

        if line_length > 26 {
            let index_j = self._read_index(line, 26).unwrap();
            self._add_bond(frame, line, index_i, index_j);
        };

        Ok(())
    }

    fn parse_helix(&self, line: &str) -> Result<(), CError> {
        if line.len() < 33 + 5 {
            println!("warning: HELIX record too short: {line}");
        };

        let chain_start = line.chars().nth(19).expect("chain start");
        let chain_end = line.chars().nth(31).expect("chain end");
        let inscode_start = line.chars().nth(25).expect("inscode start");
        let inscode_end = line.chars().nth(37).expect("inscode end");
        let resname_start = line[15..18].trim();
        let resname_end = line[27..30].trim();

        let resid_start = decode_hybrid36(4, &line[21..25])?;
        let resid_end = decode_hybrid36(4, &line[33..37])?;

        if chain_start != chain_end {
            println!("warning: HELIX chain {chain_start} and {chain_end} are not the same");
        }

        let start = FullResidueId {
            chain: chain_start,
            resid: resid_start,
            resname: resname_start.to_string(),
            insertion_code: inscode_start,
        };
        let end = FullResidueId {
            chain: chain_end,
            resid: resid_end,
            resname: resname_end.to_string(),
            insertion_code: inscode_end,
        };

        let helix_type = &line[38..40]
            .parse::<usize>()
            .inspect_err(|e| eprintln!("failed to parse helix type: {e}"))
            .unwrap();
        if *helix_type <= 10 {
            self.secinfo.borrow_mut().insert(
                start,
                (
                    end,
                    HelixType::nth(2).expect("unknown helix type").to_string(),
                ),
            );
        }

        Ok(())
    }

    pub fn new() -> Self {
        PDBFormat {
            residues: RefCell::new(BTreeMap::new()),
            atom_offsets: RefCell::new(Vec::new()),
            current_secinfo: RefCell::new(None),
            secinfo: RefCell::new(BTreeMap::new()),
            models: RefCell::new(0),
        }
    }

    const END_RECORD: &str = "END";
    const ENDMDL_RECORD: &str = "ENDMDL";
}

impl FileFormat for PDBFormat {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        self.residues.borrow_mut().clear();
        let mut frame = Frame::new();
        let mut line = String::new();

        let mut got_end = false;
        while !got_end && reader.read_line(&mut line)? > 0 {
            let record = get_record(&line);
            let name = String::new();

            match record {
                Record::HEADER => {
                    if line.len() >= 50 {
                        frame.properties.insert(
                            "classification".to_string(),
                            Property::String(line[10..50].trim().to_string()),
                        );
                    }
                    if line.len() >= 59 {
                        frame.properties.insert(
                            "deposition_date".to_string(),
                            Property::String(line[50..59].trim().to_string()),
                        );
                    }
                    if line.len() >= 66 {
                        frame.properties.insert(
                            "pdb_idcode".to_string(),
                            Property::String(line[62..66].trim().to_string()),
                        );
                    }
                }
                Record::TITLE => {
                    if line.len() < 11 {
                        continue;
                    }
                    let title = line[10..80].trim();
                    let current = frame
                        .properties
                        .get("name")
                        .and_then(|p| p.as_string())
                        .unwrap_or_default();
                    let new_title = if current.is_empty() {
                        title.to_string()
                    } else {
                        format!("{} {}", current, title)
                    };
                    frame
                        .properties
                        .insert("name".to_string(), Property::String(new_title));
                }
                Record::CRYST1 => PDBFormat::parse_cryst1(&mut frame, &line).unwrap(),
                Record::ATOM => self.parse_atom(&mut frame, &line, false).unwrap(),
                Record::HETATM => self.parse_atom(&mut frame, &line, true).unwrap(),
                Record::CONECT => self.parse_conect(&mut frame, &line).unwrap(),
                Record::MODEL => *self.models.borrow_mut() += 1,
                Record::ENDMDL => {
                    // look one line ahead to see if the next Record is an `END`
                    let mut line = String::new();
                    let bytes = reader.read_line(&mut line)?;
                    if bytes > 0 {
                        reader.seek_relative(
                            -(i64::try_from(bytes).expect("failed to convert bytes offset")),
                        )?;
                        if get_record(&line) == Record::END {
                            // If that is the case then wait for the next Record
                            continue;
                        }
                    }
                    got_end = true;
                }
                Record::HELIX => self.parse_helix(&line).unwrap(),
                // Record::SHEET => {}
                // Record::TURN => {}
                // Record::TER => {}
                Record::END => break,
                // Record::IGNORED_ => {}
                // Record::UNKNOWN_ => {}
                _ => {}
            }
            line.clear();
        }

        Ok(frame)
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        let mut line = String::new();

        while reader.read_line(&mut line)? > 0 {
            if line.starts_with(Self::ENDMDL_RECORD) {
                let mut next_line = String::new();
                let bytes = reader.read_line(&mut next_line)?;

                reader.seek_relative(
                    -(i64::try_from(bytes).expect("failed to convert bytes offset")),
                )?;
                if next_line.trim().starts_with(Self::END_RECORD) {
                    next_line.clear();
                    continue;
                }
            }

            if line.starts_with(Self::END_RECORD) {
                return Ok(Some(reader.stream_position()?));
            }

            line.clear();
        }

        Ok(None)
    }

    fn write(&self, path: &Path, frame: &Frame) -> Result<(), CError> {
        println!(
            "Writing {:?} as PDB format with {} atoms",
            path,
            frame.size()
        );
        Ok(())
    }

    fn read(&self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        if reader.fill_buf().map(|b| !b.is_empty()).unwrap() {
            Ok(Some(self.read_next(reader).unwrap()))
        } else {
            Ok(None)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use assert_approx_eq::assert_approx_eq;

    use crate::{
        angle::Angle,
        bond::Bond,
        dihedral::Dihedral,
        formats::pdb::{decode_hybrid36, encode_hybrid36},
        topology::Topology,
        trajectory::Trajectory,
        unit_cell::CellShape,
    };

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/pdb/water.pdb");
        let mut trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 100);

        let path = Path::new("./src/tests-data/pdb/2hkb.pdb");
        let trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 11);
    }

    #[test]
    fn sanity_check() {
        let path = Path::new("./src/tests-data/pdb/water.pdb");
        let mut trajectory = Trajectory::new(path).unwrap();
        assert_eq!(trajectory.size, 100);

        let mut frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 297);

        let mut positions = frame.positions();
        assert_approx_eq!(positions[0][0], 0.417, 1e-3);
        assert_approx_eq!(positions[0][1], 8.303, 1e-3);
        assert_approx_eq!(positions[0][2], 11.737, 1e-3);

        assert_approx_eq!(positions[296][0], 6.664, 1e-3);
        assert_approx_eq!(positions[296][1], 11.6148, 1e-3);
        assert_approx_eq!(positions[296][2], 12.961, 1e-3);

        let cell = frame.unit_cell;
        assert_eq!(cell.shape, CellShape::Orthorhombic);
        assert_approx_eq!(cell.lengths()[0], 15.0);
        assert_approx_eq!(cell.lengths()[1], 15.0);
        assert_approx_eq!(cell.lengths()[2], 15.0);

        trajectory.read().unwrap().unwrap();
        frame = trajectory.read().unwrap().unwrap();

        assert_eq!(frame.size(), 297);
        positions = frame.positions();

        assert_approx_eq!(positions[0][0], 0.299, 1e-4);
        assert_approx_eq!(positions[0][1], 8.310, 1e-4);
        assert_approx_eq!(positions[0][2], 11.721, 1e-4);

        assert_approx_eq!(positions[296][0], 6.798, 1e-4);
        assert_approx_eq!(positions[296][1], 11.509, 1e-4);
        assert_approx_eq!(positions[296][2], 12.704, 1e-4);
    }

    #[test]
    fn read_bonds() {
        let path = Path::new("./src/tests-data/pdb/MOF-5.pdb");
        let mut trajectory = Trajectory::new(path).unwrap();
        let mut frame = trajectory.read().unwrap().unwrap();

        let topology: &mut Topology = frame.topology_as_mut();
        assert_eq!(topology.size(), 65);

        assert_eq!(topology[0].symbol, "Zn");
        assert_eq!(topology[1].symbol, "O");

        assert_eq!(topology[0].name, "ZN");
        assert_eq!(topology[1].name, "O");

        assert_eq!(topology.bonds().len(), 68);

        assert!(topology.bonds().contains(&Bond::new(9, 38)));
        assert!(topology.bonds().contains(&Bond::new(58, 62)));
        assert!(topology.bonds().contains(&Bond::new(37, 24)));
        assert!(topology.bonds().contains(&Bond::new(27, 31)));

        println!("angles: {:?}", topology.angles());
        assert!(topology.angles().contains(&Angle::new(20, 21, 23)));
        assert!(topology.angles().contains(&Angle::new(9, 38, 44)));

        assert!(
            topology
                .dihedrals()
                .contains(&Dihedral::new(64, 62, 58, 53))
        );
        assert!(
            topology
                .dihedrals()
                .contains(&Dihedral::new(22, 21, 23, 33))
        );
    }

    macro_rules! recycle_check {
        ($width:expr, $value:expr, $hybrid:expr) => {
            assert_eq!(encode_hybrid36($width, $value), $hybrid);
            assert_eq!(
                decode_hybrid36($width, $hybrid).expect(
                    format!(
                        "decode failed with width '{0}' and hybrid '{1}'",
                        $width, $hybrid
                    )
                    .as_str()
                ),
                $value
            );
        };
    }

    #[test]
    fn hybrid_encode_decode() {
        assert_eq!(decode_hybrid36(4, "    ").unwrap(), 0);
        assert_eq!(decode_hybrid36(4, "  -0").unwrap(), 0);

        recycle_check!(4, -999, "-999");
        recycle_check!(4, -78, "-78");
        recycle_check!(4, -6, "-6");
        recycle_check!(4, 0, "0");
        recycle_check!(4, 9999, "9999");
        recycle_check!(4, 10000, "A000");
        recycle_check!(4, 10001, "A001");
        recycle_check!(4, 10002, "A002");
        recycle_check!(4, 10003, "A003");
        recycle_check!(4, 10004, "A004");
        recycle_check!(4, 10005, "A005");
        recycle_check!(4, 10006, "A006");
        recycle_check!(4, 10007, "A007");
        recycle_check!(4, 10008, "A008");
        recycle_check!(4, 10009, "A009");
        recycle_check!(4, 10010, "A00A");
        recycle_check!(4, 10011, "A00B");
        recycle_check!(4, 10012, "A00C");
        recycle_check!(4, 10013, "A00D");
        recycle_check!(4, 10014, "A00E");
        recycle_check!(4, 10015, "A00F");
        recycle_check!(4, 10016, "A00G");
        recycle_check!(4, 10017, "A00H");
        recycle_check!(4, 10018, "A00I");
        recycle_check!(4, 10019, "A00J");
        recycle_check!(4, 10020, "A00K");
        recycle_check!(4, 10021, "A00L");
        recycle_check!(4, 10022, "A00M");
        recycle_check!(4, 10023, "A00N");
        recycle_check!(4, 10024, "A00O");
        recycle_check!(4, 10025, "A00P");
        recycle_check!(4, 10026, "A00Q");
        recycle_check!(4, 10027, "A00R");
        recycle_check!(4, 10028, "A00S");
        recycle_check!(4, 10029, "A00T");
        recycle_check!(4, 10030, "A00U");
        recycle_check!(4, 10031, "A00V");
        recycle_check!(4, 10032, "A00W");
        recycle_check!(4, 10033, "A00X");
        recycle_check!(4, 10034, "A00Y");
        recycle_check!(4, 10035, "A00Z");
        recycle_check!(4, 10036, "A010");
        recycle_check!(4, 10046, "A01A");
        recycle_check!(4, 10071, "A01Z");
        recycle_check!(4, 10072, "A020");
        recycle_check!(4, 10000 + 36 * 36 - 1, "A0ZZ");
        recycle_check!(4, 10000 + 36 * 36, "A100");
        recycle_check!(4, 10000 + 36 * 36 * 36 - 1, "AZZZ");
        recycle_check!(4, 10000 + 36 * 36 * 36, "B000");
        recycle_check!(4, 10000 + 26 * 36 * 36 * 36 - 1, "ZZZZ");
        recycle_check!(4, 10000 + 26 * 36 * 36 * 36, "a000");
        recycle_check!(4, 10000 + 26 * 36 * 36 * 36 + 35, "a00z");
        recycle_check!(4, 10000 + 26 * 36 * 36 * 36 + 36, "a010");
        recycle_check!(4, 10000 + 26 * 36 * 36 * 36 + 36 * 36 - 1, "a0zz");
        recycle_check!(4, 10000 + 26 * 36 * 36 * 36 + 36 * 36, "a100");
        recycle_check!(4, 10000 + 26 * 36 * 36 * 36 + 36 * 36 * 36 - 1, "azzz");
        recycle_check!(4, 10000 + 26 * 36 * 36 * 36 + 36 * 36 * 36, "b000");
        recycle_check!(4, 10000 + 2 * 26 * 36 * 36 * 36 - 1, "zzzz");

        assert_eq!(decode_hybrid36(5, "    ").unwrap(), 0);
        assert_eq!(decode_hybrid36(5, "  -0").unwrap(), 0);
        recycle_check!(5, -9999, "-9999");
        recycle_check!(5, -123, "-123");
        recycle_check!(5, -45, "-45");
        recycle_check!(5, -6, "-6");
        recycle_check!(5, 0, "0");
        recycle_check!(5, 12, "12");
        recycle_check!(5, 345, "345");
        recycle_check!(5, 6789, "6789");
        recycle_check!(5, 99999, "99999");
        recycle_check!(5, 100000, "A0000");
        recycle_check!(5, 100010, "A000A");
        recycle_check!(5, 100035, "A000Z");
        recycle_check!(5, 100036, "A0010");
        recycle_check!(5, 100046, "A001A");
        recycle_check!(5, 100071, "A001Z");
        recycle_check!(5, 100072, "A0020");
        recycle_check!(5, 100000 + 36 * 36 - 1, "A00ZZ");
        recycle_check!(5, 100000 + 36 * 36, "A0100");
        recycle_check!(5, 100000 + 36 * 36 * 36 - 1, "A0ZZZ");
        recycle_check!(5, 100000 + 36 * 36 * 36, "A1000");
        recycle_check!(5, 100000 + 36 * 36 * 36 * 36 - 1, "AZZZZ");
        recycle_check!(5, 100000 + 36 * 36 * 36 * 36, "B0000");
        recycle_check!(5, 100000 + 2 * 36 * 36 * 36 * 36, "C0000");
        recycle_check!(5, 100000 + 26 * 36 * 36 * 36 * 36 - 1, "ZZZZZ");
        recycle_check!(5, 100000 + 26 * 36 * 36 * 36 * 36, "a0000");
        recycle_check!(5, 100000 + 26 * 36 * 36 * 36 * 36 + 36 - 1, "a000z");
        recycle_check!(5, 100000 + 26 * 36 * 36 * 36 * 36 + 36, "a0010");
        recycle_check!(5, 100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 - 1, "a00zz");
        recycle_check!(5, 100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36, "a0100");
        recycle_check!(
            5,
            100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36 - 1,
            "a0zzz"
        );
        recycle_check!(5, 100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36, "a1000");
        recycle_check!(
            5,
            100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36 * 36 - 1,
            "azzzz"
        );
        recycle_check!(
            5,
            100000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36 * 36,
            "b0000"
        );
        recycle_check!(5, 100000 + 2 * 26 * 36 * 36 * 36 * 36 - 1, "zzzzz");

        assert_eq!(encode_hybrid36(4, -99999), "****");
        assert_eq!(encode_hybrid36(4, 9999999), "****");
    }

    #[test]
    #[should_panic(expected = "the value '*0000' is not a valid hybrid 36 number")]
    fn decode_bad1() {
        decode_hybrid36(5, "*0000").unwrap();
    }

    #[test]
    #[should_panic(expected = "the value 'A*000' is not a valid hybrid 36 number")]
    fn decode_bad2() {
        decode_hybrid36(5, "A*000").unwrap();
    }

    #[test]
    #[should_panic(expected = "the value 'a*000' is not a valid hybrid 36 number")]
    fn decode_bad3() {
        decode_hybrid36(5, "a*000").unwrap();
    }

    #[test]
    #[should_panic(expected = "length of '12345' is greater than the width '2'. this is a bug")]
    fn decode_bad4() {
        decode_hybrid36(2, "12345").unwrap();
    }
}
