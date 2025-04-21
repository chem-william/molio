use crate::atom::Atom;
use crate::error::CError;
use crate::format::FileFormat;
use crate::frame::Frame;
use crate::property::Properties;
use crate::property::Property;
use crate::unit_cell::UnitCell;
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Seek};
use std::path::Path;

#[derive(Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct FullResidueId {
    /// Chain identifier
    chain: char,
    /// Residue id
    resid: i64,
    /// Residue name
    resname: String,
    /// Insertion code of the residue
    insertion_code: char,
}

#[derive(Default)]
pub struct Residue {
    name: String,
    id: Option<i64>,
    atoms: BTreeSet<usize>,
    properties: Properties,
}

impl Residue {
    fn add_atom(&mut self, index: usize) {
        self.atoms.insert(index);
    }
}

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

pub fn get_record(line: &str) -> Record {
    // let rec: String = line.chars().take(6).collect();

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

    let f = line.chars().next().unwrap();
    if let Some(f) = line.chars().next() {
        if f == '-' || f == ' ' {
            // Negative or space-prefixed numbers are not encoded → return 0
            return Ok(0);
        } else if f.is_ascii_digit() {
            // Try to parse as a number
            return line.parse::<i64>().map_err(|_e| {
                CError::GenericError(format!("expected negative number. got '{f}'"))
            });
        }
    }

    if line.is_empty() {
        return Ok(0);
    }

    if DIGITS_UPPER.contains(f) {
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

    Ok(1)
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
    pub secinfo: BTreeMap<FullResidueId, (FullResidueId, String)>,
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
            let initial_offset = decode_hybrid36(5, &line[6..11]);
            println!(
                "{} is not a valid atom id, assuming '1'",
                initial_offset.as_ref().unwrap()
            );
            self.atom_offsets.borrow_mut().push(0);
            if initial_offset.is_err() {
                println!(
                    "{} is not a valid atom id, assuming '1'",
                    initial_offset.as_ref().unwrap()
                );
                self.atom_offsets.borrow_mut().push(0);
            }

            if *initial_offset.as_ref().unwrap() <= 0 {
                println!(
                    "warning: '{}' is too small, assuming id is '1'",
                    initial_offset.as_ref().unwrap()
                );
                self.atom_offsets.borrow_mut().push(0);
            } else {
                self.atom_offsets.borrow_mut().push(
                    usize::try_from(*initial_offset.as_ref().unwrap())
                        .expect("decode_hybrid36 returned a negative number"),
                );
            };
        }

        let mut atom = Atom::default();
        let name = &line[12..16];
        if line.len() >= 78 {
            let atom_type = &line[76..78];
            atom.name = name.trim().to_string();
            atom.symbol = atom_type.to_string();
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
            if let Some(secinfo_for_residue) = self.secinfo.get(&full_residue_id) {
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
                println!("ignoring custom spce group ({space_group}), using P1 instead");
            }
        }

        Ok(())
    }
    pub fn new() -> Self {
        PDBFormat {
            residues: RefCell::new(BTreeMap::new()),
            atom_offsets: RefCell::new(Vec::new()),
            current_secinfo: RefCell::new(None),
            secinfo: BTreeMap::new(),
        }
    }
    const END_RECORD: &str = "END";
    const ENDMDL_RECORD: &str = "ENDMDL";
}

impl FileFormat for PDBFormat {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        self.residues.borrow_mut().clear();
        let mut line = String::new();
        let mut frame = Frame::new();

        while reader.read_line(&mut line)? > 0 {
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
                // Record::HETATM => {}
                // Record::CONECT => {}
                // Record::MODEL => {}
                // Record::ENDMDL => {}
                // Record::HELIX => {}
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

    use crate::{formats::pdb::decode_hybrid36, trajectory::Trajectory, unit_cell::CellShape};

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/pdb/water.pdb");
        let trajectory = Trajectory::new(path).unwrap();
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

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 297);

        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 0.417, 1e-5);
        assert_approx_eq!(positions[0][1], 8.303, 1e-5);
        assert_approx_eq!(positions[0][2], 11.737, 1e-5);

        let cell = frame.unit_cell;
        assert_eq!(cell.shape, CellShape::Orthorhombic);
        assert_approx_eq!(cell.lengths()[0], 15.0);
        assert_approx_eq!(cell.lengths()[1], 15.0);
        assert_approx_eq!(cell.lengths()[2], 15.0);
    }

    #[test]
    fn hybrid_decode() {
        assert_eq!(decode_hybrid36(4, "    ").unwrap(), 0);
        assert_eq!(decode_hybrid36(4, "  -0").unwrap(), 0);
    }
}
