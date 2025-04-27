use crate::atom::Atom;
use crate::bond::BondOrder;
use crate::error::CError;
use crate::format::FileFormat;
use crate::frame::Frame;
use crate::property::Property;
use crate::property::PropertyKind;
use crate::residue::{FullResidueId, Residue};
use crate::unit_cell::UnitCell;
use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Seek, Write};

use super::pdb_connectivity::{self, find};

/// Maximum value for a width 4 number
const MAX_HYBRID36_W4_NUMBER: i64 = 2_436_111;

/// Maximum value for a width 5 number
const MAX_HYBRID36_W5_NUMBER: i64 = 87_440_031;

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

/// See <http://www.wwpdb.org/documentation/file-format-content/format23/sect5.html>
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
            0 => Some("right-handed alpha helix"),
            1 => Some("right-handed omega helix"),
            2 => Some("right-handed pi helix"),
            3 => Some("right-handed gamma helix"),
            4 => Some("right-handed 3-10 helix"),
            5 => Some("left-handed alpha helix"),
            6 => Some("left-handed omega helix"),
            7 => Some("left-handed gamma helix"),
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
        result += i64::from(digit_to_value(c));
    }

    result
}

fn encode_pure(digits: &str, value: i64) -> String {
    if value == 0 {
        return digits.chars().nth(1).unwrap().to_string();
    }

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
        '0'..='9' => i32::from(c as u8 - b'0'),
        'A'..='Z' => i32::from(c as u8 - b'A') + 10,
        'a'..='z' => i32::from(c as u8 - b'a') + 10,
        _ => panic!("Invalid character: {c}"),
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
    let trimmed_line = line.trim();

    let f = trimmed_line.chars().next();
    if let Some(f) = f {
        if f == ' ' {
            // Space-prefixed numbers are not encoded => return 0
            return Ok(0);
        } else if f == '-' || f.is_ascii_digit() {
            // Try to parse as a number
            return trimmed_line.parse().map_err(|_e| {
                CError::GenericError(format!("expected negative number. got '{f}'"))
            });
        }
    }

    if trimmed_line.is_empty() {
        return Ok(0);
    }

    if DIGITS_UPPER.contains(f.unwrap()) {
        let is_valid = trimmed_line
            .chars()
            .all(|c| c.is_ascii_digit() || c.is_ascii_uppercase());

        if !is_valid {
            return Err(CError::GenericError(format!(
                "the value '{trimmed_line}' is not a valid hybrid 36 number"
            )));
        }

        return Ok(decode_pure(trimmed_line) - 10 * pow_int(36, width - 1) + pow_int(10, width));
    }

    if DIGITS_LOWER.contains(f.unwrap()) {
        let is_valid = trimmed_line
            .chars()
            .all(|c| c.is_ascii_digit() || c.is_ascii_lowercase());

        if !is_valid {
            return Err(CError::GenericError(format!(
                "the value '{trimmed_line}' is not a valid hybrid 36 number"
            )));
        }

        return Ok(decode_pure(trimmed_line) + 16 * pow_int(36, width - 1) + pow_int(10, width));
    }

    Err(CError::GenericError(format!(
        "the value '{line}' is not a valid hybrid 36 number"
    )))
}

pub(crate) fn encode_hybrid36(width: usize, value: i64) -> String {
    // the number is too negative to be encoded
    if value < (1 - pow_int(10, width - 1)) {
        return "*".repeat(width);
    }

    // no need to encode
    if value < pow_int(10, width) {
        return value.to_string();
    }

    // use upper case set
    let mut val = value;
    val -= pow_int(10, width);
    if val < 26 * pow_int(36, width - 1) {
        val += 10 * pow_int(36, width - 1);
        return encode_pure(DIGITS_UPPER, val);
    }

    // use lower case set
    val -= 26 * pow_int(36, width - 1);
    if val < 26 * pow_int(36, width - 1) {
        val += 10 * pow_int(36, width - 1);
        return encode_pure(DIGITS_LOWER, val);
    }

    // too large to be encoded
    "*".repeat(width)
}

pub struct PDBFormat<'a> {
    /// Residue information in the current step
    pub residues: RefCell<Vec<(FullResidueId, Residue)>>,

    /// List of all atom offsets. This maybe pushed in `parse_atom` or if a TER
    /// record is found. It is reset every time a frame is read.
    pub atom_offsets: RefCell<Vec<usize>>,

    /// This will be None when no secondary structure information should be
    /// read. Else it is set to the final residue of a secondary structure and
    /// the text description which should be set.
    pub current_secinfo: RefCell<Option<(FullResidueId, &'a str)>>,

    /// Store secondary structure information. Keys are the
    /// starting residue of the secondary structure, and values are pairs
    /// containing the ending residue and a string which is a written
    /// description of the secondary structure
    pub secinfo: RefCell<BTreeMap<FullResidueId, (FullResidueId, &'a str)>>,

    /// Number of models read/written to the file
    models: RefCell<usize>,
}

impl<'a> PDBFormat<'a> {
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
            if initial_offset.is_err() {
                eprintln!("initial offset was err");
            }

            let unwrapped = initial_offset.unwrap();
            if unwrapped <= 0 {
                eprintln!("warning: '{unwrapped}' is too small, assuming id is '1'",);
                self.atom_offsets.borrow_mut().push(0);
            } else {
                self.atom_offsets.borrow_mut().push(
                    usize::try_from(unwrapped).expect("decode_hybrid36 returned a negative number")
                        - 1,
                );
            }
        }

        let name = &line[12..16].trim();
        let mut atom = if line.len() >= 78 {
            let atom_type = &line[76..78];
            let name = (*name).to_string();
            let symbol = atom_type.trim().to_string();
            Atom::with_symbol(name, symbol)
        } else {
            // Read just the atom name and hope for the best
            let name = name.to_string();
            Atom::new(name)
        };

        let altloc = &line[16..17];
        if altloc != " " {
            atom.properties
                .insert("altloc".into(), Property::String(altloc.to_string()));
        }

        let x = Property::parse_value(line[30..38].trim(), &PropertyKind::Double)?.expect_double();
        let y = Property::parse_value(line[38..46].trim(), &PropertyKind::Double)?.expect_double();
        let z = Property::parse_value(line[46..54].trim(), &PropertyKind::Double)?.expect_double();
        frame.add_atom(atom, [x, y, z]);

        let atom_id = frame.size() - 1;
        let insertion_code = line.chars().nth(26).unwrap();
        // If there's no residue information so return early
        let Ok(resid) = decode_hybrid36(4, &line[22..26]) else {
            return Ok(());
        };

        let chain = &line.chars().nth(21).unwrap();
        let resname = line[17..20].trim().to_string();
        let full_residue_id = FullResidueId {
            chain: *chain,
            resid,
            resname: resname.clone(),
            insertion_code,
        };

        let mut residues = self.residues.borrow_mut();

        // Find all matching residues
        let mut matching_positions = residues
            .iter()
            .enumerate()
            .filter(|(_, (id, _))| *id == full_residue_id)
            .map(|(pos, _)| pos)
            .peekable();

        if matching_positions.peek().is_some() {
            let val = matching_positions.next().unwrap();
            // Add atom to the first matching residue
            residues[val].1.add_atom(atom_id);
        } else {
            let mut residue = Residue {
                name: resname,
                id: Some(resid),
                ..Default::default()
            };
            residue.add_atom(atom_id);

            if insertion_code != ' ' {
                residue.properties.insert(
                    "insertion_code".into(),
                    Property::String(insertion_code.to_string()),
                );
            }

            // Set whether or not the residue is standardized
            residue
                .properties
                .insert("is_standard_pdb".into(), Property::Bool(!is_hetatm));

            // This is saved as a string (instead of a number) on purpose
            // to match MMTF format
            residue.properties.insert(
                "chainid".into(),
                Property::String(line.chars().nth(21).unwrap().to_string()),
            );

            // PDB format makes no distinction between chainid and chainname
            residue.properties.insert(
                "chainname".into(),
                Property::String(line.chars().nth(21).unwrap().to_string()),
            );

            // segment name is not part of the standard, but something added by
            // CHARM/NAMD in the un-used character range 73-76
            if line.len() > 72 {
                let segname = line[72..76].trim();
                if !segname.is_empty() {
                    residue
                        .properties
                        .insert("segname".into(), Property::String(segname.to_string()));
                }
            }

            // Are we within a secondary information sequence?
            // we add an extra scope so that the borrow of `self.current_secinfo` is dropped
            // before we access it further down
            {
                let current_secinfo_borrow = self.current_secinfo.borrow();
                if let Some((end_residue_id, description)) = current_secinfo_borrow.as_ref() {
                    residue.properties.insert(
                        "secondary_structure".into(),
                        Property::String(description.to_string()),
                    );

                    // Are we at the end of a secondary information sequence?
                    if *end_residue_id == full_residue_id {
                        // we drop the current borrow to allow the mutable borrow
                        drop(current_secinfo_borrow);
                        *self.current_secinfo.borrow_mut() = None;
                    }
                }
            }

            // Are we at the start of a secondary information sequence?
            // second borrow of `self.secinfo` happens here
            if let Some(secinfo_for_residue) = self.secinfo.borrow().get(&full_residue_id) {
                *self.current_secinfo.borrow_mut() = Some(secinfo_for_residue.clone());
                residue.properties.insert(
                    "secondary_structure".into(),
                    Property::String(secinfo_for_residue.1.to_string()),
                );
            }

            // Find the first spot where either
            // 1) the existing resid is greater than ours → we should come before it, or
            // 2) we hit our same resid → then we'll want to append after the last equal entry.
            let insert_pos = match residues
                .iter()
                .position(|(id, _)| id.resid >= full_residue_id.resid)
            {
                Some(pos) if residues[pos].0.resid == full_residue_id.resid => {
                    // we found an existing block of equal IDs; extend to its end
                    residues
                        .iter()
                        .rposition(|(id, _)| id.resid == full_residue_id.resid)
                        .unwrap()
                        + 1
                }
                Some(pos) => {
                    // the first resid ≥ ours but not equal → insert here
                    pos
                }
                None => {
                    // all existing resids are smaller → push to the end
                    residues.len()
                }
            };

            residues.insert(insert_pos, (full_residue_id, residue));
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

        let a = fast_float::parse(line[6..15].trim())
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let b = fast_float::parse(line[15..24].trim())
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let c = fast_float::parse(line[24..33].trim())
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let alpha = fast_float::parse(line[33..40].trim())
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let beta = fast_float::parse(line[40..47].trim())
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;
        let gamma = fast_float::parse(line[47..54].trim())
            .map_err(|e| CError::GenericError(format!("could not parse float: {e}")))?;

        let unit_cell = UnitCell::new_from_lengths_angles([a, b, c], &mut [alpha, beta, gamma])?;
        frame.unit_cell = unit_cell;

        if line.len() >= 55 {
            let space_group_slice = &line[55..std::cmp::min(line.len(), 65)];
            // TODO: handle this as a warning (somehow)?
            if !space_group_slice.contains("P 1") && !space_group_slice.contains("P1") {
                eprintln!(
                    "warning: ignoring custom space group ({space_group_slice}), using P1 instead"
                );
            }
        }

        Ok(())
    }
    fn read_index(&self, line: &str, initial: usize) -> Result<usize, CError> {
        debug_assert!(line.len() >= 5); // otherwise the indexing fails
        let mut pdb_atom_id = decode_hybrid36(5, &line[initial..initial + 5])?;

        // Find the lower bound index
        let lower = self
            .atom_offsets
            .borrow()
            .binary_search(&(pdb_atom_id as usize))
            .unwrap_or_else(|insert_pos| insert_pos)
            - 1;
        pdb_atom_id -= lower as i64;

        // TODO: is this correct?
        pdb_atom_id -= 1;
        Ok(usize::try_from(pdb_atom_id).unwrap() - self.atom_offsets.borrow().first().unwrap())
    }

    fn add_bond(frame: &mut Frame, line: &str, i: usize, j: usize) {
        if i >= frame.size() || j >= frame.size() {
            eprintln!(
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

    fn parse_conect(&self, frame: &mut Frame, line: &str) {
        debug_assert_eq!(&line[..6], "CONECT");

        let line_length = line.trim().len();

        let index_i = self
            .read_index(line, 6)
            .expect("could not read index when parsing CONECT");

        if line_length > 11 {
            let index_j = self.read_index(line, 11).unwrap();
            PDBFormat::add_bond(frame, line, index_i, index_j);
        }

        if line_length > 16 {
            let index_j = self.read_index(line, 16).unwrap();
            PDBFormat::add_bond(frame, line, index_i, index_j);
        }

        if line_length > 21 {
            let index_j = self.read_index(line, 21).unwrap();
            PDBFormat::add_bond(frame, line, index_i, index_j);
        }

        if line_length > 26 {
            let index_j = self.read_index(line, 26).unwrap();
            PDBFormat::add_bond(frame, line, index_i, index_j);
        }
    }

    fn parse_helix(&self, line: &str) -> Result<(), CError> {
        if line.len() < 33 + 5 {
            eprintln!("warning: HELIX record too short: {line}");
        }

        let chain_start = line.chars().nth(19).expect("start of chain");
        let chain_end = line.chars().nth(31).expect("end of chain");
        let inscode_start = line.chars().nth(25).expect("start of insertion code");
        let inscode_end = line.chars().nth(37).expect("end of insertion code");
        let resname_start = line[15..18].trim();
        let resname_end = line[27..30].trim();

        let resid_start = decode_hybrid36(4, &line[21..25])?;
        let resid_end = decode_hybrid36(4, &line[33..37])?;

        if chain_start != chain_end {
            return Err(CError::GenericError(format!(
                "warning: HELIX chain {chain_start} and {chain_end} are not the same"
            )));
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
            .trim()
            .parse::<usize>()
            .inspect_err(|e| eprintln!("failed to parse helix type: {e}"))
            .unwrap();

        if *helix_type <= 10 {
            self.secinfo.borrow_mut().insert(
                start,
                (
                    end,
                    HelixType::nth(helix_type - 1).expect("unknown helix type"),
                ),
            );
        }

        Ok(())
    }

    fn parse_secondary(&self, line: &str, start: usize, end: usize) -> Result<(), CError> {
        if line.len() < end + 10 {
            eprintln!("warning: secondary structure record too short: '{line}'");
        }

        let resname_start = &line[start..start + 3].trim();
        let resname_end = &line[end..end + 3].trim();
        let chain_start = &line.chars().nth(start + 4).expect("start of chain");
        let chain_end = &line.chars().nth(end + 4).expect("end of chain");

        if chain_start != chain_end {
            return Err(CError::GenericError(format!(
                "warning: secondary chain {chain_start} and {chain_end} are not the same"
            )));
        }

        let resid_start = decode_hybrid36(4, &line[start..start + 4])?;
        let resid_end = decode_hybrid36(4, &line[end..end + 4])?;
        let inscode_start = &line
            .chars()
            .nth(start + 9)
            .expect("start of insertion code");
        let inscode_end = &line.chars().nth(end + 9).expect("end of insertion code");

        let start = FullResidueId {
            chain: *chain_start,
            resid: resid_start,
            resname: (*resname_start).to_string(),
            insertion_code: *inscode_start,
        };

        let end = FullResidueId {
            chain: *chain_end,
            resid: resid_end,
            resname: (*resname_end).to_string(),
            insertion_code: *inscode_end,
        };

        self.secinfo.borrow_mut().insert(start, (end, "extended"));

        Ok(())
    }

    pub fn new() -> Self {
        PDBFormat {
            residues: RefCell::new(Vec::new()),
            atom_offsets: RefCell::new(Vec::new()),
            current_secinfo: RefCell::new(None),
            secinfo: RefCell::new(BTreeMap::new()),
            models: RefCell::new(0),
        }
    }

    fn chain_ended(&self, frame: &mut Frame) {
        // We drain as a 'hack' to allow for badly formatted PDB files which restart
        // the residue ID after a TER residue in cases where they should not.
        // IE a metal Ion given the chain ID of A and residue ID of 1 even though
        // this residue already exists.
        let mut residues = self.residues.borrow_mut();
        for (_, residue) in residues.drain(..) {
            let _ = frame.add_residue(residue);
        }
    }

    const END_RECORD: &'a str = "END";
    const ENDMDL_RECORD: &'a str = "ENDMDL";

    fn link_standard_residue_bonds(frame: &mut Frame) {
        let mut link_previous_peptide = false;
        let mut link_previous_nucleic = false;
        let mut previous_residue_id = 0;
        let mut previous_carboxylic_id = 0;

        let residues = &frame.topology().residues;

        // Collect all the bonds we need to add first
        // We need at least `residues.len()` bonds to add + fudge factor
        let mut bonds_to_add = Vec::with_capacity(residues.len() * 3);

        for residue in residues {
            let residue_table = find(&residue.name);
            if residue_table.is_none() {
                continue;
            }

            let mut atom_name_to_index = HashMap::<String, usize>::new();
            for &atom in residue {
                atom_name_to_index.insert(frame[atom].name.to_string(), atom);
            }

            let amide_nitrogen = atom_name_to_index.get("N");
            let amide_carbon = atom_name_to_index.get("C");

            // Check if the residue has an ID
            let Some(resid) = residue.id else {
                eprintln!("warning: got a residue without id, this should not happen");
                return;
            };

            if link_previous_peptide {
                if let Some(&n_index) = amide_nitrogen {
                    if resid == previous_residue_id + 1 {
                        link_previous_peptide = false;
                        bonds_to_add.push((previous_carboxylic_id, n_index));
                    }
                }
            }

            if let Some(&ac_index) = amide_carbon {
                link_previous_peptide = true;
                previous_carboxylic_id = ac_index;
                previous_residue_id = resid;
            }

            let three_prime_oxygen = atom_name_to_index.get("O3'");
            let five_prime_phosphorus = atom_name_to_index.get("P");
            if link_previous_nucleic {
                if let Some(&_) = five_prime_phosphorus {
                    if resid == previous_residue_id + 1 {
                        link_previous_nucleic = false;
                        if let Some(&o_index) = three_prime_oxygen {
                            bonds_to_add.push((previous_carboxylic_id, o_index));
                        }
                    }
                }
            }

            if let Some(&o_index) = three_prime_oxygen {
                link_previous_nucleic = true;
                previous_carboxylic_id = o_index;
                previous_residue_id = resid;
            }

            // A special case missed by the standards committee??
            if atom_name_to_index.contains_key("HO5'") {
                if let (Some(&ho5), Some(&o5)) = (
                    atom_name_to_index.get("HO5'"),
                    atom_name_to_index.get("O5'"),
                ) {
                    bonds_to_add.push((ho5, o5));
                }
            }

            for link in residue_table.unwrap() {
                let first_atom = atom_name_to_index.get(pdb_connectivity::INTERNER[link.0]);
                let second_atom = atom_name_to_index.get(pdb_connectivity::INTERNER[link.1]);

                if first_atom.is_none() {
                    let first_name = pdb_connectivity::INTERNER[link.0];
                    let mut chars = first_name.chars();
                    let first_char = chars.next();
                    let second_char = chars.next();

                    if first_char != Some('H')
                        && first_name != "OXT"
                        && first_char != Some('P')
                        && first_char != Some('O')
                        && second_char != Some('P')
                    {
                        eprintln!(
                            "warning: could not find standard atom '{first_name}' in residue '{}' (resid {resid})",
                            residue.name
                        );
                    }
                    continue;
                }

                if second_atom.is_none() {
                    let second_name = pdb_connectivity::INTERNER[link.1];
                    let mut chars = second_name.chars();
                    let first_char = chars.next();
                    let second_char = chars.next();

                    if first_char != Some('H')
                        && second_name != "OXT"
                        && first_char != Some('P')
                        && first_char != Some('O')
                        && second_char != Some('P')
                    {
                        eprintln!(
                            "warning: could not find standard atom '{second_name}' in residue '{}' (resid {resid})",
                            residue.name
                        );
                    }
                    continue;
                }

                if let (Some(&first), Some(&second)) = (first_atom, second_atom) {
                    bonds_to_add.push((first, second));
                }
            }
        }
        // Now add all the bonds at once
        for (i, j) in bonds_to_add {
            frame
                .add_bond(i, j, BondOrder::Unknown)
                .expect("unable to add bond to frame");
        }
    }

    // Check the number of digits before the decimal separator to be sure than we
    // can represent them. In case of error, use the given `context` in the error
    // message
    fn check_values_size(values: [f64; 3], width: i32, context: &str) -> Result<(), CError> {
        let max_pos = f64::powi(10.0, width) - 1.0;
        let max_neg = -f64::powi(10.0, width) - 1.0;

        if values[0] > max_pos
            || values[1] > max_pos
            || values[2] > max_pos
            || values[0] < max_neg
            || values[1] < max_neg
            || values[2] < max_neg
        {
            return Err(CError::GenericError(format!(
                "value in {context} is too big for representation in PDB format"
            )));
        }

        Ok(())
    }

    fn to_pdb_index(value: i64, width: usize) -> String {
        let encoded = encode_hybrid36(width, value + 1);
        if encoded.chars().nth(0).unwrap() == '*'
            && (value == MAX_HYBRID36_W4_NUMBER || value == MAX_HYBRID36_W5_NUMBER)
        {
            let t = if width == 5 { "atom" } else { "residue" };
            eprintln!(
                "warning: the value for a {t} serial/id is too large, using '{encoded}' instead"
            );
        }

        encoded
    }

    fn get_residue_information(
        residue: Option<&Residue>,
        max_resid: &mut i64,
    ) -> ResidueInformation {
        let mut info = ResidueInformation::new();

        if residue.is_none() {
            let val = *max_resid;
            *max_resid += 1;
            info.resid = PDBFormat::to_pdb_index(val, 4);
            return info;
        }

        let residue = residue.as_ref().unwrap();

        if residue
            .get("is_standard_pdb")
            .and_then(Property::as_bool)
            .unwrap_or(false)
        {
            info.atom_hetatm = "ATOM  ".to_string();
        }

        info.resname = residue.name.clone();
        if info.resname.len() > 3 {
            eprint!(
                "warning: residue '{}' name is too long, it will be truncated",
                info.resname
            );
            info.resname = info.resname.chars().take(3).collect();
        }

        if residue.id.is_some() {
            info.resid = PDBFormat::to_pdb_index(residue.id.unwrap() - 1, 4);
        }

        info.chainid = residue
            .get("chainid")
            .and_then(Property::as_string)
            .unwrap_or(" ")
            .to_string();
        if info.chainid.len() > 1 {
            eprintln!(
                "warning: residues's chain id '{}' is too long, it will be trunctated",
                info.chainid
            );
            info.chainid = info.chainid.chars().nth(0).unwrap().to_string();
        }

        info.insertion_code = residue
            .get("insertion_code")
            .and_then(Property::as_string)
            .unwrap_or("")
            .to_string();
        if info.insertion_code.len() > 1 {
            eprintln!(
                "warning: residue's insertion code '{}' is too long, it will be truncated",
                info.insertion_code
            );
            info.insertion_code = info.insertion_code.chars().nth(0).unwrap().to_string();
        }

        info.segment = residue
            .get("segname")
            .and_then(Property::as_string)
            .unwrap_or("")
            .to_string();
        if info.segment.len() > 4 {
            eprintln!(
                "residue's segment name '{}' is too long, it will be truncated",
                info.segment
            );
            info.segment = info.segment.chars().take(4).collect();
        }

        info.composition_type = residue
            .get("composition_type")
            .and_then(Property::as_string)
            .unwrap_or("")
            .to_string();

        info
    }

    // This function adjusts a given index to account for intervening TER records.
    // It does so by determining the position of the greatest TER record in `ters`
    // and uses iterator arithmetic to calculate the adjustment. Note that `ters`
    // is expected to be sorted
    fn adjust_for_ter_residues(v: usize, ters: &[usize]) -> i64 {
        let insertion_point = ters
            .binary_search_by(|&ter| {
                if ter < v + 1 {
                    std::cmp::Ordering::Less
                } else {
                    std::cmp::Ordering::Greater
                }
            })
            .unwrap_err();
        (v + insertion_point) as i64
    }
}

impl FileFormat for PDBFormat<'_> {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        self.residues.borrow_mut().clear();
        self.atom_offsets.borrow_mut().clear();

        let mut frame = Frame::new();
        let mut line = String::new();

        let mut got_end = false;
        while !got_end && reader.read_line(&mut line)? > 0 {
            let record = get_record(&line);

            match record {
                Record::HEADER => {
                    if line.len() >= 50 {
                        frame.properties.insert(
                            "classification".into(),
                            Property::String(line[10..50].trim().to_string()),
                        );
                    }
                    if line.len() >= 59 {
                        frame.properties.insert(
                            "deposition_date".into(),
                            Property::String(line[50..59].trim().to_string()),
                        );
                    }
                    if line.len() >= 66 {
                        frame.properties.insert(
                            "pdb_idcode".into(),
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
                        format!("{current} {title}")
                    };
                    frame
                        .properties
                        .insert("name".into(), Property::String(new_title));
                }
                Record::CRYST1 => PDBFormat::parse_cryst1(&mut frame, &line).unwrap(),
                Record::ATOM => self.parse_atom(&mut frame, &line, false).unwrap(),
                Record::HETATM => self.parse_atom(&mut frame, &line, true).unwrap(),
                Record::CONECT => self.parse_conect(&mut frame, &line),
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
                    // Else we have read a frame
                    got_end = true;
                }
                Record::HELIX => self.parse_helix(&line).unwrap(),
                Record::SHEET => self.parse_secondary(&line, 17, 28).unwrap(),
                Record::TURN => self.parse_secondary(&line, 15, 26).unwrap(),
                Record::TER => {
                    if line.len() >= 12 {
                        let ter_serial =
                            decode_hybrid36(5, &line[6..11]).expect("TER record not numeric");
                        if ter_serial != 0 {
                            // This happens if the TER serial number is blank
                            self.atom_offsets.borrow_mut().push(
                                usize::try_from(ter_serial)
                                    .expect("could not parse ter_serial to usize"),
                            );
                        }
                    }
                    self.chain_ended(&mut frame);
                }
                Record::END => got_end = true,
                Record::IGNORED_ => {}
                Record::UNKNOWN_ => {
                    eprintln!("ignoring unknown record: {line}");
                }
            }
            line.clear();
        }

        if !got_end {
            eprintln!("warning: missing END record in file");
        }

        self.chain_ended(&mut frame);
        PDBFormat::link_standard_residue_bonds(&mut frame);
        Ok(frame)
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        let mut line = String::new();
        let position = reader.stream_position()?;

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
                return Ok(Some(position));
            }

            line.clear();
        }

        if position == 0 {
            Ok(Some(position))
        } else {
            Ok(None)
        }
    }

    fn write_next(&self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
        writeln!(writer, "MODEL {:>4}", *self.models.borrow() + 1)?;

        let lengths = frame.unit_cell.lengths();
        let angles = frame.unit_cell.angles();

        PDBFormat::check_values_size(lengths, 9, "cell lengths")?;
        PDBFormat::check_values_size(angles, 7, "cell angles")?;
        // Do not try to guess the space group and the z value, just use the
        // default one.
        writeln!(
            writer,
            "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} P 1           1",
            lengths[0], lengths[1], lengths[2], angles[0], angles[1], angles[2]
        )?;

        // Only use numbers bigger than the biggest residue id as "resSeq" for
        // atoms without associated residue.
        let mut max_resid = 0;
        for residue in &frame.topology().residues {
            let resid = residue.id;
            if resid.is_some() && resid.unwrap() > max_resid {
                max_resid = resid.unwrap();
            }
        }

        // Used to skip writing unnecessary connect record
        // use std::deque because std::vector<bool> have a surprising behavior due
        // to C++ standard requiring it to pack multiple bool in a byte
        let mut is_atom_record = vec![false; frame.len()];

        // Used for writing TER records.
        let mut ter_count = 0;
        let mut last_residue: Option<ResidueInformation> = None;
        let mut ter_serial_numbers: Vec<usize> = vec![];

        for (idx, pos) in frame.positions().iter().enumerate() {
            let mut altloc = frame[idx]
                .properties
                .get("altloc")
                .and_then(|prop| prop.as_string())
                .unwrap_or(" ")
                .to_string();
            if altloc.len() > 1 {
                eprintln!("warning: altloc '{altloc}' is too long, it will be truncated");
                altloc = altloc.chars().next().unwrap().to_string();
            }

            let residue = frame.topology().residue_for_atom(idx);
            let resinfo = PDBFormat::get_residue_information(residue.as_ref(), &mut max_resid);

            if resinfo.atom_hetatm == "ATOM  " {
                is_atom_record[idx] = true;
            }

            debug_assert!(resinfo.resname.len() <= 3);

            if last_residue.is_some()
                && last_residue.as_ref().unwrap().chainid != resinfo.chainid
                && last_residue.as_ref().unwrap().needs_ter_record()
            {
                let pdb_index = PDBFormat::to_pdb_index(
                    i64::try_from(idx + ter_count).expect("idx + ter_count fits in i64"),
                    5,
                );
                let unwrapped = last_residue.as_ref().unwrap();
                writeln!(
                    writer,
                    "TER   {:>5}      {:3} {}{:>4}{}",
                    pdb_index,
                    unwrapped.resname,
                    unwrapped.chainid,
                    unwrapped.resid,
                    unwrapped.insertion_code
                )?;
                ter_serial_numbers.push(idx + ter_count);
                ter_count += 1;
            }

            PDBFormat::check_values_size(*pos, 8, "atomic position")?;
            writeln!(
                writer,
                "{:<6}{:>5} {:<4}{:1}{:3} {:1}{:>4}{:1}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}      {: <4}{: >2}",
                resinfo.atom_hetatm,
                PDBFormat::to_pdb_index(
                    i64::try_from(idx + ter_count).expect("convert idx + ter_count"),
                    5
                ),
                frame[idx].name,
                altloc,
                resinfo.resname,
                resinfo.chainid,
                resinfo.resid,
                // 444,
                resinfo.insertion_code,
                pos[0],
                pos[1],
                pos[2],
                1.0,
                0.0,
                resinfo.segment,
                frame[idx].symbol
            )?;
            if residue.as_ref().is_some() {
                last_residue = Some(resinfo);
            } else {
                last_residue = None;
            }
        }

        let mut connect = vec![Vec::new(); frame.size()];
        for bond in frame.topology().bonds() {
            if is_atom_record[bond[0]] && is_atom_record[bond[1]] {
                // both must be standard residue atoms
                continue;
            }
            if bond[0] > MAX_HYBRID36_W5_NUMBER as usize
                || bond[1] > MAX_HYBRID36_W5_NUMBER as usize
            {
                eprintln!(
                    "warning: atomic index is too big for CONNECT, removing the bond between {} and {}",
                    bond[0], bond[1]
                );
            }

            connect[bond[0]].push(PDBFormat::adjust_for_ter_residues(
                bond[1],
                &ter_serial_numbers,
            ));
            connect[bond[1]].push(PDBFormat::adjust_for_ter_residues(
                bond[0],
                &ter_serial_numbers,
            ));
        }

        for (idx, c) in connect.iter().enumerate().take(frame.size()) {
            let connections = c.len();
            let lines = connections / 4 + 1;
            if connections == 0 {
                continue;
            }

            let correction = PDBFormat::adjust_for_ter_residues(idx, &ter_serial_numbers);

            for conect_line in 0..lines {
                write!(
                    writer,
                    "CONECT{:>5}",
                    PDBFormat::to_pdb_index(correction, 5)
                )?;

                let last = std::cmp::min(connections, 4 * (conect_line + 1));
                for c_ij in c.iter().take(last).skip(4 * conect_line) {
                    write!(writer, "{:>5}", PDBFormat::to_pdb_index(*c_ij, 5))?;
                }
                writeln!(writer)?;
            }
        }

        writeln!(writer, "ENDMDL")?;
        let current_models = *self.models.borrow();
        *self.models.borrow_mut() = current_models + 1;

        Ok(())
    }

    fn read(&self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        if reader.fill_buf().map(|b| !b.is_empty()).unwrap() {
            Ok(Some(self.read_next(reader).unwrap()))
        } else {
            Ok(None)
        }
    }
    /// Finalize the PDB file by writing the END record if needed
    ///
    /// This should be called when done writing to a PDB file to ensure it's properly closed
    fn finalize(&self, writer: &mut BufWriter<File>) -> Result<(), CError> {
        if *self.models.borrow() > 0 {
            writeln!(writer, "END")?;
        }
        Ok(())
    }
}

#[derive(Default)]
struct ResidueInformation {
    atom_hetatm: String,
    resname: String,
    resid: String,
    chainid: String,
    insertion_code: String,
    composition_type: String,
    segment: String,
}

impl ResidueInformation {
    pub fn new() -> Self {
        Self {
            atom_hetatm: "HETATM".to_string(),
            ..Default::default()
        }
    }

    pub fn needs_ter_record(&self) -> bool {
        let ct = &self.composition_type;
        // If it's empty, or matches "other"/"non-polymer" (in any ASCII case), we skip.
        let skip = ["other", "non-polymer"];
        !(ct.is_empty() || skip.iter().any(|&s| ct.eq_ignore_ascii_case(s)))
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use assert_approx_eq::assert_approx_eq;
    use tempfile::Builder;

    use crate::{
        angle::Angle,
        atom::Atom,
        bond::{Bond, BondOrder},
        dihedral::Dihedral,
        formats::pdb::{decode_hybrid36, encode_hybrid36},
        frame::Frame,
        property::Property,
        residue::Residue,
        topology::Topology,
        trajectory::Trajectory,
        unit_cell::{CellShape, UnitCell},
    };

    use std::borrow::ToOwned;

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/pdb/water.pdb");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 100);

        let path = Path::new("./src/tests-data/pdb/2hkb.pdb");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 11);
    }

    #[test]
    fn sanity_check() {
        let path = Path::new("./src/tests-data/pdb/water.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
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
        let mut trajectory = Trajectory::open(path).unwrap();
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

        assert!(topology.angles().contains(&Angle::new(20, 21, 23)));
        assert!(topology.angles().contains(&Angle::new(9, 38, 44)));

        assert!(topology
            .dihedrals()
            .contains(&Dihedral::new(64, 62, 58, 53)));
        assert!(topology
            .dihedrals()
            .contains(&Dihedral::new(22, 21, 23, 33)));
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
        recycle_check!(4, 10_000, "A000");
        recycle_check!(4, 10_001, "A001");
        recycle_check!(4, 10_002, "A002");
        recycle_check!(4, 10_003, "A003");
        recycle_check!(4, 10_004, "A004");
        recycle_check!(4, 10_005, "A005");
        recycle_check!(4, 10_006, "A006");
        recycle_check!(4, 10_007, "A007");
        recycle_check!(4, 10_008, "A008");
        recycle_check!(4, 10_009, "A009");
        recycle_check!(4, 10_010, "A00A");
        recycle_check!(4, 10_011, "A00B");
        recycle_check!(4, 10_012, "A00C");
        recycle_check!(4, 10_013, "A00D");
        recycle_check!(4, 10_014, "A00E");
        recycle_check!(4, 10_015, "A00F");
        recycle_check!(4, 10_016, "A00G");
        recycle_check!(4, 10_017, "A00H");
        recycle_check!(4, 10_018, "A00I");
        recycle_check!(4, 10_019, "A00J");
        recycle_check!(4, 10_020, "A00K");
        recycle_check!(4, 10_021, "A00L");
        recycle_check!(4, 10_022, "A00M");
        recycle_check!(4, 10_023, "A00N");
        recycle_check!(4, 10_024, "A00O");
        recycle_check!(4, 10_025, "A00P");
        recycle_check!(4, 10_026, "A00Q");
        recycle_check!(4, 10_027, "A00R");
        recycle_check!(4, 10_028, "A00S");
        recycle_check!(4, 10_029, "A00T");
        recycle_check!(4, 10_030, "A00U");
        recycle_check!(4, 10_031, "A00V");
        recycle_check!(4, 10_032, "A00W");
        recycle_check!(4, 10_033, "A00X");
        recycle_check!(4, 10_034, "A00Y");
        recycle_check!(4, 10_035, "A00Z");
        recycle_check!(4, 10_036, "A010");
        recycle_check!(4, 10_046, "A01A");
        recycle_check!(4, 10_071, "A01Z");
        recycle_check!(4, 10_072, "A020");
        recycle_check!(4, 10_000 + 36 * 36 - 1, "A0ZZ");
        recycle_check!(4, 10_000 + 36 * 36, "A100");
        recycle_check!(4, 10_000 + 36 * 36 * 36 - 1, "AZZZ");
        recycle_check!(4, 10_000 + 36 * 36 * 36, "B000");
        recycle_check!(4, 10_000 + 26 * 36 * 36 * 36 - 1, "ZZZZ");
        recycle_check!(4, 10_000 + 26 * 36 * 36 * 36, "a000");
        recycle_check!(4, 10_000 + 26 * 36 * 36 * 36 + 35, "a00z");
        recycle_check!(4, 10_000 + 26 * 36 * 36 * 36 + 36, "a010");
        recycle_check!(4, 10_000 + 26 * 36 * 36 * 36 + 36 * 36 - 1, "a0zz");
        recycle_check!(4, 10_000 + 26 * 36 * 36 * 36 + 36 * 36, "a100");
        recycle_check!(4, 10_000 + 26 * 36 * 36 * 36 + 36 * 36 * 36 - 1, "azzz");
        recycle_check!(4, 10_000 + 26 * 36 * 36 * 36 + 36 * 36 * 36, "b000");
        recycle_check!(4, 10_000 + 2 * 26 * 36 * 36 * 36 - 1, "zzzz");

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
        recycle_check!(5, 99_999, "99999");
        recycle_check!(5, 100_000, "A0000");
        recycle_check!(5, 100_010, "A000A");
        recycle_check!(5, 100_035, "A000Z");
        recycle_check!(5, 100_036, "A0010");
        recycle_check!(5, 100_046, "A001A");
        recycle_check!(5, 100_071, "A001Z");
        recycle_check!(5, 100_072, "A0020");
        recycle_check!(5, 100_000 + 36 * 36 - 1, "A00ZZ");
        recycle_check!(5, 100_000 + 36 * 36, "A0100");
        recycle_check!(5, 100_000 + 36 * 36 * 36 - 1, "A0ZZZ");
        recycle_check!(5, 100_000 + 36 * 36 * 36, "A1000");
        recycle_check!(5, 100_000 + 36 * 36 * 36 * 36 - 1, "AZZZZ");
        recycle_check!(5, 100_000 + 36 * 36 * 36 * 36, "B0000");
        recycle_check!(5, 100_000 + 2 * 36 * 36 * 36 * 36, "C0000");
        recycle_check!(5, 100_000 + 26 * 36 * 36 * 36 * 36 - 1, "ZZZZZ");
        recycle_check!(5, 100_000 + 26 * 36 * 36 * 36 * 36, "a0000");
        recycle_check!(5, 100_000 + 26 * 36 * 36 * 36 * 36 + 36 - 1, "a000z");
        recycle_check!(5, 100_000 + 26 * 36 * 36 * 36 * 36 + 36, "a0010");
        recycle_check!(5, 100_000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 - 1, "a00zz");
        recycle_check!(5, 100_000 + 26 * 36 * 36 * 36 * 36 + 36 * 36, "a0100");
        recycle_check!(
            5,
            100_000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36 - 1,
            "a0zzz"
        );
        recycle_check!(5, 100_000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36, "a1000");
        recycle_check!(
            5,
            100_000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36 * 36 - 1,
            "azzzz"
        );
        recycle_check!(
            5,
            100_000 + 26 * 36 * 36 * 36 * 36 + 36 * 36 * 36 * 36,
            "b0000"
        );
        recycle_check!(5, 100_000 + 2 * 26 * 36 * 36 * 36 * 36 - 1, "zzzzz");

        assert_eq!(encode_hybrid36(4, -99_999), "****");
        assert_eq!(encode_hybrid36(4, 9_999_999), "****");
    }

    #[test]
    fn read_residue_information() {
        let path = Path::new("./src/tests-data/pdb/water.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        assert_eq!(frame.topology().residues.len(), 99);
        let opt_residue = frame.topology().residue_for_atom(1);
        assert!(opt_residue.is_some());
        let residue = opt_residue.unwrap();

        assert_eq!(residue.size(), 3);
        assert!(residue.contains(0));
        assert!(residue.contains(1));
        assert!(residue.contains(2));
        assert!(residue.get("chainid").is_some());
        assert_eq!(
            residue
                .get("chainid")
                .map(ToOwned::to_owned)
                .unwrap()
                .expect_string(),
            "X"
        );

        let path = Path::new("./src/tests-data/pdb/MOF-5.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        assert_eq!(frame.topology().residues.len(), 1);
        let residue = frame.topology().residues[0].clone();
        assert_eq!(residue.size(), frame.size());
        assert_eq!(residue.name, "LIG");
    }

    #[test]
    fn read_atom_hetatm_information() {
        let path = Path::new("./src/tests-data/pdb/hemo.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();
        let residues = frame.topology().residues.clone();

        assert!(!residues[0].get("is_standard_pdb").unwrap().expect_bool());

        for res in &residues[1..] {
            assert!(res.get("is_standard_pdb").unwrap().expect_bool());
        }
        assert_eq!(frame[74].symbol, "C");
    }

    #[test]
    fn handle_multiple_ter_records() {
        let path = Path::new("./src/tests-data/pdb/4hhb.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        assert_eq!(frame[4556].name, "ND");
        assert_eq!(frame[4557].name, "FE");
        assert_eq!(
            frame.topology().bond_order(4556, 4557).unwrap(),
            BondOrder::Unknown
        );

        let topology = frame.topology();
        assert_eq!(
            topology.residues[5]
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed alpha helix"
        );
    }

    #[test]
    fn secondary_structure_with_insertion_code_test() {
        let path = Path::new("./src/tests-data/pdb/1bcu.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        // check that residues have been inserted correctly
        let topology = frame.topology();
        assert_eq!(topology.residue_for_atom(0).unwrap().name, "ALA");
        assert_eq!(
            topology
                .residue_for_atom(0)
                .unwrap()
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "B"
        );
        assert_eq!(
            topology
                .residue_for_atom(5)
                .unwrap()
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "A"
        );
        assert!(topology
            .residue_for_atom(13)
            .unwrap()
            .get("insertion_code")
            .is_none(),);

        // // Check secondary structure, no insertion code
        assert_eq!(
            topology.residues[9]
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );
        assert_eq!(
            topology.residues[10]
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );
        assert_eq!(
            topology.residues[11]
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );
        assert!(topology.residues[12].get("secondary_structure").is_none());
        assert!(topology.residues[13].get("secondary_structure").is_none());
        assert!(topology.residues[14].get("secondary_structure").is_none());
        assert!(topology.residues[15].get("secondary_structure").is_none());
        assert!(topology.residues[16].get("secondary_structure").is_none());
        assert!(topology.residues[17].get("secondary_structure").is_none());

        // First residue in a long list of residues with the same secondary structure
        let ins_check = topology.residues[18].clone();
        assert_eq!(
            ins_check
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed alpha helix"
        );
        assert_eq!(
            ins_check.get("insertion_code").unwrap().expect_string(),
            "C"
        );
        assert_eq!(ins_check.id.unwrap(), 14);
        assert_eq!(ins_check.get("chainid").unwrap().expect_string(), "L");

        assert_eq!(
            topology.residues[19]
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed alpha helix"
        );
        assert_eq!(
            topology.residues[19]
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "D"
        );
        assert_eq!(
            topology.residues[20]
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed alpha helix"
        );
        assert_eq!(
            topology.residues[20]
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "E"
        );
        assert_eq!(
            topology.residues[21]
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed alpha helix"
        );
        assert_eq!(
            topology.residues[21]
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "F"
        );
        assert_eq!(
            topology.residues[22]
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed alpha helix"
        );
        assert_eq!(
            topology.residues[22]
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "G"
        );

        assert!(topology.residues[23].get("secondary_structure").is_none());
        assert_eq!(
            topology.residues[23]
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "H"
        );
        assert_eq!(topology.residues[23].id.unwrap(), 14);
        assert_eq!(
            topology.residues[23]
                .get("chainid")
                .unwrap()
                .expect_string(),
            "L"
        );
    }

    #[test]
    fn read_protein_residues() {
        let path = Path::new("./src/tests-data/pdb/hemo.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        let topology = frame.topology();

        assert!(!topology.are_linked(&topology.residues[2], &topology.residues[3]));
        assert!(topology.are_linked(&topology.residues[3], &topology.residues[3]));
        assert!(topology.are_linked(&topology.residues[3], &topology.residues[4]));
        assert!(!topology.are_linked(&topology.residues[3], &topology.residues[5]));
        assert_eq!(topology.bonds().len(), 482);
    }

    #[test]
    fn read_nucleic_residues() {
        let path = Path::new("./src/tests-data/pdb/2hkb.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        let topology = frame.topology();

        assert!(topology.are_linked(&topology.residues[3], &topology.residues[4]));
        assert!(!topology.are_linked(&topology.residues[3], &topology.residues[5]));
        assert_eq!(topology.bonds().len(), 815);
    }

    #[test]
    fn read_atomic_insertion_codes() {
        let path = Path::new("./src/tests-data/pdb/insertion-code.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        let topology = frame.topology();

        assert_eq!(
            topology
                .residue_for_atom(0)
                .unwrap()
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "a"
        );
        assert_eq!(
            topology
                .residue_for_atom(1)
                .unwrap()
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "c"
        );
        assert_eq!(
            topology
                .residue_for_atom(2)
                .unwrap()
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "x"
        );
        assert!(frame[3].properties.get("insertion_code").is_none());
    }

    #[test]
    fn multiple_residues_with_the_same_id() {
        let path = Path::new("./src/tests-data/pdb/psfgen-output.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        let topology = frame.topology();

        assert_eq!(topology.residues.len(), 3);
        assert_eq!(topology.residues[0].name, "ALA");
        assert_eq!(topology.residues[0].id.unwrap(), 1);
        assert_eq!(
            topology.residues[0].get("segname").unwrap().expect_string(),
            "PROT"
        );

        assert_eq!(topology.residues[1].name, "GLY");
        assert_eq!(topology.residues[1].id.unwrap(), 1);
        assert_eq!(
            topology.residues[1].get("segname").unwrap().expect_string(),
            "PROT"
        );

        assert_eq!(topology.residues[2].name, "GLY");
        assert_eq!(topology.residues[2].id.unwrap(), 2);
        assert_eq!(
            topology.residues[2].get("segname").unwrap().expect_string(),
            "PROT"
        );
    }

    #[test]
    fn odd_pdb_numbering() {
        let path = Path::new("./src/tests-data/pdb/odd-start.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        assert_eq!(frame.size(), 20);
        assert_eq!(frame[0].name, "C1");
        assert_eq!(frame[19].name, "C18");
        assert_eq!(
            frame.topology().bond_order(0, 1).unwrap(),
            BondOrder::Unknown
        );
        assert_eq!(
            frame.topology().bond_order(19, 13).unwrap(),
            BondOrder::Unknown
        );
    }

    #[test]
    fn atom_id_starts_at_0() {
        let path = Path::new("./src/tests-data/pdb/atom-id-0.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 1);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 2);

        assert_eq!(frame[0].name, "C1");
        assert_eq!(frame[1].name, "C2");
        assert_eq!(frame[0].symbol, "C");
        assert_eq!(frame[0].symbol, "C");
    }

    #[test]
    fn multiple_end_records() {
        let path = Path::new("./src/tests-data/pdb/end-endmdl.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 2);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 4);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 7);
    }

    #[test]
    fn multiple_model_without_end() {
        let path = Path::new("./src/tests-data/pdb/model.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 2);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 2223);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 2223);
    }

    #[test]
    fn file_generated_by_crystal_maker() {
        let path = Path::new("./src/tests-data/pdb/crystal-maker.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 1);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 8);
    }

    #[test]
    fn short_cryst1_record() {
        let path = Path::new("./src/tests-data/pdb/short-cryst1.pdb");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 1);
    }

    #[test]
    fn short_atom_record() {
        let path = Path::new("./src/tests-data/pdb/short-atom.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 9);

        assert_eq!(frame[0].name, "O");
        assert_eq!(frame[5].name, "H");
        assert_eq!(frame[0].symbol, "O");
        assert_eq!(frame[5].symbol, "H");

        assert_approx_eq!(frame.positions()[0][0], 0.417);
        assert_approx_eq!(frame.positions()[0][1], 8.303);
        assert_approx_eq!(frame.positions()[0][2], 11.737);

        assert_approx_eq!(frame.positions()[5][0], 8.922);
        assert_approx_eq!(frame.positions()[5][1], 9.426);
        assert_approx_eq!(frame.positions()[5][2], 5.320);
    }

    #[test]
    fn bug_in_1htq() {
        // https://github.com/chemfiles/chemfiles/issues/328
        // some secondary structure residues are not in the expected order
        let path = Path::new("./src/tests-data/pdb/1htq.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();
        let topology = frame.topology();
        // The residue IDs are out of order, but still read correctly
        let first_residue = topology.residue_for_atom(2316).unwrap();
        assert_eq!(first_residue.id.unwrap(), 503);
        assert_eq!(
            first_residue
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );

        // The 'second' residue
        let second_residue = topology.residue_for_atom(2320).unwrap();
        assert_eq!(second_residue.id.unwrap(), 287);
        assert_eq!(
            second_residue
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );

        // The 'third' residue
        let third_residue = topology.residue_for_atom(2332).unwrap();
        assert_eq!(third_residue.id.unwrap(), 288);
        assert_eq!(
            third_residue
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );

        // The 'last' residue
        let final_residue = topology.residue_for_atom(2337).unwrap();
        assert_eq!(final_residue.id.unwrap(), 289);
        assert_eq!(
            final_residue
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );

        // No secondary structure after the chain
        let no_ss_residue = topology.residue_for_atom(2341).unwrap();
        assert_eq!(no_ss_residue.id.unwrap(), 290);
        assert!(no_ss_residue.get("secondary_structure").is_none());
    }

    #[test]
    fn bug_in_1avg() {
        // https://github.com/chemfiles/chemfiles/issues/342
        // some secondary structure residues are not in the expected order
        let path = Path::new("./src/tests-data/pdb/1avg.pdb");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();
        let topology = frame.topology();

        let pre_residue = topology.residue_for_atom(75).unwrap();
        assert_eq!(pre_residue.id.unwrap(), 1);
        assert_eq!(
            pre_residue.get("insertion_code").unwrap().expect_string(),
            "D"
        );
        assert!(pre_residue.get("secondary_structure").is_none());

        let first_residue = topology.residue_for_atom(79).unwrap();
        assert_eq!(first_residue.id.unwrap(), 1);
        assert_eq!(
            first_residue.get("insertion_code").unwrap().expect_string(),
            "C"
        );
        assert_eq!(
            first_residue
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );

        let second_residue = topology.residue_for_atom(88).unwrap();
        assert_eq!(second_residue.id.unwrap(), 1);
        assert_eq!(
            second_residue
                .get("insertion_code")
                .unwrap()
                .expect_string(),
            "B"
        );
        assert_eq!(
            second_residue
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );

        let third_residue = topology.residue_for_atom(93).unwrap();
        assert_eq!(third_residue.id.unwrap(), 1);
        assert_eq!(
            third_residue.get("insertion_code").unwrap().expect_string(),
            "A"
        );
        assert_eq!(
            third_residue
                .get("secondary_structure")
                .unwrap()
                .expect_string(),
            "right-handed 3-10 helix"
        );

        let fourth_residue = topology.residue_for_atom(101).unwrap();
        assert_eq!(fourth_residue.id.unwrap(), 1);
        assert!(fourth_residue.get("insertion_code").is_none());
        assert!(fourth_residue.get("secondary_structure").is_none());
    }

    #[test]
    fn file_by_ase() {
        let path = Path::new("./src/tests-data/pdb/ase.pdb");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 156);
    }

    // TODO: fix this test - requires implementing compressed reading
    // #[test]
    // fn test_left_handed_helix() {
    //     let path = Path::new("./src/tests-data/pdb/insertion-code.pdb");
    //     let mut trajectory = Trajectory::open(&path).unwrap();
    //     let frame = trajectory.read().unwrap().unwrap();

    //     let topology = frame.topology();

    //     assert_eq!(
    //         topology
    //             .residues[226]
    //             .get("secondary_structure")
    //             .unwrap()
    //             .expect_string(),
    //         "left-handed alpha helix"
    //     );

    //     assert_eq!(
    //         topology
    //             .residues[138]
    //             .get("secondary_structure")
    //             .unwrap()
    //             .expect_string(),
    //         "right-handed alpha helix"
    //     );
    // }

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

    #[test]
    fn write_file() {
        const EXPECTED: &str = r#"MODEL    1
CRYST1   22.000   22.000   22.000  90.00  90.00  90.00 P 1           1
HETATM    1 A   A        1       1.000   2.000   3.000  1.00  0.00           A
HETATM    2 B   B        2       1.000   2.000   3.000  1.00  0.00           B
HETATM    3 C            3       1.000   2.000   3.000  1.00  0.00           C
HETATM    4 D            4       1.000   2.000   3.000  1.00  0.00           D
CONECT    1    2
CONECT    2    1
ENDMDL
MODEL    2
CRYST1   22.000   22.000   22.000  90.00  90.00  90.00 P 1           1
HETATM    1 A   A        4       1.000   2.000   3.000  1.00  0.00           A
ATOM      2 B   Bfoo A   3       1.000   2.000   3.000  1.00  0.00           B
ATOM      3 C    foo A   3       1.000   2.000   3.000  1.00  0.00           C
TER       4      foo A   3
HETATM    5 D    bar C    B      1.000   2.000   3.000  1.00  0.00      SEGM D
HETATM    6 E            5       4.000   5.000   6.000  1.00  0.00           E
HETATM    7 F    baz    -2       4.000   5.000   6.000  1.00  0.00           F
HETATM    8 G            6       4.000   5.000   6.000  1.00  0.00           G
CONECT    1    2    8
CONECT    2    1    8
CONECT    3    8
CONECT    5    8
CONECT    6    7    8
CONECT    7    6    8
CONECT    8    1    2    3    5
CONECT    8    6    7
ENDMDL
END
"#;
        let named_tmpfile = Builder::new()
            .prefix("test-pdb")
            .suffix(".pdb")
            .tempfile()
            .unwrap();

        // Write the expected content into the temp file
        let mut trajectory = Trajectory::create(named_tmpfile.path()).unwrap();

        let mut frame = Frame::new();
        frame.unit_cell = UnitCell::new_from_lengths([22.0, 22.0, 22.0]);
        let atom = Atom::new("A".into());
        frame.add_atom(atom, [1.0, 2.0, 3.0]);
        let atom = Atom::new("B".into());
        frame.add_atom(atom, [1.0, 2.0, 3.0]);
        let atom = Atom::new("C".into());
        frame.add_atom(atom, [1.0, 2.0, 3.0]);
        let atom = Atom::new("D".into());
        frame.add_atom(atom, [1.0, 2.0, 3.0]);
        frame.add_bond(0, 1, BondOrder::Unknown).unwrap();
        frame[0]
            .properties
            .insert("altloc".into(), Property::String("A".into()));
        frame[1]
            .properties
            .insert("altloc".into(), Property::String("BB".into()));

        trajectory.write(&frame).unwrap();

        let atom = Atom::new("E".into());
        frame.add_atom(atom, [4.0, 5.0, 6.0]);
        let atom = Atom::new("F".into());
        frame.add_atom(atom, [4.0, 5.0, 6.0]);
        let atom = Atom::new("G".into());
        frame.add_atom(atom, [4.0, 5.0, 6.0]);

        frame.add_bond(4, 5, BondOrder::Unknown).unwrap();
        frame.add_bond(0, 6, BondOrder::Unknown).unwrap();
        frame.add_bond(1, 6, BondOrder::Unknown).unwrap();
        frame.add_bond(1, 2, BondOrder::Unknown).unwrap(); // This bond will not be printed
        frame.add_bond(2, 6, BondOrder::Unknown).unwrap();
        frame.add_bond(3, 6, BondOrder::Unknown).unwrap();
        frame.add_bond(4, 6, BondOrder::Unknown).unwrap();
        frame.add_bond(5, 6, BondOrder::Unknown).unwrap();

        let mut residue = Residue::new("foo".into(), 3);
        residue.add_atom(1);
        residue.add_atom(2);
        residue
            .properties
            .insert("chainid".into(), Property::String("A".into()));
        residue
            .properties
            .insert("is_standard_pdb".into(), Property::Bool(true));
        residue.properties.insert(
            "composition_type".into(),
            Property::String("L-PEPTIDE LINKING".into()),
        );
        frame.add_residue(residue).unwrap();

        residue = Residue::new_from_name("barbar".into()); // Name will be truncated in output
        residue.add_atom(3);
        residue
            .properties
            .insert("chainid".into(), Property::String("CB".into()));
        residue
            .properties
            .insert("insertion_code".into(), Property::String("BB".into()));
        residue
            .properties
            .insert("segname".into(), Property::String("SEGMENT".into()));
        frame.add_residue(residue).unwrap();

        residue = Residue::new("baz".into(), -2);
        residue.add_atom(5);
        frame.add_residue(residue).unwrap();

        trajectory.write(&frame).unwrap();
        trajectory.finish().unwrap();

        let mut read_trajectory = Trajectory::open(named_tmpfile.path()).unwrap();
        assert_eq!(read_trajectory.size, 2);
        let frame1 = read_trajectory.read().unwrap().unwrap();
        assert_eq!(frame1.size(), 4);
        assert_eq!(
            frame1[0].properties.get("altloc").unwrap().expect_string(),
            "A"
        );
        assert_eq!(
            frame1[1].properties.get("altloc").unwrap().expect_string(),
            "B"
        );
        assert_eq!(read_trajectory.read().unwrap().unwrap().size(), 7);

        let contents = std::fs::read_to_string(named_tmpfile.path()).unwrap();
        assert_eq!(contents, EXPECTED);
    }
}
