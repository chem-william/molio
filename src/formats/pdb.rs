use crate::error::CError;
use crate::format::FileFormat;
use crate::frame::Frame;
use crate::property::Properties;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Seek};
use std::path::Path;

struct FullResidueId {
    /// Chain identifier
    chain: char,
    /// Residue id
    resid: i64,
    /// Residue name
    resname: String,
    /// Insertion code of the residue
    insertion_code: char,
}
struct Residue {
    name: String,
    id: Option<usize>,
    atoms: BTreeSet<usize>,
    properties: Properties,
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
    let rec = &line[..6];

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
    pub residues: BTreeMap<FullResidueId, Residue>,
}

impl PDBFormat {
    // fn parse_atom_line(&self, line: &str) -> Result<Atom, CError> {
    //     if line.len() < 54 {
    //         return Err(CError::InvalidRecord {
    //             record_type: "ATOM".to_string(),
    //             reason: "line too short".to_string(),
    //         });
    //     }

    //     let serial = line[6..11].trim().parse::<usize>().map_err(|e| CError::ParseError {
    //         record_type: "ATOM".to_string(),
    //         field: "serial number".to_string(),
    //         error: e.to_string(),
    //     })?;

    //     let name = line[12..16].trim().to_string();
    //     let alt_loc = line[16..17].trim().to_string();
    //     let res_name = line[17..20].trim().to_string();
    //     let chain_id = line[21..22].trim().to_string();
    //     let res_seq = line[22..26].trim().parse::<i32>().map_err(|e| CError::ParseError {
    //         record_type: "ATOM".to_string(),
    //         field: "residue sequence number".to_string(),
    //         error: e.to_string(),
    //     })?;
    //     let i_code = line[26..27].trim().to_string();

    //     let x = line[30..38].trim().parse::<f64>().map_err(|e| CError::ParseError {
    //         record_type: "ATOM".to_string(),
    //         field: "x coordinate".to_string(),
    //         error: e.to_string(),
    //     })?;
    //     let y = line[38..46].trim().parse::<f64>().map_err(|e| CError::ParseError {
    //         record_type: "ATOM".to_string(),
    //         field: "y coordinate".to_string(),
    //         error: e.to_string(),
    //     })?;
    //     let z = line[46..54].trim().parse::<f64>().map_err(|e| CError::ParseError {
    //         record_type: "ATOM".to_string(),
    //         field: "z coordinate".to_string(),
    //         error: e.to_string(),
    //     })?;

    //     let element = if line.len() >= 78 {
    //         line[76..78].trim().to_string()
    //     } else {
    //         "".to_string()
    //     };

    //     Ok(Atom {
    //         serial,
    //         name,
    //         alt_loc,
    //         res_name,
    //         chain_id,
    //         res_seq,
    //         i_code,
    //         x,
    //         y,
    //         z,
    //         element,
    //     })
    // }

    // fn parse_cryst1_line(&self, line: &str) -> Result<UnitCell, CError> {
    //     if line.len() < 54 {
    //         return Err(CError::InvalidRecord {
    //             record_type: "CRYST1".to_string(),
    //             reason: "line too short".to_string(),
    //         });
    //     }

    //     let a = line[6..15].trim().parse::<f64>().map_err(|e| CError::ParseError {
    //         record_type: "CRYST1".to_string(),
    //         field: "a parameter".to_string(),
    //         error: e.to_string(),
    //     })?;
    //     let b = line[15..24].trim().parse::<f64>().map_err(|e| CError::ParseError {
    //         record_type: "CRYST1".to_string(),
    //         field: "b parameter".to_string(),
    //         error: e.to_string(),
    //     })?;
    //     let c = line[24..33].trim().parse::<f64>().map_err(|e| CError::ParseError {
    //         record_type: "CRYST1".to_string(),
    //         field: "c parameter".to_string(),
    //         error: e.to_string(),
    //     })?;
    //     let alpha = line[33..40].trim().parse::<f64>().map_err(|e| CError::ParseError {
    //         record_type: "CRYST1".to_string(),
    //         field: "alpha angle".to_string(),
    //         error: e.to_string(),
    //     })?;
    //     let beta = line[40..47].trim().parse::<f64>().map_err(|e| CError::ParseError {
    //         record_type: "CRYST1".to_string(),
    //         field: "beta angle".to_string(),
    //         error: e.to_string(),
    //     })?;
    //     let gamma = line[47..54].trim().parse::<f64>().map_err(|e| CError::ParseError {
    //         record_type: "CRYST1".to_string(),
    //         field: "gamma angle".to_string(),
    //         error: e.to_string(),
    //     })?;

    //     Ok(UnitCell {
    //         a,
    //         b,
    //         c,
    //         alpha,
    //         beta,
    //         gamma,
    //     })
    // }

    const END_RECORD: &str = "END";
    const ENDMDL_RECORD: &str = "ENDMDL";
}

impl FileFormat for PDBFormat {
    fn read_next(&self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        self.residues.clear();
        let mut line = String::new();

        while reader.read_line(line)? > 0 {}

        // let mut frame = Frame::new();
        // let mut line = String::new();

        // while reader.read_line(&mut line)? > 0 {
        //     let trimmed = line.trim();
        //     if trimmed.is_empty() {
        //         continue;
        //     }

        //     if trimmed.starts_with("CRYST1") {
        //         let unit_cell = Self::parse_cryst1_line(self, trimmed)?;
        //         // frame.unit_cell = unit_cell;
        //     } else if trimmed.starts_with("ATOM") || trimmed.starts_with("HETATM") {
        //         let atom = Self::parse_atom_line(self, trimmed)?;
        //         frame.add_atom(atom);
        //     } else if trimmed.starts_with("END") {
        //         break;
        //     }

        //     line.clear();
        // }

        // Ok(frame)
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
        todo!();
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::{frame::Frame, trajectory::Trajectory, unit_cell::UnitCell};
    use assert_approx_eq::assert_approx_eq;

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
    fn hybrid_decode() {
        assert_eq!(decode_hybrid36(4, "    ").unwrap(), 0);
        assert_eq!(decode_hybrid36(4, "  -0").unwrap(), 0);
    }
}
