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
}
