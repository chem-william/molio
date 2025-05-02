// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::atom::Atom;
use crate::bond::BondOrder;
use crate::property::Property;
use crate::residue::Residue;
use crate::topology::Topology;
use crate::{error::CError, format::FileFormat, frame::Frame};
use log::warn;
use purr::feature::Aromatic;
use purr::graph::Builder;
use purr::read::read;
use std::cell::RefCell;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Seek},
};

/// Currently, we're not handling 'CurlySMILES'
#[derive(Default)]
pub struct SMIFormat {
    /// Residue information in the current step
    pub residues: RefCell<Vec<Residue>>,
}

// TODO: support dative bonds (<- and ->)
// otherwise, we can't parse something like
// N->Co(<-N)(<-N)(<-N)
// CCl.[O-]>C(Cl)Cl>CO.[Cl-]
// probably requires upstream changes to Purr/Balsa unless we roll our own or
// find another lib for parsing SMILES
impl FileFormat for SMIFormat {
    fn read_next(&mut self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        self.residues.borrow_mut().clear();

        let mut frame = Frame::new();
        let mut topology = Topology::default();
        let mut groupid = 0;
        self.residues.borrow_mut().push(Residue {
            name: format!("group {groupid}"),
            ..Default::default()
        });

        let mut line = String::new();
        reader.read_line(&mut line)?;
        let mut smiles = line.trim();
        while smiles.is_empty() {
            line.clear();
            reader.read_line(&mut line)?;
            smiles = line.trim();
        }

        let mut builder = Builder::new();
        let _ = read(smiles, &mut builder, None);

        dbg!(&line);
        let built_smiles = builder.build();
        if built_smiles.is_err() {
            let e = line.trim();
            warn!("could not parse '{e}'. skipping for now");
            eprintln!("could not parse '{e}'. skipping for now");
            return Ok(Frame::new());
        }
        let parsed_smiles = built_smiles.unwrap();

        topology.atoms.reserve(parsed_smiles.len());
        for (i, purr) in parsed_smiles.iter().enumerate() {
            // Add the atom itself.
            topology.add_atom(Atom::new(purr.kind.to_string()));
            match purr.kind {
                purr::feature::AtomKind::Aromatic(_) => topology[0]
                    .properties
                    .insert("is_aromatic".to_string(), Property::Bool(true)),
                _ => None,
            };

            // Get a mutable ref to the newly added atom:
            // let atom = topology.atoms.last_mut();

            // For each neighbor index < i, add a bond.
            for j in &purr.bonds {
                if j.tid < i {
                    // TODO: convert the bond order from Purr
                    topology.add_bond(j.tid, i, BondOrder::Single)?;
                }
            }

            // Record this atom in the current residue.
            let topo_size = topology.size();
            self.residues
                .borrow_mut()
                .last_mut()
                .expect("no residue")
                .add_atom(topo_size - 1);
        }

        // dbg!(&line);
        frame.resize(topology.size())?;
        frame.set_topology(topology)?;
        Ok(frame)
    }

    fn read(&mut self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        Ok(Some(self.read_next(reader)?))
    }

    fn write_next(
        &mut self,
        writer: &mut BufWriter<File>,
        frame: &crate::frame::Frame,
    ) -> Result<(), crate::error::CError> {
        todo!()
    }

    fn forward(&self, reader: &mut BufReader<File>) -> Result<Option<u64>, CError> {
        let pos: u64;

        let mut line = String::new();
        loop {
            match reader.read_line(&mut line)? {
                0 => return Ok(None),                    // EOF
                _ if line.trim().is_empty() => continue, // skip blank
                _ => {
                    pos = reader.stream_position()?; // position after this line
                    break;
                }
            }
        }

        Ok(Some(pos))
    }
    // optional<uint64_t> SMIFormat::forward() {
    //     auto position = file_.tellpos();

    //     auto line = file_.readline();
    //     while (trim(line).empty()) {
    //         if (file_.eof()) {
    //             return nullopt;
    //         }
    //         line = file_.readline();
    //     }

    //     return position;
    // }

    fn finalize(
        &self,
        writer: &mut std::io::BufWriter<std::fs::File>,
    ) -> Result<(), crate::error::CError> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        bond::BondOrder,
        frame::Frame,
        topology,
        trajectory::{self, Trajectory},
    };
    use std::path::Path;

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/smi/test.smi");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 8);
    }

    #[test]
    fn check_nsteps_with_newlines() {
        let path = Path::new("./src/tests-data/smi/spaces.smi");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 8);
    }

    #[test]
    fn read_next_frame() {
        let path = Path::new("./src/tests-data/smi/test.smi");
        let mut trajectory = Trajectory::open(path).unwrap();

        // Check to make sure things aren't exploding
        // reading: C1CC2C1CC2
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 6);
        assert_eq!(frame.topology().bonds().len(), 7);

        // reading: c1ccccc1	Benzene
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 6);
        assert_eq!(frame.topology().bonds().len(), 6);

        // reading: C(Cl)(Cl)(Cl)
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 4);
        let topology = frame.topology();
        let bonds = topology.bonds();
        assert_eq!(bonds.len(), 3);
        assert!(bonds[0][0] == 0 && bonds[0][1] == 1);
        assert!(bonds[1][0] == 0 && bonds[1][1] == 2);
        assert!(bonds[2][0] == 0 && bonds[2][1] == 3);
        assert_eq!(frame[0].symbol, "C");
        assert_eq!(frame[1].symbol, "Cl");
        assert_eq!(frame[2].symbol, "Cl");
        assert_eq!(frame[3].symbol, "Cl");
    }

    #[test]
    fn read_specific_step() {
        let path = Path::new("./src/tests-data/smi/test.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read_at(1).unwrap().unwrap();
        assert_eq!(frame.size(), 6);
        let topology = frame.topology();
        assert_eq!(topology.bonds().len(), 6);

        // TODO: requires correct parsing of dative bonds (<- and ->)
        // let frame = trajectory.read_at(7).unwrap().unwrap();
        // assert_eq!(frame.size(), 9);
        // let topology = frame.topology();
        // assert_eq!(topology.bonds().len(), 6);

        let frame = trajectory.read_at(5).unwrap().unwrap();
        assert_eq!(frame.size(), 6);
    }

    #[test]
    fn read_specific_step_with_whitespace() {
        let path = Path::new("./src/tests-data/smi/spaces.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read_at(1).unwrap().unwrap();
        assert_eq!(frame.size(), 6);
        let topology = frame.topology();
        assert_eq!(topology.bonds().len(), 6);

        // TODO: requires correct parsing of dative bonds (<- and ->)
        // let frame = trajectory.read_at(7).unwrap().unwrap();
        // assert_eq!(frame.size(), 9);
        // let topology = frame.topology();
        // assert_eq!(topology.bonds().len(), 6);

        let frame = trajectory.read_at(5).unwrap().unwrap();
        assert_eq!(frame.size(), 6);
        let topology = frame.topology();
        assert_eq!(topology.bonds().len(), 6);

        // Check that calling trajectory.read() repeatedly is the same as frame.read_at()
        let path = Path::new("./src/tests-data/smi/spaces.smi");
        let mut trajectory = Trajectory::open(path).unwrap();

        trajectory.read().unwrap().unwrap();
        trajectory.read().unwrap().unwrap();
        trajectory.read().unwrap().unwrap();
        trajectory.read().unwrap().unwrap();
        trajectory.read().unwrap().unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        assert_eq!(frame.size(), 6);
        let topology = frame.topology();
        assert_eq!(topology.bonds().len(), 6);
    }

    #[test]
    fn read_entire_file() {
        let path = Path::new("./src/tests-data/smi/rdkit_problems.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 70);

        let mut frame = Frame::new();
        while let Some(next_frame) = trajectory.read().unwrap() {
            frame = next_frame;
        }
        // TODO: until we have good support for SMILES, we rewind to the latest SMILES'
        // in the test file
        let frame = trajectory.read_at(67).unwrap().unwrap();
        assert_eq!(frame.size(), 14);
        assert_eq!(frame[0].symbol, "[Db]");
        assert_eq!(frame[13].symbol, "[Og]");
    }
}
