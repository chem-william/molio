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
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Seek},
};
use yowl::feature::BondKind;
use yowl::feature::{AtomKind, Symbol};
use yowl::graph::Builder;
use yowl::read::read;
use yowl::read::Trace;

/// Currently, we're not handling 'CurlySMILES'
#[derive(Default)]
pub struct SMIFormat {
    /// Residue information in the current step
    pub residues: Vec<Residue>,

    first_atom: bool,
    previous_atom: usize,
    current_atom: usize,
    current_bond_order: BondOrder,
}

impl From<BondKind> for BondOrder {
    fn from(value: BondKind) -> Self {
        match value {
            BondKind::Elided => BondOrder::Single,
            BondKind::Single => BondOrder::Single,
            BondKind::Double => BondOrder::Double,
            BondKind::Triple => BondOrder::Triple,
            BondKind::Quadruple => BondOrder::Quintuplet,
            BondKind::Aromatic => BondOrder::Aromatic,
            BondKind::Up => BondOrder::Up,
            BondKind::Down => BondOrder::Down,
        }
    }
}

impl SMIFormat {
    fn add_atom<'a>(&'a mut self, topology: &'a mut Topology, atom_name: &'a str) -> &'a mut Atom {
        topology.add_atom(Atom::new(atom_name.to_string()));

        if !self.first_atom {
            self.current_atom += 1;
        }

        self.first_atom = false;
        self.previous_atom = self.current_atom;
        self.current_bond_order = BondOrder::Single;
        self.residues
            .last_mut()
            .expect("at least one residue")
            .add_atom(topology.size() - 1);
        let new_atom_idx = topology.size() - 1;
        let new_atom = topology
            .atoms
            .get_mut(new_atom_idx)
            .expect("we just added the atom");

        new_atom
    }
}

// TODO: support dative bonds (<- and ->)
// otherwise, we can't parse something like
// N->Co(<-N)(<-N)(<-N)
// CCl.[O-]>C(Cl)Cl>CO.[Cl-]
impl FileFormat for SMIFormat {
    fn read_next(&mut self, reader: &mut BufReader<File>) -> Result<Frame, CError> {
        self.residues.clear();

        let mut frame = Frame::new();
        let mut topology = Topology::default();
        // TODO: this only needs to be mutable if we handle "<" and ">"
        let mut groupid = 1;
        self.current_atom = 0;
        self.previous_atom = 0;
        self.first_atom = true;
        self.current_bond_order = BondOrder::Single;
        self.residues.push(Residue {
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

        let mut builder = Builder::default();
        let mut trace = Trace::default();
        read(smiles, &mut builder, Some(&mut trace)).expect("Failed to parse SMILES");
        let built_smiles = builder.build();

        if built_smiles.is_err() {
            let e = line.trim();
            warn!("could not parse '{e}'. skipping for now");
            eprintln!("could not parse '{e}'. skipping for now");
            return Ok(Frame::new());
        }
        let mut parsed_smiles = built_smiles.unwrap();

        topology.atoms.reserve(parsed_smiles.len());
        for (i, yowl_atom) in parsed_smiles.iter_mut().enumerate() {
            // we need to invert the configuration on each atom to match the SMILES
            // instead of how the graph was built.
            // The internal viewpoint of the graph in `yowl` is child -> parent which
            // is the mirror of the SMILES' parent -> child
            yowl_atom.kind.invert_configuration();

            if !self.first_atom {
                // XXX: `yowl` does not emit information about '.' (that is, splits), we
                // manually parse part of the SMILES to add the residues
                let prev_range_end = trace.atom(i - 1).unwrap().end;
                let curr_range_start = trace.atom(i).unwrap().start;
                if curr_range_start - prev_range_end > 1 {
                    let check_hole = &smiles[prev_range_end + 1..curr_range_start];
                    if check_hole == "." {
                        self.residues.push(Residue {
                            name: format!("group {groupid}"),
                            ..Default::default()
                        });
                    }
                }
            }
            match yowl_atom.kind {
                AtomKind::Symbol(Symbol::Star) => {
                    self.add_atom(&mut topology, "*");
                }
                AtomKind::Symbol(Symbol::Aliphatic(element)) => {
                    self.add_atom(&mut topology, element.symbol());
                }
                AtomKind::Symbol(Symbol::Aromatic(element)) => {
                    let new_atom = self.add_atom(&mut topology, element.symbol());
                    new_atom
                        .properties
                        .insert("is_aromatic".to_string(), Property::Bool(true));
                }
                AtomKind::Bracket {
                    isotope,
                    symbol,
                    configuration,
                    hcount,
                    charge,
                    map,
                } => {
                    let new_atom = match symbol {
                        Symbol::Star => self.add_atom(&mut topology, "*"),
                        Symbol::Aliphatic(element) => {
                            self.add_atom(&mut topology, element.symbol())
                        }
                        Symbol::Aromatic(element) => {
                            let new_atom = self.add_atom(&mut topology, element.symbol());
                            new_atom
                                .properties
                                .insert("is_aromatic".to_string(), Property::Bool(true));
                            new_atom
                        }
                    };

                    if let Some(isotope) = isotope {
                        new_atom.mass = f64::from(isotope.mass_number());
                    }

                    if let Some(hcount) = hcount {
                        new_atom.properties.insert(
                            "hydrogen_count".to_string(),
                            Property::Double(f64::from(u8::from(&hcount))),
                        );
                    }

                    if let Some(charge) = charge {
                        new_atom.charge = f64::from(i8::from(charge));
                    }

                    if let Some(map) = map {
                        new_atom.properties.insert(
                            "smiles_class".to_string(),
                            Property::String(map.to_string()),
                        );
                    }

                    if let Some(configuration) = configuration {
                        new_atom.properties.insert(
                            "chirality".to_string(),
                            Property::String(configuration.to_string()),
                        );
                    }
                }
            }

            // For each neighbor index < i, add a bond.
            for j in &yowl_atom.bonds {
                if j.tid < i {
                    topology.add_bond(j.tid, i, j.kind.into())?;
                }
            }
        }

        for residue in &mut self.residues {
            topology
                .add_residue(residue.clone())
                .expect("able to add residues to topology");
        }

        if smiles.split_whitespace().count() > 1 {
            let name = smiles
                .split_whitespace()
                .skip(1)
                .collect::<Vec<_>>()
                .join(" ");
            frame
                .properties
                .insert("name".to_string(), Property::String(name.to_string()));
        }

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

    fn finalize(
        &self,
        writer: &mut std::io::BufWriter<std::fs::File>,
    ) -> Result<(), crate::error::CError> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::{bond::BondOrder, frame::Frame, property::Property, trajectory::Trajectory};
    use std::{hint::black_box, path::Path};

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

    // we've removed the following test case as `yowl` does not support single-quotation
    // marks in the SMILES string
    // ['Db']['Sg']['Bh']['Hs']['Mt']['Ds']['Rg']['Cn']['Nh']['Fl']['Mc']['Lv']['Ts']['Og']
    #[test]
    fn read_entire_file() {
        let path = Path::new("./src/tests-data/smi/rdkit_problems.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 69);

        let mut frame = Frame::new();
        while let Some(next_frame) = trajectory.read().unwrap() {
            frame = next_frame;
        }
        assert_eq!(frame.size(), 14);
        assert_eq!(frame[0].symbol, "Db");
        assert_eq!(frame[13].symbol, "Og");
    }

    #[test]
    fn check_parsing_results() {
        let path = Path::new("./src/tests-data/smi/details.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 1);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 5);
        assert_eq!(frame[0].charge, 0.0);
        assert_eq!(frame[0].symbol, "O");
        assert_eq!(frame[4].charge, -1.0);
        assert_eq!(frame[4].symbol, "O");
    }

    #[test]
    fn ugly_smiles_strings() {
        let path = Path::new("./src/tests-data/smi/ugly.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 3);

        // C1(CC1CC1CC1)
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 7);
        let bonds = frame.topology().bonds();
        assert_eq!(bonds.len(), 8);
        assert!((bonds[0][0] == 0 && bonds[0][1] == 1));
        assert!((bonds[1][0] == 0 && bonds[1][1] == 2));
        assert!((bonds[2][0] == 1 && bonds[2][1] == 2));
        assert!((bonds[3][0] == 2 && bonds[3][1] == 3));
        assert!((bonds[4][0] == 3 && bonds[4][1] == 4));
        assert!((bonds[5][0] == 4 && bonds[5][1] == 5));
        assert!((bonds[6][0] == 4 && bonds[6][1] == 6));
        assert!((bonds[7][0] == 5 && bonds[7][1] == 6));

        // C1.C1CC1CC1
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 6);
        let bonds = frame.topology().bonds();
        assert_eq!(bonds.len(), 6);
        assert!((bonds[0][0] == 0 && bonds[0][1] == 1));
        assert!((bonds[1][0] == 1 && bonds[1][1] == 2));
        assert!((bonds[2][0] == 2 && bonds[2][1] == 3));
        assert!((bonds[3][0] == 3 && bonds[3][1] == 4));
        assert!((bonds[4][0] == 3 && bonds[4][1] == 5));
        assert!((bonds[5][0] == 4 && bonds[5][1] == 5));
        let topology = frame.topology();
        assert_eq!(topology.residues.len(), 2);
        assert!(topology.are_linked(&topology.residues[0], &topology.residues[1]));

        // C1CC11CC1
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 5);
        let bonds = frame.topology().bonds();
        assert_eq!(bonds.len(), 6);

        assert!((bonds[0][0] == 0 && bonds[0][1] == 1));
        assert!((bonds[1][0] == 0 && bonds[1][1] == 2));
        assert!((bonds[2][0] == 1 && bonds[2][1] == 2));
        assert!((bonds[3][0] == 2 && bonds[3][1] == 3));
        assert!((bonds[4][0] == 2 && bonds[4][1] == 4));
        assert!((bonds[5][0] == 3 && bonds[5][1] == 4));
    }

    #[test]
    fn rdkit_problems() {
        let path = Path::new("./src/tests-data/smi/rdkit_problems.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 69);

        // C1CC2C1CC2
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 6);
        let bonds = frame.topology().bonds();
        assert_eq!(bonds.len(), 7);
        assert!((bonds[0][0] == 0 && bonds[0][1] == 1));
        assert!((bonds[1][0] == 0 && bonds[1][1] == 3));
        assert!((bonds[2][0] == 1 && bonds[2][1] == 2));
        assert!((bonds[3][0] == 2 && bonds[3][1] == 3));
        assert!((bonds[4][0] == 2 && bonds[4][1] == 5));
        assert!((bonds[5][0] == 3 && bonds[5][1] == 4));
        assert!((bonds[6][0] == 4 && bonds[6][1] == 5));

        // [CH2+]C[CH+2]
        let frame = trajectory.read_at(6).unwrap().unwrap();
        assert_eq!(
            *frame[0].properties.get("hydrogen_count").unwrap(),
            Property::Double(2.0)
        );
        assert_eq!(frame[0].charge, 1.0);
        assert_eq!(
            *frame[2].properties.get("hydrogen_count").unwrap(),
            Property::Double(1.0)
        );
        assert_eq!(frame[2].charge, 2.0);

        // C1CC=1
        let frame = trajectory.read_at(8).unwrap().unwrap();
        let bond_orders = frame.topology().bond_orders();
        assert_eq!(bond_orders[0], BondOrder::Single);
        assert_eq!(bond_orders[1], BondOrder::Double);

        // C=1CC1
        let frame = trajectory.read_at(9).unwrap().unwrap();
        let bond_orders = frame.topology().bond_orders();
        assert_eq!(bond_orders[0], BondOrder::Single);
        assert_eq!(bond_orders[1], BondOrder::Double);
    }

    #[test]
    fn chirality() {
        let path = Path::new("./src/tests-data/smi/chiral.smi");
        let mut trajectory = Trajectory::open(path).unwrap();

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(
            *frame[1].properties.get("chirality").unwrap(),
            Property::String("@TB1".to_string())
        );
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(
            *frame[1].properties.get("chirality").unwrap(),
            Property::String("@TB15".to_string())
        );
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(
            *frame[1].properties.get("chirality").unwrap(),
            Property::String("@@".to_string())
        );
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(
            *frame[1].properties.get("chirality").unwrap(),
            Property::String("@OH15".to_string())
        );
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(
            *frame[1].properties.get("chirality").unwrap(),
            Property::String("@@".to_string())
        );
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(
            *frame[1].properties.get("chirality").unwrap(),
            Property::String("@".to_string())
        );
    }

    #[test]
    fn other_tests() {
        let path = Path::new("./src/tests-data/smi/test.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        trajectory.read().unwrap().unwrap();
        let mut frame = trajectory.read().unwrap().unwrap();

        assert_eq!(
            *frame[0].properties.get("is_aromatic").unwrap(),
            Property::Bool(true)
        );
        assert_eq!(
            *frame.properties.get("name").unwrap(),
            Property::String("Benzene".to_string())
        );

        while let Some(next_frame) = trajectory.read().unwrap() {
            frame = next_frame;
        }

        black_box(frame);
    }

    #[test]
    fn issue_303() {
        let path = Path::new("./src/tests-data/smi/issue_303.smi");
        let mut trajectory = Trajectory::open(path).unwrap();

        // In issue 303 (from the original c++ chemfiles lib), this failed due to the "%11" marker
        trajectory.read().unwrap().unwrap();

        // No explicit hydrogens, so the size should be 26 atoms
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 26);

        // Converting the original SDF file using MarvinSketch preverses the explicit hydrogens
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 30);

        // For the next test, too many bonds were parsed
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.topology().bonds().len(), 34);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.topology().bonds().len(), 182);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.topology().bonds().len(), 171);
    }

    #[test]
    #[should_panic(expected = "Character(2)")]
    fn bad_element() {
        let path = Path::new("./src/tests-data/smi/bad_element.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        trajectory.read().unwrap().unwrap();
    }

    #[test]
    #[should_panic(expected = "Character(2)")]
    fn bad_paren() {
        let path = Path::new("./src/tests-data/smi/bad_paren.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        trajectory.read().unwrap().unwrap();
    }

    #[test]
    #[should_panic(expected = "EndOfLine")]
    fn bad_percentage() {
        let path = Path::new("./src/tests-data/smi/bad_percentage_sign.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        trajectory.read().unwrap().unwrap();
    }

    #[test]
    #[should_panic(expected = "Character(2)")]
    fn bad_ring() {
        let path = Path::new("./src/tests-data/smi/bad_ring.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        trajectory.read().unwrap().unwrap();
    }

    #[test]
    #[should_panic(expected = "Character(3)")]
    fn bad_symbol() {
        let path = Path::new("./src/tests-data/smi/bad_symbol.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        trajectory.read().unwrap().unwrap();
    }

    #[test]
    #[should_panic(expected = "Character(2)")]
    fn misplaced_property() {
        let path = Path::new("./src/tests-data/smi/misplaced_property.smi");
        let mut trajectory = Trajectory::open(path).unwrap();
        trajectory.read().unwrap().unwrap();
    }
}
