// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::atom::Atom;
use crate::bond::BondOrder;
use crate::property::{Property, PropertyKind};
use crate::residue::Residue;
use crate::topology::Topology;
use crate::{error::CError, format::FileFormat, frame::Frame};
use log::warn;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Seek, Write},
};
use yowl::feature::BondKind;
use yowl::feature::{AtomKind, Symbol};
use yowl::graph::Builder;
use yowl::read::read;
use yowl::read::Trace;

impl fmt::Display for BondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            BondOrder::Unknown => "~",
            BondOrder::Single => "",
            BondOrder::Double => "=",
            BondOrder::Triple => "#",
            BondOrder::Quadruple => "$",
            BondOrder::Quintuplet => "",
            BondOrder::Down => "\\",
            BondOrder::Up => "/",
            BondOrder::DativeR => "",
            BondOrder::DativeL => "",
            BondOrder::Amide => "",
            BondOrder::Aromatic => ":",
        };
        write!(f, "{s}")
    }
}

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
    fn is_aliphatic_organic(s: &str) -> bool {
        matches!(
            s,
            "B" | "C" | "N" | "O" | "S" | "P" | "F" | "Cl" | "Br" | "I" | "H"
        )
    }

    fn is_chirality_tag(s: &str) -> bool {
        matches!(s, "TH" | "SP" | "TB" | "OH" | "AL")
    }

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

    // Helper to check if all atoms have been hit/processed.
    fn all_hit(hit: &[bool]) -> bool {
        hit.iter().all(|&b| b)
    }

    fn find_rings(adj_list: &[Vec<usize>]) -> HashMap<usize, usize> {
        let mut ring_atoms = HashMap::new();

        let n_atoms = adj_list.len();

        let mut hit_atoms = vec![false; n_atoms];
        let mut ring_bonds: HashSet<(usize, usize)> = HashSet::new();

        while !SMIFormat::all_hit(&hit_atoms) {
            // Find index of first `false` in `hit_atoms`.
            let current_atom = hit_atoms
                .iter()
                .position(|&b| !b)
                .expect("we just checked that not all are hit");
            // Mark it as processed.
            hit_atoms[current_atom] = true;

            // If this atom has no neighbors, it cannot be part of any ring → skip.
            if adj_list[current_atom].is_empty() {
                continue;
            }

            // We'll perform a DFS starting from `current_atom`, treating `(current, previous)` as the state.
            // Initialize stack with (start, previous = start).
            let mut stack: Vec<(usize, usize)> = Vec::new();
            stack.push((current_atom, current_atom));

            // As long as there are atoms to process in this connected component...
            while let Some((cur, prev)) = stack.pop() {
                // If you want to mirror naming: `current_atom` ← `cur`, `previous_atom` ← `prev`.
                let current = cur;
                let previous = prev;

                // For each neighbor of `current`, in reverse order to match C++'s use of a reverse‐iterator:
                //   C++ did `for (auto neighbor_iter = current_atom_bonds.rbegin(); ...)`.
                //
                // Rust's `Vec` has `iter().rev()` to iterate in reverse.
                for &neighbor in adj_list[current].iter().rev() {
                    // Skip backtracking to the immediate parent.
                    if neighbor == previous {
                        continue;
                    }

                    // If we have already “hit” this neighbor and have **not** yet marked `current` as hit,
                    // that means we've found a ring bond coming “backward” to an already‐seen vertex.
                    //
                    // In the original C++, they checked `if (hit_atoms[neighbor] && !hit_atoms[current_atom])`
                    // to ensure they only count each ring once, in the moment they first discover that bond.
                    if hit_atoms[neighbor] && !hit_atoms[current] {
                        // Canonicalize the unordered pair as (min, max).
                        let bond = if neighbor < current {
                            (neighbor, current)
                        } else {
                            (current, neighbor)
                        };

                        // If this bond has already been recorded, skip it.
                        if ring_bonds.contains(&bond) {
                            continue;
                        }

                        // Otherwise, insert it and increment ring‐count for `neighbor`.
                        ring_bonds.insert(bond);

                        // Increment ring count for `neighbor`:
                        let count = ring_atoms.entry(neighbor).or_insert(0);
                        *count += 1;

                        // We do NOT re‐push neighbor, because that “branch” of search is done—it closed a loop.
                        continue;
                    }

                    // If the neighbor has not been hit yet, push it onto the stack to continue DFS.
                    if !hit_atoms[neighbor] {
                        stack.push((neighbor, current));
                    }
                }

                // Now that all neighbors of `current` have been considered, mark `current` as hit.
                hit_atoms[current] = true;
            }
        }

        ring_atoms
    }

    /// Returns `Some(int_part)` if `value` has no fractional part;
    /// otherwise returns `None`.
    fn as_integer(value: f64) -> Option<i32> {
        if value.fract() == 0.0 && value != 0.0 {
            // Casting is safe because we know there’s no fractional portion.
            Some(value as i32)
        } else {
            None
        }
    }

    fn write_atom_smiles(writer: &mut BufWriter<File>, atom: &Atom) -> Result<(), CError> {
        let mut needs_brackets = false;
        let mut symbol = atom.symbol.clone();

        // The mass must be an integer is the only check we need as all atoms in the
        // periodic table have non-integer masses. Therefore, if the the mass is
        // an integer, then we know the user has set an isotope.
        let mass_int = SMIFormat::as_integer(atom.mass);
        if mass_int.is_some() {
            needs_brackets = true;
        }

        // If any of these values are set / not a default value
        let chirality = atom
            .properties
            .get("chirality")
            .cloned()
            .unwrap_or(Property::String("".to_string()));
        if !chirality.expect_string().is_empty() {
            needs_brackets = true;
        }

        let smi_class = if let Some(smi_class_prop) = atom.properties.get("smiles_class") {
            if smi_class_prop.kind() == PropertyKind::Double {
                if let Some(val) = SMIFormat::as_integer(smi_class_prop.expect_double()) {
                    if val >= 0 {
                        needs_brackets = true;
                        Some(smi_class_prop)
                    } else {
                        warn!("the 'smiles_class' property must be an integer >= 0: {val}");
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };

        let explicit_h = if let Some(explicit_h_prop) = atom.properties.get("hydrogen_count") {
            if explicit_h_prop.kind() == PropertyKind::Double {
                if let Some(val) = SMIFormat::as_integer(explicit_h_prop.expect_double()) {
                    if val >= 0 {
                        needs_brackets = true;
                        Some(val)
                    } else {
                        warn!("the 'hydrogen_count' property must be an integer >= 0: {val}");
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };

        // Charge is similar to mass. It must be an integer, otherwise it is ignored
        let charge_int = if let Some(charge_int) = SMIFormat::as_integer(atom.charge) {
            // Use |= because we don't want to unset if charge_int is 0
            needs_brackets |= charge_int != 0;
            Some(charge_int)
        } else {
            None
        };

        // Before the atom is printed, we must know if we are prining a property set.
        if !SMIFormat::is_aliphatic_organic(symbol.as_str()) {
            needs_brackets = true;
        }

        let is_aromatic = atom
            .properties
            .get("is_aromatic")
            .unwrap_or(&Property::Bool(false))
            .expect_bool();
        if is_aromatic {
            // We need to print brackets for aromatic Te, Se, As, and Si
            // Note, we allow aromatic Boron to be bare (no brackets) in following
            // ChemAxon's standard.
            if symbol.len() > 1 {
                needs_brackets = true;
            }
            symbol = symbol.to_ascii_lowercase();
        }

        if needs_brackets {
            write!(writer, "[")?;
        }

        // Mass must be first, before the element is printed
        if mass_int.is_some() {
            write!(writer, "{}", mass_int.unwrap())?;
        }

        write!(writer, "{symbol}")?;

        if smi_class.is_some() {
            write!(writer, ":{}", smi_class.unwrap().expect_double())?;
        }

        let mut is_good_tag = false;
        let chirality_string = chirality.expect_string();
        match chirality_string.len() {
            0 => is_good_tag = true,
            2 => {
                if chirality == Property::String("CW".to_string()) {
                    is_good_tag = true;
                    write!(writer, "@@")?;
                }
            }
            3 => {
                if chirality == Property::String("CCW".to_string()) {
                    is_good_tag = true;
                    write!(writer, "@")?;
                }
            }
            7 => {
                if chirality_string.starts_with("CC")
                    && SMIFormat::is_chirality_tag(&chirality_string[4..6])
                    && chirality_string
                        .as_bytes()
                        .get(6)
                        .map(|b| b.is_ascii_digit())
                        .unwrap_or(false)
                {
                    is_good_tag = true;
                    write!(writer, "@{}", &chirality_string[4..])?;
                }
            }
            8 => {
                if chirality_string.starts_with("CCW")
                    && SMIFormat::is_chirality_tag(&chirality_string[4..6])
                    && chirality_string
                        .as_bytes()
                        .get(6..8)
                        .map(|b| b[0].is_ascii_digit() && b[1].is_ascii_digit())
                        .unwrap_or(false)
                {
                    is_good_tag = true;
                    write!(writer, "@{}", &chirality_string[4..])?;
                }
            }
            _ => {}
        }

        if !is_good_tag {
            warn!("invalid chirality tag '{chirality_string}'");
        }

        if let Some(val) = explicit_h {
            if val == 1 {
                write!(writer, "H")?;
            } else {
                write!(writer, "H{val}")?;
            }
        }

        if let Some(val) = charge_int {
            if val < 0 {
                write!(writer, "-")?;
            } else {
                write!(writer, "+")?;
            }

            let abs_val = val.abs();
            if abs_val != 1 && abs_val != 0 {
                write!(writer, "{val}")?;
            }
        }

        if needs_brackets {
            write!(writer, "]")?;
        }

        Ok(())
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
        let groupid = 1;
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

        let mut parts = smiles.split_whitespace();
        let smiles_part = parts.next().unwrap_or("");
        let name_part = parts.collect::<Vec<_>>().join(" ");

        if !name_part.is_empty() {
            frame
                .properties
                .insert("name".to_string(), Property::String(name_part));
        }

        let smiles = smiles_part;

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

        frame.resize(topology.size())?;
        frame.set_topology(topology)?;
        Ok(frame)
    }

    fn read(&mut self, reader: &mut BufReader<File>) -> Result<Option<Frame>, CError> {
        Ok(Some(self.read_next(reader)?))
    }

    fn write_next(&mut self, writer: &mut BufWriter<File>, frame: &Frame) -> Result<(), CError> {
        if frame.size() == 0 {
            writeln!(writer)?;
        };

        let mut adj_list: Vec<Vec<usize>> = Vec::with_capacity(frame.size());
        adj_list.extend((0..frame.size()).map(|_| Vec::new()));
        for bond in frame.topology().bonds() {
            adj_list[bond[0]].push(bond[1]);
            adj_list[bond[1]].push(bond[0]);
        }

        let ring_atoms = SMIFormat::find_rings(&adj_list);
        let mut written = vec![false; frame.size()];
        let mut branch_stack = 0;

        let mut ring_stack: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut ring_count = 0;

        let mut first_atom = true;

        while !SMIFormat::all_hit(&written) {
            if !first_atom {
                write!(writer, ".")?;
            }
            // Find index of first `false` in `hit_atoms`.
            let start_atom = written
                .iter()
                .position(|&b| !b)
                .expect("we just checked that not all are hit");
            // We have found an atom that has not yet been printed! Now we must print it out
            // along with all of its connections (the entire component).

            // A structure to store a depth first like search through the component.
            let mut atoms_to_process: Vec<(usize, usize, bool)> =
                vec![(start_atom, start_atom, false)];

            while let Some((previous_atom, current_atom, needs_branch)) = atoms_to_process.pop() {
                // Skip if we've already written this atom
                if written[current_atom] {
                    continue;
                }

                let current_atom_bonds = &adj_list[current_atom];
                written[current_atom] = true;

                // If this node needs a branching '('
                if needs_branch {
                    write!(writer, "(")?;
                    branch_stack += 1;
                }

                // Print the bond symbol unless this is the very first atom in this chain
                if current_atom != previous_atom {
                    let bo = frame.topology().bond_order(previous_atom, current_atom)?;
                    write!(writer, "{bo}")?;
                }
                SMIFormat::write_atom_smiles(writer, &frame[current_atom])?;

                // Prevent printing of additional '('
                let mut ring_start = 0;

                // We already know all rings in the structure, so the number of potential rings
                // for the current atom can be printed.
                let any_rings = ring_atoms.get(&current_atom);
                if let Some(&count) = any_rings {
                    for _ in 0..count {
                        ring_count += 1;
                        ring_start += 1;

                        if ring_count >= 10 {
                            write!(writer, "%{ring_count}")?;
                        } else {
                            write!(writer, "{ring_count}")?;
                        }

                        ring_stack.entry(current_atom).or_default().push(ring_count);
                    }
                }

                // Avoid the printing of branch begin/end
                let mut ring_end = 0;

                // Find all ring connections first
                for &neighbor in current_atom_bonds {
                    // avoid 'trivial' rings
                    if neighbor == previous_atom {
                        continue;
                    }

                    // We must have a ring to terminate
                    if written[neighbor] {
                        if let Some(rings) = ring_stack.get_mut(&neighbor) {
                            if !rings.is_empty() {
                                // Remove the first element
                                let ring_num = rings.remove(0);

                                write!(
                                    writer,
                                    "{}",
                                    frame.topology().bond_order(current_atom, neighbor)?
                                )?;

                                // Print the ring index (with "%" if ≥10)
                                if ring_num >= 10 {
                                    write!(writer, "%{ring_num}")?;
                                } else {
                                    write!(writer, "{ring_num}")?;
                                }

                                // If that Vec is now empty, drop the key entirely
                                if rings.is_empty() {
                                    ring_stack.remove(&neighbor);
                                }

                                ring_end += 1;
                            }
                        }
                    }
                }

                // Handle branching
                let mut neighbors_printed = 0;
                for neighbor in current_atom_bonds.iter().rev() {
                    let neighbor = *neighbor;

                    // prevent back tracking
                    if neighbor == previous_atom {
                        continue;
                    }

                    // This got taken care of by printing a ring
                    if written[neighbor] {
                        continue;
                    }

                    // To print a start bracket, we need to be branching (> 2 non-ring bonds)
                    // and we don't want to branch the last neighbor printed
                    let needs_to_branch = neighbors_printed != 0 && neighbors_printed > ring_start;

                    // Depth First Search like recursion
                    atoms_to_process.push((current_atom, neighbor, needs_to_branch));

                    // we printed a neighbor, if there's more than 1 neighbor, then we need to
                    // branch for all but the last neighbor
                    neighbors_printed += 1;
                }

                // End of branch
                if current_atom_bonds.len() - ring_end == 1 && branch_stack != 0 {
                    write!(writer, ")")?;
                    branch_stack -= 1;
                }
            }
            first_atom = false;
        }

        let name = frame.properties.get("name");
        if name.is_some() && name.unwrap().as_string().is_some() {
            write!(writer, "\t{}", name.unwrap().expect_string())?;
        }

        writeln!(writer)?;
        Ok(())
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

    fn finalize(&self, _writer: &mut BufWriter<File>) -> Result<(), CError> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        atom::Atom, bond::BondOrder, frame::Frame, property::Property, trajectory::Trajectory,
    };
    use std::{hint::black_box, path::Path};

    #[test]
    fn check_nsteps() {
        let path = Path::new("./src/tests-data/smi/test.smi");
        let trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 6);
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
        assert_eq!(trajectory.size, 66);

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
        assert_eq!(trajectory.size, 66);

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

    #[test]
    fn write_smi_file() {
        const EXPECTED_CONTENT: &str = r#"C(C)(C)(C)C
C
C~N
C~N(P)=O
C~N(P(#F)$B)=O
C1~N(P(#F:1)$B)=O
C12~N(P(#F:1)$B/2)=O	test
C12(~N(P(#F:1)$B/2)=O)~I	test
C12(~N(P(#F:1)$B/2)(=O)~S)~I	test
[WH5+3].[35Cl-][c:1@H][te@SP3]\[C@@]
O.O.O
"#;
        let named_tmpfile = tempfile::Builder::new()
            .prefix("test-smi")
            .suffix(".smi")
            .tempfile()
            .unwrap();
        let mut trajectory = Trajectory::create(named_tmpfile.path()).unwrap();
        let mut frame = Frame::new();
        for _ in 0..5 {
            frame.add_atom(Atom::new("C".to_string()), [0.0, 0.0, 0.0]);
        }

        frame.add_bond(0, 1, BondOrder::Single).unwrap();
        frame.add_bond(0, 2, BondOrder::Single).unwrap();
        frame.add_bond(0, 3, BondOrder::Single).unwrap();
        frame.add_bond(0, 4, BondOrder::Single).unwrap();

        trajectory.write(&frame).unwrap();

        let mut frame = Frame::new();
        frame.add_atom(Atom::new("C".to_string()), [0.0, 0.0, 0.0]);
        trajectory.write(&frame).unwrap();

        frame.add_atom(Atom::new("N".to_string()), [0.0, 0.0, 0.0]);
        frame.add_bond(0, 1, BondOrder::Unknown).unwrap();
        trajectory.write(&frame).unwrap();

        frame.add_atom(Atom::new("P".to_string()), [0.0, 0.0, 0.0]);
        frame.add_atom(Atom::new("O".to_string()), [0.0, 0.0, 0.0]);
        frame.add_bond(1, 2, BondOrder::Single).unwrap();
        frame.add_bond(1, 3, BondOrder::Double).unwrap();
        trajectory.write(&frame).unwrap();

        frame.add_atom(Atom::new("F".to_string()), [0.0, 0.0, 0.0]);
        frame.add_atom(Atom::new("B".to_string()), [0.0, 0.0, 0.0]);
        frame.add_bond(2, 4, BondOrder::Triple).unwrap();
        frame.add_bond(2, 5, BondOrder::Quadruple).unwrap();
        trajectory.write(&frame).unwrap();

        frame.add_bond(0, 4, BondOrder::Aromatic).unwrap();
        trajectory.write(&frame).unwrap();

        frame.add_bond(0, 5, BondOrder::Up).unwrap();
        frame
            .properties
            .insert("name".to_string(), Property::String("test".to_string()));
        trajectory.write(&frame).unwrap();

        frame.add_atom(Atom::new("I".to_string()), [0.0, 0.0, 0.0]);
        frame.add_bond(0, 6, BondOrder::Unknown).unwrap();
        trajectory.write(&frame).unwrap();

        frame.add_atom(Atom::new("S".to_string()), [0.0, 0.0, 0.0]);
        frame.add_bond(1, 7, BondOrder::Unknown).unwrap();
        trajectory.write(&frame).unwrap();

        // Reinitialize
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("W".to_string()), [0.0, 0.0, 0.0]);
        frame[0].charge = 3.0;
        frame[0]
            .properties
            .insert("hydrogen_count".to_string(), Property::Double(5.0));
        frame[0].properties.insert(
            "chirality".to_string(),
            Property::String("CCW TX99".to_string()),
        );

        frame.add_atom(Atom::new("Cl".to_string()), [0.0, 0.0, 0.0]);
        frame[1].charge = -1.0;
        frame[1].mass = 35.0;
        frame[1]
            .properties
            .insert("hydrogen_count".to_string(), Property::Double(-1.0)); // warning
        frame[1].properties.insert(
            "smiles_class".to_string(),
            Property::String("35-chloride".to_string()),
        ); // warning
        frame[1]
            .properties
            .insert("chirality".to_string(), Property::String("CXX".to_string())); // warning

        frame.add_atom(Atom::new("C".to_string()), [0.0, 0.0, 0.0]);
        frame[2]
            .properties
            .insert("is_aromatic".to_string(), Property::Bool(true));
        frame[2]
            .properties
            .insert("smiles_class".to_string(), Property::Double(1.0));
        frame[2]
            .properties
            .insert("hydrogen_count".to_string(), Property::Double(1.0));
        frame[2]
            .properties
            .insert("chirality".to_string(), Property::String("CCW".to_string()));

        frame.add_atom(Atom::new("Te".to_string()), [0.0, 0.0, 0.0]);
        frame[3]
            .properties
            .insert("is_aromatic".to_string(), Property::Bool(true));
        frame[3].properties.insert(
            "chirality".to_string(),
            Property::String("CCW SP3".to_string()),
        );

        frame.add_atom(Atom::new("C".to_string()), [0.0, 0.0, 0.0]);
        frame[4]
            .properties
            .insert("chirality".to_string(), Property::String("CW".to_string()));

        frame.add_bond(1, 2, BondOrder::Single).unwrap(); // in chemfiles, this was DATIVE_R
        frame.add_bond(2, 3, BondOrder::Single).unwrap(); // in chemfiles, this was DATIVE_L
        frame.add_bond(3, 4, BondOrder::Down).unwrap();

        trajectory.write(&frame).unwrap();

        // Reinitialize and test for discrete molecules
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("O".to_string()), [0.0, 0.0, 0.0]);
        frame.add_atom(Atom::new("O".to_string()), [0.0, 0.0, 0.0]);
        frame.add_atom(Atom::new("O".to_string()), [0.0, 0.0, 0.0]);
        trajectory.write(&frame).unwrap();

        trajectory.finish().unwrap();

        let contents = std::fs::read_to_string(named_tmpfile.path()).unwrap();
        assert_eq!(EXPECTED_CONTENT, contents);
    }
}
