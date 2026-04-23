// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use std::path::Path;
use tng_rs::TngError;

use crate::{error::CError, frame::Frame, topology::Topology, unit_cell::UnitCell};

#[derive(Debug, Default)]
struct TNGFile {
    traj: tng_rs::trajectory::Trajectory,
}

impl TNGFile {
    fn open(path: &Path) -> Result<TNGFile, CError> {
        // TODO: make this inputtable
        let mode = 'r';

        let mut traj = tng_rs::trajectory::Trajectory::new();
        traj.util_trajectory_open(path, mode)
            .map_err(|e| CError::GenericError(e.to_string()))?;

        match mode {
            'r' => {
                if traj.file_headers_read(true).is_err() {
                    traj.util_trajectory_close()
                        .map_err(|e| CError::GenericError(e.to_string()))?;
                    return Err(CError::GenericError(format!(
                        "could not open file at '{}'",
                        path.to_string_lossy()
                    )));
                }
            }
            'w' | 'a' => {
                traj.last_program_name_set("molio");
                // TODO: for now we just use user-molio instead of grabbing the real name
                traj.last_user_name_set("user-molio");
                // TODO: for now we just use host-molio instead of grabbing the real name
                traj.last_computer_name_set("host-molio");

                if mode == 'w' {
                    traj.set_first_program_name("molio");
                    traj.first_user_name_set("molio");
                    traj.first_computer_name_set("host-molio");
                }

                traj.file_headers_write(true)
                    .map_err(|e| CError::GenericError(e.to_string()))?;
            }
            _ => {
                return Err(CError::GenericError(format!(
                    "unsupported mode. Got '{mode}"
                )));
            }
        }

        Ok(TNGFile { traj })
    }
}

impl Drop for TNGFile {
    fn drop(&mut self) {
        self.traj
            .util_trajectory_close()
            .expect("be able to Drop tng trajectory");
    }
}

#[derive(Debug, Default)]
pub struct TNGFormat {
    /// Associated TNG file.
    tng: TNGFile,
    /// Index of the next frame to read.
    index: usize,
    /// Scale factor for all length-dependent data:
    /// positions, velocities, forces, and box shape
    distance_scale_factor: f64,
    /// All the simulation steps in this file.
    simulation_steps: Vec<i64>,
    // As frames can have a varying amount of atoms, this number is frame-dependent
    // Therefore, it can't be set upon opening the file
    natoms: Option<usize>,
}
impl TNGFormat {
    pub(crate) fn open(path: &Path) -> Result<Self, CError> {
        let mut tng = TNGFile::open(path)?;

        let exp = tng.traj.distance_unit_exponential_get();
        // Calculate the scale factor from a given length scale to angstrom
        let distance_scale_factor = 10.0f64.powf(exp as f64 + 10.0);

        // Work around a bug in tng_num_frames_get (https://redmine.gromacs.org/issues/2937). However,
        // that link is dead so it's unclear what the bug is.
        // Manually query all frames in the trajectory, and get the corresponding
        // TNG frame number in `simulation_steps`
        let mut current_step = -1;
        let mut simulation_steps = Vec::new();

        loop {
            // Look for all frames with at least position data.
            let block_ids = [tng_rs::gen_block::BlockID::TrajPositions];
            let result = tng
                .traj
                .util_trajectory_next_frame_present_data_blocks_find(current_step, 1, &block_ids);
            match result {
                Ok((next_step, _n_data_blocks_in_next_frame)) => {
                    current_step = next_step;
                    simulation_steps.push(current_step);
                }
                Err(TngError::Critical(_)) => {
                    return Err(CError::GenericError(
                        "a critical error happened while reading TNG file".to_string(),
                    ));
                }
                Err(_) => {
                    // We found the end of the file
                    break;
                }
            }
        }

        Ok(TNGFormat {
            tng,
            index: 0,
            distance_scale_factor,
            simulation_steps,
            natoms: None,
        })
    }

    pub fn create(_path: &Path) -> Result<Self, CError> {
        unimplemented!("create (and write) is still unimplemented for TNG")
    }

    pub fn read(&mut self) -> Result<Frame, CError> {
        let mut frame = Frame::new();

        frame.properties.insert(
            "simulation_step".into(),
            crate::property::Property::Double(self.simulation_steps[self.index] as f64),
        );
        let natoms = self.tng.traj.num_particles_get() as usize;
        debug_assert!(natoms > 0);
        self.natoms = Some(natoms);
        frame.resize(natoms)?;

        if let Ok(time) = self
            .tng
            .traj
            .util_time_of_frame_get(self.simulation_steps[self.index])
        {
            // TNG stores time in seconds
            // convert to picoseconds
            frame.properties.insert(
                "time".into(),
                crate::property::Property::Double(time * 1e12),
            );
        } else {
            frame
                .properties
                .insert("time".into(), crate::property::Property::Double(0.0));
        }

        self.read_positions(&mut frame)?;
        self.read_velocities(&mut frame)?;
        self.read_cell(&mut frame)?;
        self.read_topology(&mut frame)?;

        self.index += 1;

        Ok(frame)
    }

    pub fn read_at(&mut self, index: usize) -> Result<Frame, CError> {
        self.index = index;
        self.read()
    }

    pub fn write(&mut self, _frame: &Frame) -> Result<(), CError> {
        unimplemented!("write is not yet implemented for TNG files")
    }

    pub fn finish(&mut self) -> Result<(), CError> {
        unimplemented!("write is not yet implemented for TNG files")
    }

    pub fn len(&self) -> Result<usize, CError> {
        Ok(self.simulation_steps.len())
    }

    pub fn is_empty(&self) -> bool {
        self.len().expect("should always be available") == 0
    }

    fn read_topology(&self, frame: &mut Frame) -> Result<(), CError> {
        let mut topology = Topology::default();
        topology.reserve(self.natoms.expect("Some atoms to be present"));

        let moltypes = self.tng.traj.num_molecule_types_get() as usize;
        let molecules_count = self.tng.traj.molecule_cnt_list_get();

        // Read all molecules types
        for (moltype, &molecules_count_moltype) in molecules_count.iter().enumerate().take(moltypes)
        {
            let molecule = self
                .tng
                .traj
                .molecule_of_index_get(moltype as i64)
                .map_err(|e| CError::GenericError(e.to_string()))?;

            // Loop over all the molecules of a given type
            for _ in 0..molecules_count_moltype {
                // For each type, get all the residues in the molecule
                let n_residues = self.tng.traj.molecule_num_residues_get(molecule) as usize;
                for resid in 0..n_residues {
                    let tng_residue = self
                        .tng
                        .traj
                        .molecule_residue_of_index_get(molecule, resid)
                        .map_err(|e| CError::GenericError(e.to_string()))?;
                    let resname = self
                        .tng
                        .traj
                        .residue_name_get(tng_residue, 32)
                        .map_err(|e| CError::GenericError(e.to_string()))?;

                    let mut residue = crate::residue::Residue::new_from_name(resname);

                    // And finally get all the atoms in the residue
                    let n_atoms = self.tng.traj.residue_num_atoms_get(tng_residue) as usize;
                    for id in 0..n_atoms {
                        let tng_atom = self
                            .tng
                            .traj
                            .residue_atom_of_index_get(tng_residue, id)
                            .map_err(|e| CError::GenericError(e.to_string()))?;
                        let name = self
                            .tng
                            .traj
                            .atom_name_get(tng_atom, 32)
                            .map_err(|e| CError::GenericError(e.to_string()))?;
                        let atom_type = self
                            .tng
                            .traj
                            .atom_type_get(tng_atom, 32)
                            .map_err(|e| CError::GenericError(e.to_string()))?;

                        residue.add_atom(topology.len());
                        topology.add_atom(crate::atom::Atom::with_symbol(name, atom_type));
                    }
                    topology.add_residue(residue)?;
                }
            }
        }

        if let Some((n_bonds, from_atoms, to_atoms)) = self.tng.traj.molsystem_bonds_get() {
            for i in 0..n_bonds {
                topology.add_bond(
                    from_atoms[i] as usize,
                    to_atoms[i] as usize,
                    crate::bond::BondOrder::Unknown,
                )?;
            }
        }

        frame.set_topology(topology)?;

        Ok(())
    }

    fn read_cell(&mut self, frame: &mut Frame) -> Result<(), CError> {
        let result = self.tng.traj.util_box_shape_read_range(
            self.simulation_steps[self.index],
            self.simulation_steps[self.index],
        );

        match result {
            Ok(_) => {
                // Continue
            }
            Err(TngError::Critical(_)) => {
                return Err(CError::GenericError(
                    "a critical error happened while reading the cell shape".to_string(),
                ));
            }
            Err(_) => {
                // No unit cell in this frame
                frame.set_unit_cell(UnitCell::new());
                return Ok(());
            }
        }

        let (box_shape, _stride) = result.expect("we should've returned if it was Err");
        let mut matrix: nalgebra::Matrix3<f64> = nalgebra::Matrix3::from_row_slice(&box_shape);
        matrix *= self.distance_scale_factor;
        frame.set_unit_cell(UnitCell::new_from_matrix(matrix)?);

        Ok(())
    }

    fn read_velocities(&mut self, frame: &mut Frame) -> Result<(), CError> {
        let result = self.tng.traj.util_vel_read_range(
            self.simulation_steps[self.index],
            self.simulation_steps[self.index],
        );
        match result {
            Ok(_) => {
                // Continue
            }
            Err(TngError::Critical(_)) => {
                return Err(CError::GenericError(
                    "a critical error happened while reading velocities".to_string(),
                ));
            }
            Err(_) => {
                // No velocity in this frame
                return Ok(());
            }
        }

        frame.add_velocities();
        let (velocities, _stride) = result.expect("we should've returned if it was Err");
        for (i, vel) in frame
            .velocities_mut()
            .expect("we just added velocities")
            .iter_mut()
            .enumerate()
        {
            vel[0] = velocities[3 * i] * self.distance_scale_factor;
            vel[1] = velocities[3 * i + 1] * self.distance_scale_factor;
            vel[2] = velocities[3 * i + 2] * self.distance_scale_factor;
        }

        Ok(())
    }

    fn read_positions(&mut self, frame: &mut Frame) -> Result<(), CError> {
        let (positions, _stride) = self
            .tng
            .traj
            .util_pos_read_range(
                self.simulation_steps[self.index],
                self.simulation_steps[self.index],
            )
            .map_err(|e| CError::GenericError(e.to_string()))?;

        for (i, pos) in frame.positions_mut().iter_mut().enumerate() {
            pos[0] = positions[3 * i] * self.distance_scale_factor;
            pos[1] = positions[3 * i + 1] * self.distance_scale_factor;
            pos[2] = positions[3 * i + 2] * self.distance_scale_factor;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use assert_approx_eq::assert_approx_eq;

    use crate::{bond::Bond, trajectory::Trajectory, unit_cell::CellShape};

    #[test]
    fn read_files() {
        let path = Path::new("./src/tests-data/tng/example.tng");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.len(), 10);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.len(), 15);
        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 10.0, 1e-5);
        assert_approx_eq!(positions[0][1], 10.0, 1e-5);
        assert_approx_eq!(positions[0][2], 10.0, 1e-5);
        assert_approx_eq!(positions[11][0], 85.0, 1e-5);
        assert_approx_eq!(positions[11][1], 330.0, 1e-5);
        assert_approx_eq!(positions[11][2], 340.0, 1e-5);

        let cell = frame.cell();
        assert_eq!(cell.shape, CellShape::Infinite);

        // Skip a frame
        let _ = trajectory.read().unwrap().unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        assert_eq!(frame.len(), 15);
        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 10.1562, 1e-4);
        assert_approx_eq!(positions[0][1], 10.2344, 1e-4);
        assert_approx_eq!(positions[0][2], 10.3125, 1e-4);
        assert_approx_eq!(positions[11][0], 85.0, 1e-5);
        assert_approx_eq!(positions[11][1], 330.0, 1e-5);
        assert_approx_eq!(positions[11][2], 340.0, 1e-5);
    }

    #[test]
    fn read_cell() {
        let path = Path::new("./src/tests-data/tng/water.tng");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.len(), 29700);

        let cell = frame.cell();
        assert_eq!(cell.shape, CellShape::Orthorhombic);
        assert_approx_eq!(cell.lengths()[0], 15.0);
        assert_approx_eq!(cell.lengths()[1], 15.0);
        assert_approx_eq!(cell.lengths()[2], 15.0);
    }

    #[test]
    fn read_triclinic_cell() {
        let path = Path::new("./src/tests-data/tng/1vln-triclinic.tng");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.len(), 14520);

        let cell = frame.cell();
        assert_eq!(cell.shape, CellShape::Triclinic);

        assert_approx_eq!(cell.lengths()[0], 78.8, 1e-5);
        assert_approx_eq!(cell.lengths()[1], 79.3, 1e-5);
        assert_approx_eq!(cell.lengths()[2], 133.3, 1e-5);

        assert_approx_eq!(cell.angles()[0], 97.1, 1e-5);
        assert_approx_eq!(cell.angles()[1], 90.2, 1e-5);
        assert_approx_eq!(cell.angles()[2], 97.5, 1e-5);
    }

    #[test]
    fn read_velocities() {
        let path = Path::new("./src/tests-data/tng/1aki.tng");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.len(), 6);
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.len(), 38376);
        assert_eq!(frame.get_frame_index(), 0);
        assert_eq!(
            frame
                .properties
                .get("simulation_step")
                .unwrap()
                .as_double()
                .unwrap(),
            0.0
        );
        assert_approx_eq!(
            frame.properties.get("time").unwrap().as_double().unwrap(),
            0.0,
            1e-4
        );

        let cell = frame.cell();
        assert_eq!(cell.shape, CellShape::Orthorhombic);
        assert_approx_eq!(cell.lengths()[0], 73.39250);
        assert_approx_eq!(cell.lengths()[1], 73.39250);
        assert_approx_eq!(cell.lengths()[2], 73.39250);

        let velocities = frame.velocities().unwrap();
        assert_approx_eq!(velocities[450][0], -1.44889, 1e-4);
        assert_approx_eq!(velocities[450][1], 6.50066e-1, 1e-4);
        assert_approx_eq!(velocities[450][2], -7.64032, 1e-4);
        assert_approx_eq!(velocities[4653][0], -16.5949, 1e-4);
        assert_approx_eq!(velocities[4653][1], -4.62240, 1e-4);
        assert_approx_eq!(velocities[4653][2], -7.01133, 1e-4);

        let frame = trajectory.read_at(5).unwrap().unwrap();
        assert_eq!(frame.len(), 38376);
        assert_eq!(frame.get_frame_index(), 5);
        assert_eq!(
            frame
                .properties
                .get("simulation_step")
                .unwrap()
                .as_double()
                .unwrap(),
            50.0
        );
        assert_approx_eq!(
            frame.properties.get("time").unwrap().as_double().unwrap(),
            0.1,
            1e-4
        );

        let cell = frame.cell();
        assert_eq!(cell.shape, CellShape::Orthorhombic);
        assert_approx_eq!(cell.lengths()[0], 73.39250);
        assert_approx_eq!(cell.lengths()[1], 73.39250);
        assert_approx_eq!(cell.lengths()[2], 73.39250);

        let velocities = frame.velocities().unwrap();
        assert_approx_eq!(velocities[450][0], 8.23913, 1e-4);
        assert_approx_eq!(velocities[450][1], 2.99123, 1e-4);
        assert_approx_eq!(velocities[450][2], 10.5270, 1e-4);
        assert_approx_eq!(velocities[4653][0], -48.8318, 1e-4);
        assert_approx_eq!(velocities[4653][1], -5.90270, 1e-4);
        assert_approx_eq!(velocities[4653][2], -6.86679, 1e-4);
    }

    #[test]
    fn read_topology() {
        let path = Path::new("./src/tests-data/tng/example.tng");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();
        let topology = frame.topology();

        assert_eq!(topology.len(), 15);
        assert_eq!(topology[0].name, "O");
        assert_eq!(topology[0].symbol, "O");
        assert_eq!(topology[1].name, "HO1");
        assert_eq!(topology[1].symbol, "H");
        assert_eq!(topology[2].name, "HO2");
        assert_eq!(topology[2].symbol, "H");

        assert_eq!(topology.residues.len(), 5);
        let residue = &topology.residues[0];
        assert_eq!(residue.len(), 3);
        assert!(residue.contains(0));
        assert!(residue.contains(1));
        assert!(residue.contains(2));

        let bonds = topology.bonds();
        let expected = vec![
            [0, 1],
            [0, 2],
            [3, 4],
            [3, 5],
            [6, 7],
            [6, 8],
            [9, 10],
            [9, 11],
            [12, 13],
            [12, 14],
        ];
        assert_eq!(bonds.len(), expected.len());
        for bond in &expected {
            assert!(bonds.contains(&Bond::new(bond[0], bond[1])));
        }
    }

    #[test]
    fn non_consecutive_frame_indices() {
        // cf https://github.com/chemfiles/chemfiles/issues/242
        let path = Path::new("./src/tests-data/tng/cobrotoxin.tng");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.len(), 3);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.get_frame_index(), 0);
        assert_eq!(
            frame
                .properties
                .get("simulation_step")
                .unwrap()
                .as_double()
                .unwrap(),
            0.0
        );
        assert_approx_eq!(
            frame.properties.get("time").unwrap().as_double().unwrap(),
            0.0,
            1e-4
        );

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.get_frame_index(), 1);
        assert_eq!(
            frame
                .properties
                .get("simulation_step")
                .unwrap()
                .as_double()
                .unwrap(),
            25000.0
        );
        assert_approx_eq!(
            frame.properties.get("time").unwrap().as_double().unwrap(),
            50.0,
            1e-4
        );

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.get_frame_index(), 2);
        assert_eq!(
            frame
                .properties
                .get("simulation_step")
                .unwrap()
                .as_double()
                .unwrap(),
            50000.0
        );
        assert_approx_eq!(
            frame.properties.get("time").unwrap().as_double().unwrap(),
            100.0,
            1e-4
        );

        let frame = trajectory.read_at(0).unwrap().unwrap();
        assert_eq!(frame.len(), 19385);
        let positions = frame.positions();

        assert_approx_eq!(positions[5569][0], 14.94, 1e-5);
        assert_approx_eq!(positions[5569][1], 4.03, 1e-5);
        assert_approx_eq!(positions[5569][2], 19.89, 1e-5);

        assert_approx_eq!(positions[11675][0], 44.75, 1e-5);
        assert_approx_eq!(positions[11675][1], 16.05, 1e-5);
        assert_approx_eq!(positions[11675][2], 6.1, 1e-5);
    }
}
