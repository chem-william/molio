// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::{error::CError, frame::Frame, property::Property, unit_cell::UnitCell};
use core::f64;
use log::{debug, warn};
use netcdf3::{DataSet, DataType, FileReader, FileWriter, Variable};
use std::path::{Path, PathBuf};

#[derive(Debug)]
struct VariableWithScale {
    var: Variable,
    scale: f64,
}
#[derive(Debug, Default)]
struct Variables {
    coordinates: Option<VariableWithScale>,
    velocities: Option<VariableWithScale>,
    cell_lengths: Option<VariableWithScale>,
    cell_angles: Option<VariableWithScale>,
    time: Option<VariableWithScale>,
}

#[derive(Debug)]
struct BufferedFrame {
    positions: Vec<f32>,
    velocities: Option<Vec<f32>>,
    cell_lengths: [f32; 3],
    cell_angles: [f32; 3],
    time: Option<f32>,
}

/// AMBER does not use the text codec traits. It is a dedicated binary
/// format entry point with its own owned state.
#[derive(Debug, Default)]
pub struct AMBERTrajFormat {
    reader: Option<FileReader>,
    variables: Variables,
    index: usize,
    file_title: String,
    n_atoms: usize,
    // Writing state
    write_path: Option<PathBuf>,
    buffered_frames: Vec<BufferedFrame>,
    has_velocities: bool,
    initialized: bool,
}

fn read_array(
    file_reader: &mut FileReader,
    index: usize,
    variable: &VariableWithScale,
    array: &mut [[f64; 3]],
) -> Result<(), CError> {
    match variable.var.data_type() {
        DataType::I8 | DataType::U8 | DataType::I16 | DataType::I32 => {
            return Err(CError::GenericError(
                "invalid type for variable, expected f32/f64 data".to_string(),
            ));
        }
        DataType::F32 => {
            let buffer = file_reader.read_record_f32(variable.var.name(), index)?;
            for (idx, item) in array.iter_mut().enumerate() {
                item[0] = variable.scale * f64::from(buffer[3 * idx + 0]);
                item[1] = variable.scale * f64::from(buffer[3 * idx + 1]);
                item[2] = variable.scale * f64::from(buffer[3 * idx + 2]);
            }
        }
        DataType::F64 => {
            let buffer = file_reader.read_record_f64(variable.var.name(), index)?;
            for (idx, item) in array.iter_mut().enumerate() {
                item[0] = variable.scale * buffer[3 * idx + 0];
                item[1] = variable.scale * buffer[3 * idx + 1];
                item[2] = variable.scale * buffer[3 * idx + 2];
            }
        }
    };

    Ok(())
}

fn validate_common(file_reader: &FileReader, convention: &str) -> Result<(), CError> {
    let data_set = file_reader.data_set();

    fn check_attr(data_set: &DataSet, attr: &str, expected_attr: &str) -> Result<(), CError> {
        let read_attr = data_set.get_global_attr(attr).ok_or(CError::GenericError(
            "expected a 'Conventions' attribute to be defined".to_string(),
        ))?;

        let read_attr = read_attr.get_as_string();

        if let Some(read_attr) = read_attr {
            assert_eq!(read_attr, expected_attr, "expected '{expected_attr}'");
        } else {
            return Err(CError::GenericError(
                "could not read attr as string".to_string(),
            ));
        };

        Ok(())
    }
    check_attr(data_set, "Conventions", convention)?;
    check_attr(data_set, "ConventionVersion", "1.0")?;

    if let Some(spatial) = file_reader.data_set().get_dim("spatial") {
        assert_eq!(
            spatial.size(),
            3,
            "{}",
            format!(
                "'spatial' dimension must have a size of 3, got {}",
                spatial.size()
            )
        );
    } else {
        return Err(CError::GenericError(format!("missing 'spatial' dimension")));
    }

    file_reader
        .data_set()
        .get_dim("atom")
        .ok_or(CError::GenericError("missing 'atom' dimension".to_string()))?;

    if let Some(cell_spatial) = file_reader.data_set().get_dim("cell_spatial") {
        assert_eq!(
            cell_spatial.size(),
            3,
            "{}",
            format!(
                "'cell_spatial' dimension must have a size of 3, got {}",
                cell_spatial.size()
            )
        );
    };
    if let Some(cell_angular) = file_reader.data_set().get_dim("cell_angular") {
        assert_eq!(
            cell_angular.size(),
            3,
            "{}",
            format!(
                "'cell_angular' dimension must have a size of 3, got {}",
                cell_angular.size()
            )
        );
    };

    Ok(())
}

fn scale_for_distance(units: &str) -> f64 {
    match units.to_lowercase().as_str() {
        "angstroms" | "angstrom" | "a" => 1.0,
        "meters" | "meter" | "m" => 1e10,
        "centimeters" | "centimeter" | "cm" => 1e8,
        "micrometers" | "micrometer" | "µm" | "um" => 1e4,
        "nanometers" | "nanometer" | "nm" => 1e4,
        "bohrs" | "bohr" => 0.52918,
        unknown => {
            warn!("unknown unit ({unknown}) for distances");
            1.0
        }
    }
}

fn scale_for_velocity(units: &str) -> f64 {
    let lower_case = units.to_lowercase();

    let mut splitted = lower_case.split('/');
    if let Some(scale) = splitted.next() {
        if let Some(time_unit) = splitted.next().as_mut() {
            let scale = scale_for_distance(scale);
            let time_unit = scale_for_time(*time_unit);
            return scale / time_unit;
        };
    }
    warn!("unknown unit ({units}) for velocities");
    1.0
}

fn scale_for_time(units: &str) -> f64 {
    match units.to_lowercase().as_str() {
        // nothing to do for picoseconds
        "picoseconds" | "picosecond" | "ps" => 1.0,
        "femtoseconds" | "femtosecond" | "fs" => 1e-3,
        "nanoseconds" | "nanosecond" | "ns" => 1e3,
        "microseconds" | "microsecond" | "ms" => 1e6,
        "seconds" | "second" | "s" => 1e12,
        unknown => {
            warn!("unknown unit ({unknown}) for time");
            1.0
        }
    }
}

fn scale_factor(variable: &Variable) -> Result<f64, CError> {
    if let Some(scale) = variable
        .get_attr_f32("scale_factor")
        .and_then(|scale| scale.first().copied())
    {
        return Ok(f64::from(scale));
    }

    if let Some(scale) = variable
        .get_attr_f64("scale_factor")
        .and_then(|scale| scale.first().copied())
    {
        return Ok(scale);
    }

    if variable.get_attr("scale_factor").is_some() {
        return Err(CError::GenericError(format!(
            "scale_factor attribute for '{}' must be a floating point value",
            variable.name()
        )));
    }

    Ok(1.0)
}

impl AMBERTrajFormat {
    fn read_cell(&mut self) -> Result<Option<UnitCell>, CError> {
        if self.variables.cell_lengths.is_none() || self.variables.cell_angles.is_none() {
            return Ok(None);
        };
        let reader = self
            .reader
            .as_mut()
            .expect("reader should've been initialized");

        let mut lengths = [0.0; 3];
        let cell_lengths = self
            .variables
            .cell_lengths
            .as_ref()
            .expect("we just checked for None");
        match self
            .variables
            .cell_lengths
            .as_ref()
            .expect("we just checked for None")
            .var
            .data_type()
        {
            DataType::I8 | DataType::U8 | DataType::I16 | DataType::I32 => unreachable!(),
            DataType::F32 => {
                let buffer = reader.read_record_f32("cell_lengths", self.index)?;

                lengths[0] = cell_lengths.scale * f64::from(buffer[0]);
                lengths[1] = cell_lengths.scale * f64::from(buffer[1]);
                lengths[2] = cell_lengths.scale * f64::from(buffer[2]);
            }
            DataType::F64 => {
                let buffer = reader.read_record_f64("cell_lengths", self.index)?;
                lengths[0] = cell_lengths.scale * buffer[0];
                lengths[1] = cell_lengths.scale * buffer[1];
                lengths[2] = cell_lengths.scale * buffer[2];
            }
        };

        let mut angles = [0.0; 3];
        let cell_angles = self
            .variables
            .cell_angles
            .as_ref()
            .expect("we just checked for None");
        match self
            .variables
            .cell_angles
            .as_ref()
            .expect("we just checked for None")
            .var
            .data_type()
        {
            DataType::I8 | DataType::U8 | DataType::I16 | DataType::I32 => unreachable!(),
            DataType::F32 => {
                let buffer = reader.read_record_f32("cell_angles", self.index)?;

                angles[0] = cell_angles.scale * f64::from(buffer[0]);
                angles[1] = cell_angles.scale * f64::from(buffer[1]);
                angles[2] = cell_angles.scale * f64::from(buffer[2]);
            }
            DataType::F64 => {
                let buffer = reader.read_record_f64("cell_angles", self.index)?;
                angles[0] = cell_angles.scale * buffer[0];
                angles[1] = cell_angles.scale * buffer[1];
                angles[2] = cell_angles.scale * buffer[2];
            }
        };

        Ok(Some(UnitCell::new_from_lengths_angles(
            lengths,
            &mut angles,
        )?))
    }

    pub fn open(path: &Path) -> Result<Self, CError> {
        let mut file_reader = FileReader::open(path).unwrap();

        validate_common(&file_reader, "AMBER")?;

        let file_title =
            if let Some(file_title) = file_reader.data_set().get_global_attr_as_string("title") {
                file_title
            } else {
                "".to_string()
            };

        let n_atoms = file_reader
            .data_set()
            .get_dim("atom")
            .expect("we just just validated that 'atom' exists")
            .size();

        // Get the variables actually defined in the file
        let coordinates = if let Some(coordinates) = file_reader.data_set().get_var("coordinates") {
            let mut coords = Some(VariableWithScale {
                var: coordinates.clone(),
                scale: scale_factor(&coordinates)?,
            });
            if let Some(units) = coordinates.get_attr_as_string("units") {
                coords.as_mut().expect("we just init'ed").scale *=
                    scale_for_distance(units.as_str());
            } else {
                debug!("not scaling coordinates");
            };
            coords
        } else {
            warn!("the coordinates variable is not defined in this file.");
            None
        };

        let velocities = if let Some(velocities) = file_reader.data_set().get_var("velocities") {
            let mut velo = Some(VariableWithScale {
                var: velocities.clone(),
                scale: scale_factor(&velocities)?,
            });
            if let Some(units) = velocities.get_attr_as_string("units") {
                velo.as_mut().expect("we just init'ed").scale *= scale_for_velocity(units.as_str());
            } else {
                debug!("not scaling velocities");
            };
            velo
        } else {
            warn!("the velocities variable is not defined in this file.");
            None
        };

        let read_cell_lengths = file_reader.data_set().get_var("cell_lengths");
        let read_cell_angles = file_reader.data_set().get_var("cell_angles");
        let mut cell_lengths = None;
        let mut cell_angles = None;
        if let Some(read_cell_angles) = read_cell_angles
            && let Some(read_cell_lengths) = read_cell_lengths
        {
            cell_lengths = Some(VariableWithScale {
                var: read_cell_lengths.clone(),
                scale: scale_factor(&read_cell_lengths)?,
            });
            if let Some(units) = read_cell_lengths.get_attr_as_string("units") {
                cell_lengths.as_mut().expect("we just init'ed it").scale *=
                    scale_for_distance(units.as_str());
            } else {
                debug!("not scaling cell lengths");
            };

            cell_angles = Some(VariableWithScale {
                var: read_cell_angles.clone(),
                scale: scale_factor(&read_cell_angles)?,
            });
            if let Some(units) = read_cell_angles.get_attr_as_string("units") {
                let scaling_factor = match units.to_lowercase().as_str() {
                    "" | "degrees" | "degree" => 1.0,
                    "radians" | "radian" => 180.0 / f64::consts::PI,
                    unknown => {
                        warn!("unknown unit ({unknown}) for angles");
                        1.0
                    }
                };
                cell_angles.as_mut().expect("we just init'ed it").scale *= scaling_factor;
            } else {
                debug!("not scaling cell angles");
            };
        } else if read_cell_lengths.is_some() {
            if read_cell_angles.is_none() {
                return Err(CError::GenericError(
                    "cell_lengths requires cell_angles to be defined".to_string(),
                ));
            }
            if coordinates.is_none() {
                return Err(CError::GenericError(
                    "cell_lengths requires coordinates to be defined.".to_string(),
                ));
            }
        };

        let time = if let Some(time) = file_reader.data_set().get_var("time") {
            let mut t = Some(VariableWithScale {
                var: time.clone(),
                scale: scale_factor(&time)?,
            });
            if let Some(units) = time.get_attr_as_string("units") {
                t.as_mut().expect("we just init'ed").scale *= scale_for_time(units.as_str());
            } else {
                debug!("not scaling time");
            };
            t
        } else {
            warn!("the timevariable is not defined in this file.");
            None
        };

        let mut array = vec![[0.0; 3]; n_atoms];
        read_array(
            &mut file_reader,
            0,
            coordinates.as_ref().unwrap(),
            array.as_mut_slice(),
        )?;

        let variables = Variables {
            coordinates,
            velocities,
            cell_lengths,
            cell_angles,
            time,
        };
        let format = AMBERTrajFormat {
            reader: Some(file_reader),
            index: 0,
            file_title,
            variables,
            n_atoms,
            ..Default::default()
        };

        // TODO: append mode, line 143 chemfiles/src/formats/AmberNetCDF.cpp

        Ok(format)
    }

    pub fn create(path: &Path) -> Result<Self, CError> {
        Ok(Self {
            write_path: Some(path.to_path_buf()),
            ..Default::default()
        })
    }

    pub fn read(&mut self) -> Result<Frame, CError> {
        let frame = self.read_at(self.index);
        self.index += 1;

        frame
    }

    pub fn read_at(&mut self, index: usize) -> Result<Frame, CError> {
        self.index = index;

        let mut frame = Frame::new();
        if let Some(unitcell) = self.read_cell()? {
            frame.set_unitcell(unitcell);
        };

        if !self.file_title.is_empty() {
            frame.properties.insert(
                "name".to_string(),
                crate::property::Property::String(self.file_title.clone()),
            );
        }

        frame.resize(self.n_atoms)?;

        let reader = self
            .reader
            .as_mut()
            .expect("reader should've been initialized");
        if let Some(coordinates) = self.variables.coordinates.as_ref() {
            read_array(reader, self.index, &coordinates, frame.positions_mut())?;
        }

        if let Some(velocities) = self.variables.velocities.as_ref() {
            frame.add_velocities();
            read_array(
                reader,
                self.index,
                &velocities,
                frame.velocities_mut().expect("just resized velocities"),
            )?;
        }

        if let Some(time) = self.variables.time.as_ref() {
            let time_value;
            match time.var.data_type() {
                DataType::I8 | DataType::U8 | DataType::I16 | DataType::I32 => {
                    return Err(CError::GenericError(
                        "invalid type for time variable".to_string(),
                    ));
                }
                DataType::F32 => {
                    let value = reader.read_record_f32("time", self.index)?;
                    time_value = time.scale * f64::from(value[0]);
                }
                DataType::F64 => {
                    let value = reader.read_record_f64("time", self.index)?;
                    time_value = time.scale * value[0];
                }
            }
            frame
                .properties
                .insert("time".to_string(), Property::Double(time_value));
        }

        Ok(frame)
    }

    pub fn len(&self) -> Result<usize, CError> {
        self.reader
            .as_ref()
            .expect("we should have init reader")
            .data_set()
            .num_records()
            .ok_or(CError::GenericError(
                "should not have unlimited-size dimension".to_string(),
            ))
    }

    pub fn is_empty(&self) -> Result<bool, CError> {
        self.len().map(|n| n == 0)
    }

    pub fn write(&mut self, frame: &Frame) -> Result<(), CError> {
        if !self.initialized {
            self.n_atoms = frame.size();
            self.has_velocities = frame.velocities().is_some();
            if let Some(Property::String(name)) = frame.properties.get("name") {
                self.file_title = name.clone();
            }
            self.initialized = true;
        }

        if frame.size() != self.n_atoms {
            return Err(CError::GenericError(format!(
                "this file can only write frames with {} atoms, but the frame contains {} atoms",
                self.n_atoms,
                frame.size()
            )));
        }

        let mut positions = Vec::with_capacity(3 * self.n_atoms);
        for pos in frame.positions() {
            positions.push(pos[0] as f32);
            positions.push(pos[1] as f32);
            positions.push(pos[2] as f32);
        }

        let velocities = if self.has_velocities {
            frame.velocities().map(|vels| {
                let mut v = Vec::with_capacity(3 * self.n_atoms);
                for vel in vels {
                    v.push(vel[0] as f32);
                    v.push(vel[1] as f32);
                    v.push(vel[2] as f32);
                }
                v
            })
        } else {
            None
        };

        let cell = frame.cell();
        let lengths = cell.lengths();
        let angles = cell.angles();

        let time = frame
            .properties
            .get("time")
            .and_then(|p| p.as_double())
            .map(|t| t as f32);

        self.buffered_frames.push(BufferedFrame {
            positions,
            velocities,
            cell_lengths: [lengths[0] as f32, lengths[1] as f32, lengths[2] as f32],
            cell_angles: [angles[0] as f32, angles[1] as f32, angles[2] as f32],
            time,
        });

        self.index += 1;
        Ok(())
    }

    pub fn finish(&mut self) -> Result<(), CError> {
        let path = match self.write_path.as_ref() {
            Some(p) => p,
            None => return Ok(()),
        };

        if self.buffered_frames.is_empty() {
            return Ok(());
        }

        let num_frames = self.buffered_frames.len();
        let has_time = self.buffered_frames.iter().any(|f| f.time.is_some());

        // Build the DataSet locally — ephemeral, like the C++ builder
        let mut data_set = DataSet::new();
        data_set.add_global_attr_string("Conventions", "AMBER")?;
        data_set.add_global_attr_string("ConventionVersion", "1.0")?;
        data_set.add_global_attr_string("program", "molio")?;
        data_set.add_global_attr_string("programVersion", env!("CARGO_PKG_VERSION"))?;

        if !self.file_title.is_empty() {
            data_set.add_global_attr_string("title", &self.file_title)?;
        }

        data_set.set_unlimited_dim("frame", num_frames)?;
        data_set.add_fixed_dim("spatial", 3)?;
        data_set.add_fixed_dim("atom", self.n_atoms)?;
        data_set.add_fixed_dim("cell_spatial", 3)?;
        data_set.add_fixed_dim("cell_angular", 3)?;
        data_set.add_fixed_dim("label", 5)?;

        // Fixed variables
        data_set.add_var_u8("spatial", &["spatial"])?;
        data_set.add_var_u8("cell_spatial", &["cell_spatial"])?;
        data_set.add_var_u8("cell_angular", &["cell_angular", "label"])?;

        // Record variables
        if has_time {
            data_set.add_var_f32("time", &["frame"])?;
            data_set.add_var_attr_string("time", "units", "picosecond")?;
        }

        data_set.add_var_f32("coordinates", &["frame", "atom", "spatial"])?;
        data_set.add_var_attr_string("coordinates", "units", "angstrom")?;

        data_set.add_var_f32("cell_lengths", &["frame", "cell_spatial"])?;
        data_set.add_var_attr_string("cell_lengths", "units", "angstrom")?;

        data_set.add_var_f32("cell_angles", &["frame", "cell_angular"])?;
        data_set.add_var_attr_string("cell_angles", "units", "degree")?;

        if self.has_velocities {
            data_set.add_var_f32("velocities", &["frame", "atom", "spatial"])?;
            data_set.add_var_attr_string("velocities", "units", "angstrom/picosecond")?;
        }

        // Create FileWriter locally — both data_set and writer live in this scope
        let mut writer = FileWriter::open(path)?;
        writer.set_def(&data_set, netcdf3::Version::Offset64Bit, 0)?;

        // Write fixed variables
        writer.write_var_u8("spatial", b"xyz")?;
        writer.write_var_u8("cell_spatial", b"abc")?;
        writer.write_var_u8("cell_angular", b"alphabeta gamma")?;

        // Write record variables
        for (i, frame) in self.buffered_frames.iter().enumerate() {
            writer.write_record_f32("coordinates", i, &frame.positions)?;

            writer.write_record_f32("cell_lengths", i, &frame.cell_lengths)?;
            writer.write_record_f32("cell_angles", i, &frame.cell_angles)?;

            if let Some(ref velocities) = frame.velocities {
                writer.write_record_f32("velocities", i, velocities)?;
            }

            if let Some(time) = frame.time {
                writer.write_record_f32("time", i, &[time])?;
            }
        }

        writer.close()?;
        self.buffered_frames.clear();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use assert_approx_eq::assert_approx_eq;
    use tempfile::Builder;

    use crate::{
        frame::Frame,
        trajectory::Trajectory,
        unit_cell::{CellShape, UnitCell},
    };

    #[test]
    fn read_files_in_netcdf_format() {
        let path = Path::new("./src/tests-data/netcdf/water.nc");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        assert_eq!(frame.size(), 297);

        assert_eq!(frame.properties.get("name"), None);

        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 0.4172191, 1e-4);
        assert_approx_eq!(positions[0][1], 8.303366, 1e-4);
        assert_approx_eq!(positions[0][2], 11.73717, 1e-4);

        assert_approx_eq!(positions[296][0], 6.664049, 1e-4);
        assert_approx_eq!(positions[296][1], 11.61418, 1e-4);
        assert_approx_eq!(positions[296][2], 12.96149, 1e-4);

        assert_approx_eq!(
            frame.properties.get("time").unwrap().as_double().unwrap(),
            2.02
        );
    }

    #[test]
    fn read_more_than_one_frame() {
        let path = Path::new("./src/tests-data/netcdf/water.nc");
        let mut trajectory = Trajectory::open(path).unwrap();
        trajectory.read().unwrap().unwrap();
        trajectory.read().unwrap().unwrap();
        let mut frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 297);
        assert_eq!(frame.properties.get("name"), None);

        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 0.2990952, 1e-4);
        assert_approx_eq!(positions[0][1], 8.31003, 1e-4);
        assert_approx_eq!(positions[0][2], 11.72146, 1e-4);

        assert_approx_eq!(positions[296][0], 6.797599, 1e-4);
        assert_approx_eq!(positions[296][1], 11.50882, 1e-4);
        assert_approx_eq!(positions[296][2], 12.70423, 1e-4);

        assert_approx_eq!(
            frame.properties.get("time").unwrap().as_double().unwrap(),
            2.04
        );

        while let Some(next_frame) = trajectory.read().unwrap() {
            frame = next_frame;
        }

        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 0.3185586, 1e-4);
        assert_approx_eq!(positions[0][1], 8.776042, 1e-4);
        assert_approx_eq!(positions[0][2], 11.8927, 1e-4);

        assert_approx_eq!(positions[296][0], 7.089802, 1e-4);
        assert_approx_eq!(positions[296][1], 10.35007, 1e-4);
        assert_approx_eq!(positions[296][2], 12.8159, 1e-4);

        assert_approx_eq!(
            frame.properties.get("time").unwrap().as_double().unwrap(),
            3.01
        );
    }

    #[test]
    fn no_cell() {
        let path = Path::new("./src/tests-data/netcdf/no-cell.nc");
        let mut trajectory = Trajectory::open(path).unwrap();
        let frame = trajectory.read().unwrap().unwrap();

        assert_eq!(frame.size(), 1989);
        assert_eq!(
            frame.properties.get("name").unwrap().as_string().unwrap(),
            "Cpptraj Generated trajectory"
        );
        assert_eq!(*frame.cell(), UnitCell::new());
    }

    #[test]
    fn scale_factor() {
        let path = Path::new("./src/tests-data/netcdf/scaled_traj.nc");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.size, 26);

        let frame = trajectory.read_at(12).unwrap().unwrap();
        assert_eq!(frame.size(), 1938);
        assert_eq!(frame.properties.get("name"), None);

        let cell = frame.cell();
        assert_eq!(cell.shape, CellShape::Orthorhombic);
        assert_approx_eq!(cell.lengths()[0], 1.765 * 60.9682, 1e-4);
        assert_approx_eq!(cell.lengths()[1], 1.765 * 60.9682, 1e-4);
        assert_approx_eq!(cell.lengths()[2], 1.765 * 0.0, 1e-4);

        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 1.39 * 0.455, 1e-4);
        assert_approx_eq!(positions[0][1], 1.39 * 0.455, 1e-4);
        assert_approx_eq!(positions[0][2], 0.0 * 0.455, 1e-4);
        assert_approx_eq!(positions[296][0], 29.1 * 0.455, 1e-4);
        assert_approx_eq!(positions[296][1], 37.41 * 0.455, 1e-4);
        assert_approx_eq!(positions[296][2], 0.0 * 0.455, 1e-4);

        let velocities = frame.velocities().unwrap();
        assert_approx_eq!(velocities[1400][0], 0.6854072 * -0.856, 1e-4);
        assert_approx_eq!(velocities[1400][1], 0.09196011 * -0.856, 1e-4);
        assert_approx_eq!(velocities[1400][2], 2.260214 * -0.856, 1e-4);

        assert_approx_eq!(velocities[1600][0], -0.3342645 * -0.856, 1e-4);
        assert_approx_eq!(velocities[1600][1], 0.322594 * -0.856, 1e-4);
        assert_approx_eq!(velocities[1600][2], -2.446901 * -0.856, 1e-4);
    }

    fn check_frame(frame: &Frame) {
        let name = frame.properties.get("name").unwrap().as_string().unwrap();
        assert_eq!(name, "Test Title 123");

        let positions = frame.positions();
        assert_approx_eq!(positions[0][0], 0.0, 1e-6);
        assert_approx_eq!(positions[0][1], 0.0, 1e-6);
        assert_approx_eq!(positions[0][2], 0.0, 1e-6);

        assert_approx_eq!(positions[1][0], 1.0, 1e-6);
        assert_approx_eq!(positions[1][1], 2.0, 1e-6);
        assert_approx_eq!(positions[1][2], 3.0, 1e-6);

        assert_approx_eq!(positions[2][0], 2.0, 1e-6);
        assert_approx_eq!(positions[2][1], 4.0, 1e-6);
        assert_approx_eq!(positions[2][2], 6.0, 1e-6);

        assert_approx_eq!(positions[3][0], 3.0, 1e-6);
        assert_approx_eq!(positions[3][1], 6.0, 1e-6);
        assert_approx_eq!(positions[3][2], 9.0, 1e-6);

        let velocities = frame.velocities().unwrap();
        assert_approx_eq!(velocities[0][0], -3.0, 1e-6);
        assert_approx_eq!(velocities[0][1], -2.0, 1e-6);
        assert_approx_eq!(velocities[0][2], -1.0, 1e-6);

        assert_approx_eq!(velocities[1][0], -3.0, 1e-6);
        assert_approx_eq!(velocities[1][1], -2.0, 1e-6);
        assert_approx_eq!(velocities[1][2], -1.0, 1e-6);

        assert_approx_eq!(velocities[2][0], -3.0, 1e-6);
        assert_approx_eq!(velocities[2][1], -2.0, 1e-6);
        assert_approx_eq!(velocities[2][2], -1.0, 1e-6);

        assert_approx_eq!(velocities[3][0], -3.0, 1e-6);
        assert_approx_eq!(velocities[3][1], -2.0, 1e-6);
        assert_approx_eq!(velocities[3][2], -1.0, 1e-6);

        let cell = frame.cell();
        let lengths = cell.lengths();
        let angles = cell.angles();
        assert_approx_eq!(lengths[0], 2.0, 1e-6);
        assert_approx_eq!(lengths[1], 3.0, 1e-6);
        assert_approx_eq!(lengths[2], 4.0, 1e-6);

        assert_approx_eq!(angles[0], 80.0, 1e-6);
        assert_approx_eq!(angles[1], 90.0, 1e-6);
        assert_approx_eq!(angles[2], 120.0, 1e-6);
    }

    #[test]
    fn write_files_in_netcdf_format() {
        let mut frame = Frame::from_unitcell(
            UnitCell::new_from_lengths_angles([2.0, 3.0, 4.0], &mut [80.0, 90.0, 120.0]).unwrap(),
        );
        frame.properties.insert(
            "name".to_string(),
            crate::property::Property::String("Test Title 123".to_string()),
        );
        frame
            .properties
            .insert("time".to_string(), crate::property::Property::Double(2.0));
        frame.add_velocities();
        for i in 0..4 {
            frame.add_atom_with_velocity(
                crate::atom::Atom::new("X".to_string()),
                [1.0 * i as f64, 2.0 * i as f64, 3.0 * i as f64],
                [-3.0, -2.0, -1.0],
            );
        }
        let named_tmpfile = Builder::new()
            .prefix("netcdf-test-write")
            .suffix(".nc")
            .tempfile()
            .unwrap();
        let mut trajectory = Trajectory::create(named_tmpfile.path()).unwrap();
        trajectory.write(&frame).unwrap();
        trajectory.write(&frame).unwrap();
        drop(trajectory);

        let mut trajectory = Trajectory::open(named_tmpfile.path()).unwrap();
        check_frame(&trajectory.read().unwrap().unwrap());
        check_frame(&trajectory.read().unwrap().unwrap());
    }
}
