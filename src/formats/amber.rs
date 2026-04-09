// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::{error::CError, frame::Frame, property::Property, unit_cell::UnitCell};
use core::f64;
use log::{debug, warn};
use netcdf3::{DataSet, DataType, FileReader, Variable};
use std::path::Path;

#[derive(Debug)]
struct VariableWithScale {
    var: Variable,
    scale: f64,
}
#[derive(Debug)]
struct Variables {
    coordinates: Option<VariableWithScale>,
    velocities: Option<VariableWithScale>,
    cell_lengths: Option<VariableWithScale>,
    cell_angles: Option<VariableWithScale>,
    time: Option<VariableWithScale>,
}

/// AMBER does not use the text codec traits. It is a dedicated binary
/// format entry point with its own owned state.
#[derive(Debug)]
pub struct AMBERTrajFormat {
    file_reader: FileReader,
    index: usize,
    file_title: String,
    variables: Variables,
    n_atoms: usize,
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
            let buffer = file_reader
                .read_record_f32(variable.var.name(), index)
                .map_err(|e| CError::GenericError(e.to_string()))?;
            for (idx, item) in array.iter_mut().enumerate() {
                item[0] = f64::from(buffer[3 * idx + 0]);
                item[1] = f64::from(buffer[3 * idx + 1]);
                item[2] = f64::from(buffer[3 * idx + 2]);
            }
        }
        DataType::F64 => {
            let buffer = file_reader
                .read_record_f64(variable.var.name(), index)
                .map_err(|e| CError::GenericError(e.to_string()))?;
            for (idx, item) in array.iter_mut().enumerate() {
                item[0] = buffer[3 * idx + 0];
                item[1] = buffer[3 * idx + 1];
                item[2] = buffer[3 * idx + 2];
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

impl AMBERTrajFormat {
    fn not_implemented() -> CError {
        CError::UnsupportedFileFormat("AMBER is not implemented yet".to_string())
    }

    fn read_cell(&mut self) -> Result<Option<UnitCell>, CError> {
        if self.variables.cell_lengths.is_none() || self.variables.cell_angles.is_none() {
            return Ok(None);
        };

        let mut lengths = [0.0; 3];
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
                let buffer = self
                    .file_reader
                    .read_record_f32("cell_lengths", self.index)
                    .map_err(|e| CError::GenericError(e.to_string()))?;

                lengths[0] = f64::from(buffer[0]);
                lengths[1] = f64::from(buffer[1]);
                lengths[2] = f64::from(buffer[2]);
            }
            DataType::F64 => {
                let buffer = self
                    .file_reader
                    .read_record_f64("cell_lengths", self.index)
                    .map_err(|e| CError::GenericError(e.to_string()))?;
                lengths[0] = buffer[0];
                lengths[1] = buffer[1];
                lengths[2] = buffer[2];
            }
        };

        let mut angles = [0.0; 3];
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
                let buffer = self
                    .file_reader
                    .read_record_f32("cell_angles", self.index)
                    .map_err(|e| CError::GenericError(e.to_string()))?;

                angles[0] = f64::from(buffer[0]);
                angles[1] = f64::from(buffer[1]);
                angles[2] = f64::from(buffer[2]);
            }
            DataType::F64 => {
                let buffer = self
                    .file_reader
                    .read_record_f64("cell_angles", self.index)
                    .map_err(|e| CError::GenericError(e.to_string()))?;
                angles[0] = buffer[0];
                angles[1] = buffer[1];
                angles[2] = buffer[2];
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
                scale: 1.0,
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
                scale: 1.0,
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
                scale: 1.0,
            });
            if let Some(units) = read_cell_lengths.get_attr_as_string("units") {
                cell_lengths.as_mut().expect("we just init'ed it").scale *=
                    scale_for_distance(units.as_str());
            } else {
                debug!("not scaling cell lengths");
            };

            cell_angles = Some(VariableWithScale {
                var: read_cell_angles.clone(),
                scale: 1.0,
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
                scale: 1.0,
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
        // dbg!(&variables);
        let format = AMBERTrajFormat {
            file_reader,
            index: 0,
            file_title,
            variables,
            n_atoms,
        };

        // TODO: append mode, line 143 chemfiles/src/formats/AmberNetCDF.cpp

        Ok(format)
    }

    pub fn create(_path: &Path) -> Result<Self, CError> {
        Err(Self::not_implemented())
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

        if let Some(coordinates) = self.variables.coordinates.as_ref() {
            read_array(
                &mut self.file_reader,
                self.index,
                &coordinates,
                frame.positions_mut(),
            )?;
        }

        if let Some(velocities) = self.variables.velocities.as_ref() {
            read_array(
                &mut self.file_reader,
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
                    let value = self
                        .file_reader
                        .read_record_f32("time", 0)
                        .map_err(|e| CError::GenericError(e.to_string()))?;
                    time_value = time.scale * f64::from(value[0]);
                }
                DataType::F64 => {
                    let value = self
                        .file_reader
                        .read_record_f64("time", 0)
                        .map_err(|e| CError::GenericError(e.to_string()))?;
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
        self.file_reader
            .data_set()
            .num_records()
            .ok_or(CError::GenericError(
                "should not have unlimited-size dimension".to_string(),
            ))
    }

    pub fn is_empty(&self) -> Result<bool, CError> {
        self.len().map(|n| n == 0)
    }

    pub fn write(&mut self, _frame: &Frame) -> Result<(), CError> {
        Err(Self::not_implemented())
    }

    pub fn finish(&mut self) -> Result<(), CError> {
        Err(Self::not_implemented())
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use assert_approx_eq::assert_approx_eq;

    use crate::trajectory::Trajectory;

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
}
