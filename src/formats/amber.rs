// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::{
    atom::Atom, bond::BondOrder, error::CError, frame::Frame, property::Property, residue::Residue,
    topology::Topology, unit_cell::UnitCell,
};
use core::f64;
use log::{debug, warn};
use netcdf3::{DataSet, DataType, FileReader, FileWriter, Variable, Version};
use std::{
    fs::{File, OpenOptions},
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    str::FromStr,
};

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
    convention: Convention,
    file_path: Option<PathBuf>,
    // Writing state
    write_path: Option<PathBuf>,
    buffered_frames: Vec<BufferedFrame>,
    has_velocities: bool,
    initialized: bool,
    mode: FileMode,
    version: Option<Version>,

    topology: Option<Topology>,
    topology_checked: bool,
}

// TODO: probably this needs to go somewhere more general so it can be reused with other binary formats
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub(crate) enum FileMode {
    #[default]
    Read,
    Append,
}

/// Which AMBER `NetCDF` flavour the file uses.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub(crate) enum Convention {
    /// AMBER trajectory (multi-frame, unlimited `frame` dimension).
    #[default]
    Amber,
    /// AMBER restart (single snapshot, no `frame` dimension).
    Restart,
}

impl Convention {
    fn as_str(self) -> &'static str {
        match self {
            Convention::Amber => "AMBER",
            Convention::Restart => "AMBERRESTART",
        }
    }
}

enum NumericData {
    F32(Vec<f32>),
    F64(Vec<f64>),
}

fn read_values(
    file_reader: &mut FileReader,
    index: usize,
    convention: Convention,
    variable: &Variable,
    invalid_type_message: &str,
) -> Result<NumericData, CError> {
    let name = variable.name();
    match variable.data_type() {
        DataType::I8 | DataType::U8 | DataType::I16 | DataType::I32 => {
            Err(CError::GenericError(invalid_type_message.to_string()))
        }
        DataType::F32 => {
            let values = match convention {
                Convention::Amber => file_reader.read_record_f32(name, index)?,
                Convention::Restart => file_reader.read_var_f32(name)?,
            };
            Ok(NumericData::F32(values))
        }
        DataType::F64 => {
            let values = match convention {
                Convention::Amber => file_reader.read_record_f64(name, index)?,
                Convention::Restart => file_reader.read_var_f64(name)?,
            };
            Ok(NumericData::F64(values))
        }
    }
}

fn read_triplet(
    file_reader: &mut FileReader,
    index: usize,
    convention: Convention,
    variable: &VariableWithScale,
) -> Result<[f64; 3], CError> {
    match read_values(
        file_reader,
        index,
        convention,
        &variable.var,
        "invalid type for variable, expected f32/f64 data",
    )? {
        NumericData::F32(values) => Ok([
            variable.scale * f64::from(values[0]),
            variable.scale * f64::from(values[1]),
            variable.scale * f64::from(values[2]),
        ]),
        NumericData::F64(values) => Ok([
            variable.scale * values[0],
            variable.scale * values[1],
            variable.scale * values[2],
        ]),
    }
}

fn read_scalar(
    file_reader: &mut FileReader,
    index: usize,
    convention: Convention,
    variable: &VariableWithScale,
    invalid_type_message: &str,
) -> Result<f64, CError> {
    match read_values(
        file_reader,
        index,
        convention,
        &variable.var,
        invalid_type_message,
    )? {
        NumericData::F32(values) => Ok(variable.scale * f64::from(values[0])),
        NumericData::F64(values) => Ok(variable.scale * values[0]),
    }
}

fn read_array(
    file_reader: &mut FileReader,
    index: usize,
    convention: Convention,
    variable: &VariableWithScale,
    array: &mut [[f64; 3]],
) -> Result<(), CError> {
    match read_values(
        file_reader,
        index,
        convention,
        &variable.var,
        "invalid type for variable, expected f32/f64 data",
    )? {
        NumericData::F32(values) => {
            for (chunk, item) in values.chunks_exact(3).zip(array.iter_mut()) {
                item[0] = variable.scale * f64::from(chunk[0]);
                item[1] = variable.scale * f64::from(chunk[1]);
                item[2] = variable.scale * f64::from(chunk[2]);
            }
        }
        NumericData::F64(values) => {
            for (chunk, item) in values.chunks_exact(3).zip(array.iter_mut()) {
                item[0] = variable.scale * chunk[0];
                item[1] = variable.scale * chunk[1];
                item[2] = variable.scale * chunk[2];
            }
        }
    }

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
            if read_attr != expected_attr {
                return Err(CError::GenericError(format!(
                    "expected '{expected_attr}', got {read_attr}"
                )));
            }
        } else {
            return Err(CError::GenericError(
                "could not read attr as string".to_string(),
            ));
        }

        Ok(())
    }
    check_attr(data_set, "Conventions", convention)?;
    check_attr(data_set, "ConventionVersion", "1.0")?;

    if let Some(spatial) = file_reader.data_set().get_dim("spatial") {
        if spatial.size() != 3 {
            return Err(CError::GenericError(format!(
                "'spatial' dimension must have a size of 3, got {}",
                spatial.size()
            )));
        }
    } else {
        return Err(CError::GenericError(
            "missing 'spatial' dimension".to_string(),
        ));
    }

    file_reader
        .data_set()
        .get_dim("atom")
        .ok_or(CError::GenericError("missing 'atom' dimension".to_string()))?;

    if let Some(cell_spatial) = file_reader.data_set().get_dim("cell_spatial")
        && cell_spatial.size() != 3
    {
        return Err(CError::GenericError(format!(
            "'cell_spatial' dimension must have a size of 3, got {}",
            cell_spatial.size()
        )));
    }
    if let Some(cell_angular) = file_reader.data_set().get_dim("cell_angular")
        && cell_angular.size() != 3
    {
        return Err(CError::GenericError(format!(
            "'cell_angular' dimension must have a size of 3, got {}",
            cell_angular.size()
        )));
    }

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
    if let Some(scale) = splitted.next()
        && let Some(time_unit) = splitted.next().as_mut()
    {
        let scale = scale_for_distance(scale);
        let time_unit = scale_for_time(time_unit);
        return scale / time_unit;
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
        "microseconds" | "microsecond" | "µs" | "us" => 1e6,
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

#[derive(Copy, Debug, Clone)]
#[allow(dead_code)]
enum FortranFormat {
    Char { width: u16, repeat: u16 },
    Int { width: u16, repeat: u16 },
    Real { width: u16, repeat: u16 },
}

#[derive(Copy, Clone, Debug)]
enum Section {
    Title,

    /// Information about how many parameters are present in all of the sections.
    Pointers,

    /// The atom name for every atom.
    AtomName,

    /// The charge for every atom. There are `NATOM` floating point numbers in this section.
    Charge,

    /// The atomic number of every atom. There are `NATOM` integers in this section
    AtomicNumber,
    Mass,
    AtomTypeIndex,
    NumberExcludedAtoms,
    NonbondedParmIndex,

    /// Contains the residue name for every residue. There are `NRES` 4-character strings in this section
    ResidueLabel,

    /// Lists the first atom in each residue. There are `ǸRES` integers in this section.
    ResiduePointer,

    /// Lists all of the bond force constants (k in Eq. 2) in units kcal mol−1 Å^{−2} for each unique bond type.
    BondForceConstant,
    BondEquilValue,
    AngleForceConstant,
    AngleEquilValue,
    DihedralForceConstant,
    DihedralPeriodicity,
    DihedralPhase,
    SceeScaleFactor,
    ScnbScaleFactor,

    /// This section is currently unused, and while ‘future use’ is planned, this assertion has lain dormant for some time.
    Solty,

    /// the LJ A-coefficients (ai,j in Eq. 5) for all pairs of distinct LJ types.
    LennardJonesAcoef,

    /// the LJ B-coefficients (bi,j in Eq. 5) for all pairs of distinct LJ types.
    LennardJonesBcoef,

    /// list of every bond in the system in which at least one atom is Hydrogen.
    /// Each bond is identified by 3 integers—the two
    /// atoms involved in the bond and the index into the `BOND FORCE CONSTANT` and `BOND EQUIL VALUE`
    /// There are 3 x `NBONH` integers in this section
    BondsIncHydrogen,

    /// There are 3 x `NBONA` integers in this section
    BondsWithoutHydrogen,
    AnglesIncHydrogen,
    AnglesWithoutHydrogen,
    DihedralsIncHydrogen,
    DihedralsWithoutHydrogen,
    ExcludedAtomsList,
    HbondAcoef,
    HbondBcoef,
    HBCut,
    AmberAtomType,
    TreeChainClassification,
    JoinArray,
    Irotat,
    SolventPointers,
    AtomsPerMolecule,
    BoxDimensions,
    CapInfo,
    CapInfo2,
    RadiusSet,
    Radii,
    Ipol,
    Screen,
    Polarizability,
}

impl TryFrom<&str> for Section {
    type Error = CError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "TITLE" => Ok(Section::Title),
            "POINTERS" => Ok(Section::Pointers),
            "ATOM_NAME" => Ok(Section::AtomName),
            "CHARGE" => Ok(Section::Charge),
            "ATOMIC_NUMBER" => Ok(Section::AtomicNumber),
            "MASS" => Ok(Section::Mass),
            "ATOM_TYPE_INDEX" => Ok(Section::AtomTypeIndex),
            "NUMBER_EXCLUDED_ATOMS" => Ok(Section::NumberExcludedAtoms),
            "NONBONDED_PARM_INDEX" => Ok(Section::NonbondedParmIndex),
            "RESIDUE_LABEL" => Ok(Section::ResidueLabel),
            "RESIDUE_POINTER" => Ok(Section::ResiduePointer),
            "BOND_FORCE_CONSTANT" => Ok(Section::BondForceConstant),
            "BOND_EQUIL_VALUE" => Ok(Section::BondEquilValue),
            "ANGLE_FORCE_CONSTANT" => Ok(Section::AngleForceConstant),
            "ANGLE_EQUIL_VALUE" => Ok(Section::AngleEquilValue),
            "DIHEDRAL_FORCE_CONSTANT" => Ok(Section::DihedralForceConstant),
            "DIHEDRAL_PERIODICITY" => Ok(Section::DihedralPeriodicity),
            "DIHEDRAL_PHASE" => Ok(Section::DihedralPhase),
            "SCEE_SCALE_FACTOR" => Ok(Section::SceeScaleFactor),
            "SCNB_SCALE_FACTOR" => Ok(Section::ScnbScaleFactor),
            "SOLTY" => Ok(Section::Solty),
            "LENNARD_JONES_ACOEF" => Ok(Section::LennardJonesAcoef),
            "LENNARD_JONES_BCOEF" => Ok(Section::LennardJonesBcoef),
            "BONDS_INC_HYDROGEN" => Ok(Section::BondsIncHydrogen),
            "BONDS_WITHOUT_HYDROGEN" => Ok(Section::BondsWithoutHydrogen),
            "ANGLES_INC_HYDROGEN" => Ok(Section::AnglesIncHydrogen),
            "ANGLES_WITHOUT_HYDROGEN" => Ok(Section::AnglesWithoutHydrogen),
            "DIHEDRALS_INC_HYDROGEN" => Ok(Section::DihedralsIncHydrogen),
            "DIHEDRALS_WITHOUT_HYDROGEN" => Ok(Section::DihedralsWithoutHydrogen),
            "EXCLUDED_ATOMS_LIST" => Ok(Section::ExcludedAtomsList),
            "HBOND_ACOEF" => Ok(Section::HbondAcoef),
            "HBOND_BCOEF" => Ok(Section::HbondBcoef),
            "HBCUT" => Ok(Section::HBCut),
            "AMBER_ATOM_TYPE" => Ok(Section::AmberAtomType),
            "TREE_CHAIN_CLASSIFICATION" => Ok(Section::TreeChainClassification),
            "JOIN_ARRAY" => Ok(Section::JoinArray),
            "IROTAT" => Ok(Section::Irotat),
            "SOLVENT_POINTERS" => Ok(Section::SolventPointers),
            "ATOMS_PER_MOLECULE" => Ok(Section::AtomsPerMolecule),
            "BOX_DIMENSIONS" => Ok(Section::BoxDimensions),
            "CAP_INFO" => Ok(Section::CapInfo),
            "CAP_INFO2" => Ok(Section::CapInfo2),
            "RADIUS_SET" => Ok(Section::RadiusSet),
            "RADII" => Ok(Section::Radii),
            "IPOL" => Ok(Section::Ipol),
            "SCREEN" => Ok(Section::Screen),
            "POLARIZABILITY" => Ok(Section::Polarizability),
            unknown_section => Err(CError::GenericError(format!(
                "this section has not yet been implemented: {unknown_section}"
            ))),
        }
    }
}

#[derive(Copy, Clone, Debug, Default)]
enum PeriodicBox {
    #[default]
    None,
    Orthorhombic,
    TruncatedOctahedron,
}

#[derive(Copy, Clone, Debug, Default)]
/// The following descriptions are directly from the specification: https://ambermd.org/prmtop.pdf
#[allow(non_snake_case, dead_code)]
struct Pointers {
    /// Number of atoms
    NATOM: u64,

    /// Number of distinct Lennard-Jones atom types
    NTYPES: u64,

    /// Number of bonds containing Hydrogen
    NBONH: u64,

    /// Number of bonds not containing Hydrogen
    MBONA: u64,

    ///  Number of angles containing Hydrogen
    NTHETH: u64,

    /// Number of angles not containing Hydrogen
    MTHETA: u64,

    /// Number of torsions containing Hydrogen
    NPHIH: u64,

    /// Number of torsions not containing Hydrogen
    MPHIA: u64,

    /// Not currently used for anything
    NHPARM: u64,

    /// Used to determine if this is a LES-compatible prmtop
    NPARM: u64,

    /// Number of excluded atoms (length of total exclusion list)
    NNB: u64,

    /// Number of residues
    NRES: u64,

    /// MBONA + number of constraint bonds. AMBER codes no longer support constraints in the topology file
    NBONA: u64,

    /// MTHETA + number of constraint angles. AMBER codes no longer support constraints in the topology file
    NTHETA: u64,

    /// MPHIA + number of constraint torsions. AMBER codes no longer support constraints in the topology file
    NPHIA: u64,

    /// Number of unique bond types.
    NUMBND: u64,

    /// Number of unique angle types.
    NUMANG: u64,

    /// Number of unique torsion types.
    NPTRA: u64,

    /// Number of SOLTY terms. Currently unused.
    NATYP: u64,

    /// Number of distinct 10-12 hydrogen bond pair types. Modern AMBER force fields do not use a 10-12 potential
    NPHB: u64,

    /// Set to 1 if topology contains residue perturbation information. No AMBER codes support perturbed topologies anymore
    IFPERT: u64,

    /// Number of perturbed bonds. No AMBER codes support perturbed topologies anymore
    NBPER: u64,

    /// Number of perturbed angles. No AMBER codes support perturbed topologies anymore
    NGPER: u64,

    /// Number of perturbed torsions. No AMBER codes support perturbed topologies anymore
    NDPER: u64,

    /// Number of bonds in which both atoms are being perturbed. No AMBER codes support perturbed topologies anymore
    MBPER: u64,

    /// Number of angles in which all 3 atoms are being perturbed. No AMBER codes support perturbed topologies anymore
    MGPER: u64,

    /// Number of torsions in which all 4 atoms are being perturbed. No AMBER codes support perturbed topologies anymore
    MDPER: u64,

    /// Flag indicating whether a periodic box is present
    IFBOX: PeriodicBox,

    /// Number of atoms in the largest residue
    NMXRS: u64,

    /// Set to `true` if a solvent CAP is being used
    IFCAP: bool,

    /// Number of extra points in the topology file
    NUMEXTRA: u64,

    /// Number of PIMD slices or number of beads. It might not be present.
    NCOPY: Option<u64>,
}

fn symbol_from_atomic_number(n: i64) -> &'static str {
    // TODO: maybe use mendeleev?
    const SYMBOLS: &[&str] = &[
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", // 0-10
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", // 11-20
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", // 21-30
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", // 31-40
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", // 41-50
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", // 51-60
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", // 61-70
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", // 71-80
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", // 81-90
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", // 91-100
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", // 101-110
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", // 111-118
    ];
    if n >= 0 {
        SYMBOLS.get(n as usize).copied().unwrap_or("")
    } else {
        ""
    }
}

fn symbol_from_mass(mass: f64) -> &'static str {
    // Standard atomic weights (IUPAC 2021); find the element with the closest mass.
    // TODO: maybe use mendeleev?
    const ENTRIES: &[(f64, &str)] = &[
        (1.008, "H"),
        (4.003, "He"),
        (6.941, "Li"),
        (9.012, "Be"),
        (10.811, "B"),
        (12.011, "C"),
        (14.007, "N"),
        (15.999, "O"),
        (18.998, "F"),
        (20.180, "Ne"),
        (22.990, "Na"),
        (24.305, "Mg"),
        (26.982, "Al"),
        (28.086, "Si"),
        (30.974, "P"),
        (32.060, "S"),
        (35.450, "Cl"),
        (39.948, "Ar"),
        (39.098, "K"),
        (40.078, "Ca"),
        (44.956, "Sc"),
        (47.867, "Ti"),
        (50.942, "V"),
        (51.996, "Cr"),
        (54.938, "Mn"),
        (55.845, "Fe"),
        (58.933, "Co"),
        (58.693, "Ni"),
        (63.546, "Cu"),
        (65.380, "Zn"),
        (69.723, "Ga"),
        (72.630, "Ge"),
        (74.922, "As"),
        (78.971, "Se"),
        (79.904, "Br"),
        (83.798, "Kr"),
        (85.468, "Rb"),
        (87.620, "Sr"),
        (88.906, "Y"),
        (91.224, "Zr"),
        (92.906, "Nb"),
        (95.950, "Mo"),
        (98.000, "Tc"),
        (101.070, "Ru"),
        (102.906, "Rh"),
        (106.420, "Pd"),
        (107.868, "Ag"),
        (112.414, "Cd"),
        (114.818, "In"),
        (118.710, "Sn"),
        (121.760, "Sb"),
        (127.600, "Te"),
        (126.904, "I"),
        (131.293, "Xe"),
        (132.905, "Cs"),
        (137.327, "Ba"),
        (138.905, "La"),
        (140.116, "Ce"),
        (140.908, "Pr"),
        (144.242, "Nd"),
        (145.000, "Pm"),
        (150.360, "Sm"),
        (151.964, "Eu"),
        (157.250, "Gd"),
        (158.925, "Tb"),
        (162.500, "Dy"),
        (164.930, "Ho"),
        (167.259, "Er"),
        (168.934, "Tm"),
        (173.054, "Yb"),
        (174.967, "Lu"),
        (178.490, "Hf"),
        (180.948, "Ta"),
        (183.840, "W"),
        (186.207, "Re"),
        (190.230, "Os"),
        (192.217, "Ir"),
        (195.084, "Pt"),
        (196.967, "Au"),
        (200.592, "Hg"),
        (204.383, "Tl"),
        (207.200, "Pb"),
        (208.980, "Bi"),
        (209.000, "Po"),
        (210.000, "At"),
        (222.000, "Rn"),
        (223.000, "Fr"),
        (226.000, "Ra"),
        (227.000, "Ac"),
        (232.038, "Th"),
        (231.036, "Pa"),
        (238.029, "U"),
        (237.000, "Np"),
        (244.000, "Pu"),
        (243.000, "Am"),
        (247.000, "Cm"),
        (247.000, "Bk"),
        (251.000, "Cf"),
        (252.000, "Es"),
        (257.000, "Fm"),
        (258.000, "Md"),
        (259.000, "No"),
        (262.000, "Lr"),
        (267.000, "Rf"),
        (268.000, "Db"),
        (271.000, "Sg"),
        (272.000, "Bh"),
        (270.000, "Hs"),
        (276.000, "Mt"),
        (281.000, "Ds"),
        (280.000, "Rg"),
        (285.000, "Cn"),
        (284.000, "Nh"),
        (289.000, "Fl"),
        (288.000, "Mc"),
        (293.000, "Lv"),
        (292.000, "Ts"),
        (294.000, "Og"),
    ];
    const THRESHOLD: f64 = 0.2;
    ENTRIES
        .iter()
        .min_by(|a, b| (a.0 - mass).abs().partial_cmp(&(b.0 - mass).abs()).unwrap())
        .and_then(|(m, sym)| {
            let diff = (m - mass).abs();
            if diff > THRESHOLD { None } else { Some(*sym) }
        })
        .map_or("", |sym| sym)
}

impl AMBERTrajFormat {
    fn read_cell(&mut self) -> Result<Option<UnitCell>, CError> {
        if self.variables.cell_lengths.is_none() || self.variables.cell_angles.is_none() {
            return Ok(None);
        }
        let convention = self.convention;
        let index = self.index;
        let reader = self
            .reader
            .as_mut()
            .expect("reader should've been initialized");

        let cell_lengths = self
            .variables
            .cell_lengths
            .as_ref()
            .expect("we just checked for None");
        let lengths = read_triplet(reader, index, convention, cell_lengths)?;

        let cell_angles = self
            .variables
            .cell_angles
            .as_ref()
            .expect("we just checked for None");
        let angles = read_triplet(reader, index, convention, cell_angles)?;

        Ok(Some(UnitCell::new_from_lengths_angles(
            lengths.into(),
            angles.into(),
        )?))
    }

    pub(crate) fn open(
        path: &Path,
        mode: FileMode,
        convention: Convention,
    ) -> Result<Self, CError> {
        // For append mode, if the file is missing or not yet a valid NetCDF-3 file,
        // treat it the same as creating a new file.
        let file_reader = match FileReader::open(path) {
            Ok(r) => r,
            Err(_) if mode == FileMode::Append => {
                return Ok(Self {
                    write_path: Some(path.to_path_buf()),
                    mode,
                    convention,
                    ..Default::default()
                });
            }
            Err(e) => return Err(CError::GenericError(e.to_string())),
        };

        validate_common(&file_reader, convention.as_str())?;

        let file_title = file_reader
            .data_set()
            .get_global_attr_as_string("title")
            .unwrap_or_default();

        let n_atoms = file_reader
            .data_set()
            .get_dim("atom")
            .expect("we just just validated that 'atom' exists")
            .size();

        // Get the variables actually defined in the file
        let coordinates = if let Some(coordinates) = file_reader.data_set().get_var("coordinates") {
            let mut coords = Some(VariableWithScale {
                var: coordinates.clone(),
                scale: scale_factor(coordinates)?,
            });
            if let Some(units) = coordinates.get_attr_as_string("units") {
                coords.as_mut().expect("we just init'ed").scale *=
                    scale_for_distance(units.as_str());
            } else {
                debug!("not scaling coordinates");
            }
            coords
        } else {
            warn!("the coordinates variable is not defined in this file.");
            None
        };

        let velocities = if let Some(velocities) = file_reader.data_set().get_var("velocities") {
            let mut velo = Some(VariableWithScale {
                var: velocities.clone(),
                scale: scale_factor(velocities)?,
            });
            if let Some(units) = velocities.get_attr_as_string("units") {
                velo.as_mut().expect("we just init'ed").scale *= scale_for_velocity(units.as_str());
            } else {
                debug!("not scaling velocities");
            }
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
                scale: scale_factor(read_cell_lengths)?,
            });
            if let Some(units) = read_cell_lengths.get_attr_as_string("units") {
                cell_lengths.as_mut().expect("we just init'ed it").scale *=
                    scale_for_distance(units.as_str());
            } else {
                debug!("not scaling cell lengths");
            }

            cell_angles = Some(VariableWithScale {
                var: read_cell_angles.clone(),
                scale: scale_factor(read_cell_angles)?,
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
            }
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
        }

        let time = if let Some(time) = file_reader.data_set().get_var("time") {
            let mut t = Some(VariableWithScale {
                var: time.clone(),
                scale: scale_factor(time)?,
            });
            if let Some(units) = time.get_attr_as_string("units") {
                t.as_mut().expect("we just init'ed").scale *= scale_for_time(units.as_str());
            } else {
                debug!("not scaling time");
            }
            t
        } else {
            warn!("the timevariable is not defined in this file.");
            None
        };

        let variables = Variables {
            coordinates,
            velocities,
            cell_lengths,
            cell_angles,
            time,
        };
        let index = match convention {
            Convention::Amber if mode == FileMode::Append => file_reader
                .data_set()
                .num_records()
                .expect("data_set to have records"),
            Convention::Restart | Convention::Amber => 0,
        };
        let version = Some(file_reader.version());
        let has_velocities = variables.velocities.is_some();
        let write_path = if mode == FileMode::Append {
            Some(path.to_path_buf())
        } else {
            None
        };
        let format = AMBERTrajFormat {
            reader: Some(file_reader),
            index,
            file_title,
            variables,
            n_atoms,
            mode,
            convention,
            file_path: Some(path.to_path_buf()),
            version,
            has_velocities,
            write_path,
            ..Default::default()
        };

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

    fn topology(&mut self) -> Result<Option<&Topology>, CError> {
        if !self.topology_checked {
            self.topology = self.read_topology_file()?;
            self.topology_checked = true;
        }

        Ok(self.topology.as_ref())
    }

    pub fn read_at(&mut self, index: usize) -> Result<Frame, CError> {
        self.index = index;

        let mut frame = Frame::new();
        if let Some(unitcell) = self.read_cell()? {
            frame.set_unit_cell(unitcell);
        }

        if !self.file_title.is_empty() {
            frame.properties.insert(
                "name".to_string(),
                crate::property::Property::String(self.file_title.clone()),
            );
        }

        frame.resize(self.n_atoms)?;

        let convention = self.convention;
        let index = self.index;
        let reader = self
            .reader
            .as_mut()
            .expect("reader should've been initialized");
        if let Some(coordinates) = self.variables.coordinates.as_ref() {
            read_array(
                reader,
                index,
                convention,
                coordinates,
                frame.positions_mut(),
            )?;
        }

        if let Some(velocities) = self.variables.velocities.as_ref() {
            frame.add_velocities();
            read_array(
                reader,
                index,
                convention,
                velocities,
                frame.velocities_mut().expect("just resized velocities"),
            )?;
        }

        if let Some(time) = self.variables.time.as_ref() {
            let time_value = read_scalar(
                reader,
                index,
                convention,
                time,
                "invalid type for time variable",
            )?;
            frame
                .properties
                .insert("time".to_string(), Property::Double(time_value));
        }

        if let Some(topology) = self.topology()? {
            frame.set_topology(topology.clone())?;
        }

        Ok(frame)
    }

    fn read_format(substring: &str) -> Option<FortranFormat> {
        let bytes = substring.as_bytes();
        let mut i = 0;

        let mut repeat: u16 = 0;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            repeat = repeat
                .checked_mul(10)?
                .checked_add(u16::from(bytes[i] - b'0'))?;
            i += 1;
        }
        if i == 0 || i >= bytes.len() {
            return None;
        }

        let kind = bytes[i] as char;
        i += 1;

        let start = i;
        let mut width: u16 = 0;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            width = width
                .checked_mul(10)?
                .checked_add(u16::from(bytes[i] - b'0'))?;
            i += 1;
        }

        if kind == 'e' || kind == 'E' {
            // TODO: handle the rest of the float here?
        } else if start == i || i != bytes.len() {
            return None;
        }

        match kind {
            'a' | 'A' => Some(FortranFormat::Char { width, repeat }),
            'i' | 'I' => Some(FortranFormat::Int { width, repeat }),
            'e' | 'E' => Some(FortranFormat::Real { width, repeat }),
            _newkind => None,
        }
    }

    fn read_section<T: FromStr>(
        line_buffer: &mut String,
        reader: &mut BufReader<File>,
        vec_size: usize,
        width: usize,
    ) -> Result<(bool, bool, Vec<T>), CError> {
        let mut data_vec: Vec<T> = Vec::with_capacity(vec_size);

        while reader
            .read_line(line_buffer)
            .map_err(|e| CError::GenericError(e.to_string()))?
            > 0
        {
            if line_buffer.contains("%FLAG") {
                return Ok((false, true, data_vec));
            }

            let trimmed = line_buffer.trim_end_matches(['\n', '\r']);
            for i in (0..trimmed.len()).step_by(width) {
                if i + width > trimmed.len() {
                    break;
                }
                data_vec.push(
                    trimmed[i..i + width]
                        .trim_ascii()
                        .parse::<T>()
                        .map_err(|_| CError::GenericError("failed to parse field".to_string()))?,
                );
            }
            line_buffer.clear();
        }

        if line_buffer.contains("%FLAG") {
            Ok((false, true, data_vec))
        } else {
            Ok((true, false, data_vec))
        }
    }

    fn read_topology_file(&self) -> Result<Option<Topology>, CError> {
        // This implementation is based on the documentation at https://ambermd.org/prmtop.pdf
        let mut topology = Topology::default();
        topology.reserve(self.n_atoms);

        let prmtop_path = match self.file_path.as_ref().map(|p| p.with_extension("parm7")) {
            Some(p) if p.exists() => p,
            _ => return Ok(None),
        };
        let mut reader = BufReader::with_capacity(64 * 1024, File::open(&prmtop_path)?);

        let mut line = String::new();
        // The first line is the version string
        reader
            .read_line(&mut line)
            .map_err(|e| CError::GenericError(e.to_string()))?;
        if !line.starts_with("%VERSION") {
            return Err(CError::GenericError(
                "expected '%VERSION' header in parm7 file".to_string(),
            ));
        }
        line.clear();
        let mut end = false;
        // If the topology file has a %FLAG TITLE then it is an AMBER topology.
        // If it has a %FLAG CTITLE instead, then it is a chamber topology
        let mut already_read_flag_header = false;
        let mut residue_labels = None;
        let mut residue_pointers = None;
        let mut pointers = Pointers::default();
        let mut atom_names = None;
        let mut charges = None;
        let mut atomic_numbers = None;
        let mut bonds_h = None;
        let mut bonds_noh = None;
        let mut masses = None;
        while !end {
            // First line is the type of section
            if !already_read_flag_header {
                let bytes = reader
                    .read_line(&mut line)
                    .map_err(|e| CError::GenericError(e.to_string()))?;
                if bytes == 0 {
                    break;
                }
            }
            if !line.starts_with("%FLAG ") {
                return Err(CError::GenericError(format!(
                    "expected '%FLAG' line, got: {line}"
                )));
            }
            let section_name = line[6..].trim_ascii();
            let Ok(section) = section_name.try_into() else {
                debug!("skipping unknown parm7 section: {section_name}");
                end = Self::skip_section(&mut line, &mut reader, &mut already_read_flag_header);
                continue;
            };
            line.clear();

            // Second line is the format of the section.
            reader
                .read_line(&mut line)
                .map_err(|e| CError::GenericError(e.to_string()))?;
            let fortran_format = if let Some(fmt_str) = line.strip_prefix("%FORMAT(") {
                let fmt_str = &fmt_str[..fmt_str.trim_ascii_end().len() - 1];
                if let Some(f) = Self::read_format(fmt_str) {
                    f
                } else {
                    debug!("unrecognised Fortran format '{fmt_str}', skipping section");
                    end = Self::skip_section(&mut line, &mut reader, &mut already_read_flag_header);
                    continue;
                }
            } else {
                return Err(CError::GenericError(format!(
                    "expected '%FORMAT(...)' line, got: {line}"
                )));
            };

            // now the data of the section starts
            line.clear();
            match section {
                // for now, we skip the title
                Section::Title => {
                    reader
                        .read_line(&mut line)
                        .map_err(|e| CError::GenericError(e.to_string()))?;
                    line.clear();
                }
                Section::Pointers => {
                    let mut values = Vec::with_capacity(32);

                    if let FortranFormat::Int { width, repeat } = fortran_format {
                        let width = width as usize;
                        let repeat = repeat as usize;

                        for _ in 0..3 {
                            reader
                                .read_line(&mut line)
                                .map_err(|e| CError::GenericError(e.to_string()))?;
                            let trimmed_line = line.trim_end();
                            for i in 0..repeat {
                                values.push(
                                    trimmed_line[i * width..i * width + width]
                                        .trim_ascii()
                                        .parse::<u64>()
                                        .map_err(|_| {
                                            CError::GenericError(
                                                "failed to parse POINTERS value".to_string(),
                                            )
                                        })?,
                                );
                            }
                            line.clear();
                        }
                        // special-case the last line of pointers: it might hold one or two numbers
                        reader
                            .read_line(&mut line)
                            .map_err(|e| CError::GenericError(e.to_string()))?;
                        let trimmed_line = line.trim_end();
                        for i in 0..2 {
                            if let Some(value) = trimmed_line
                                .get(i * width..i * width + width)
                                .and_then(|s| s.trim_ascii().parse::<u64>().ok())
                            {
                                values.push(value);
                            }
                        }
                    }
                    line.clear();

                    let ifbox = match values[27] {
                        0 => PeriodicBox::None,
                        1 => PeriodicBox::Orthorhombic,
                        2 => PeriodicBox::TruncatedOctahedron,
                        _ => PeriodicBox::None,
                    };

                    pointers = Pointers {
                        NATOM: values[0],
                        NTYPES: values[1],
                        NBONH: values[2],
                        MBONA: values[3],
                        NTHETH: values[4],
                        MTHETA: values[5],
                        NPHIH: values[6],
                        MPHIA: values[7],
                        NHPARM: values[8],
                        NPARM: values[9],
                        NNB: values[10],
                        NRES: values[11],
                        NBONA: values[12],
                        NTHETA: values[13],
                        NPHIA: values[14],
                        NUMBND: values[15],
                        NUMANG: values[16],
                        NPTRA: values[17],
                        NATYP: values[18],
                        NPHB: values[19],
                        IFPERT: values[20],
                        NBPER: values[21],
                        NGPER: values[22],
                        NDPER: values[23],
                        MBPER: values[24],
                        MGPER: values[25],
                        MDPER: values[26],
                        IFBOX: ifbox,
                        NMXRS: values[28],
                        IFCAP: values[29] != 0,
                        NUMEXTRA: values[30],
                        NCOPY: values.get(31).copied(),
                    };
                }
                Section::AtomName => {
                    if let FortranFormat::Char { width, .. } = fortran_format {
                        let result = Self::read_section::<String>(
                            &mut line,
                            &mut reader,
                            pointers.NATOM as usize,
                            width as usize,
                        )?;
                        end = result.0;
                        already_read_flag_header = result.1;
                        atom_names = Some(result.2);
                    }
                }
                Section::Charge => {
                    if let FortranFormat::Real { width, .. } = fortran_format {
                        let result = Self::read_section::<f64>(
                            &mut line,
                            &mut reader,
                            pointers.NATOM as usize,
                            width as usize,
                        )?;
                        end = result.0;
                        already_read_flag_header = result.1;
                        charges = Some(result.2);
                    }
                }
                Section::AtomicNumber => {
                    if let FortranFormat::Int { width, .. } = fortran_format {
                        let result = Self::read_section::<i64>(
                            &mut line,
                            &mut reader,
                            pointers.NATOM as usize,
                            width as usize,
                        )?;
                        end = result.0;
                        already_read_flag_header = result.1;
                        atomic_numbers = Some(result.2);
                    }
                }
                Section::Mass => {
                    if let FortranFormat::Real { width, .. } = fortran_format {
                        let result = Self::read_section::<f64>(
                            &mut line,
                            &mut reader,
                            pointers.NATOM as usize,
                            width as usize,
                        )?;
                        end = result.0;
                        already_read_flag_header = result.1;
                        masses = Some(result.2);
                    }
                }
                Section::ResidueLabel => {
                    if let FortranFormat::Char { width, .. } = fortran_format {
                        let result = Self::read_section::<String>(
                            &mut line,
                            &mut reader,
                            pointers.NRES as usize,
                            width as usize,
                        )?;
                        end = result.0;
                        already_read_flag_header = result.1;
                        residue_labels = Some(result.2);
                    }
                }
                Section::ResiduePointer => {
                    if let FortranFormat::Int { width, .. } = fortran_format {
                        let result = Self::read_section::<u64>(
                            &mut line,
                            &mut reader,
                            pointers.NRES as usize,
                            width as usize,
                        )?;
                        end = result.0;
                        already_read_flag_header = result.1;
                        residue_pointers = Some(result.2);
                    }
                }
                Section::BondsIncHydrogen => {
                    if let FortranFormat::Int { width, .. } = fortran_format {
                        let result = Self::read_section::<usize>(
                            &mut line,
                            &mut reader,
                            pointers.NBONH as usize * 3,
                            width as usize,
                        )?;
                        end = result.0;
                        already_read_flag_header = result.1;
                        bonds_h = Some(result.2);
                    }
                }
                Section::BondsWithoutHydrogen => {
                    if let FortranFormat::Int { width, .. } = fortran_format {
                        let result = Self::read_section::<usize>(
                            &mut line,
                            &mut reader,
                            pointers.MBONA as usize * 3,
                            width as usize,
                        )?;
                        end = result.0;
                        already_read_flag_header = result.1;
                        bonds_noh = Some(result.2);
                    }
                }
                Section::AtomTypeIndex
                | Section::NumberExcludedAtoms
                | Section::NonbondedParmIndex
                | Section::BondForceConstant
                | Section::BondEquilValue
                | Section::AngleForceConstant
                | Section::AngleEquilValue
                | Section::DihedralForceConstant
                | Section::DihedralPeriodicity
                | Section::DihedralPhase
                | Section::SceeScaleFactor
                | Section::ScnbScaleFactor
                | Section::Solty
                | Section::LennardJonesAcoef
                | Section::LennardJonesBcoef
                | Section::AnglesIncHydrogen
                | Section::AnglesWithoutHydrogen
                | Section::DihedralsIncHydrogen
                | Section::DihedralsWithoutHydrogen
                | Section::ExcludedAtomsList
                | Section::HbondAcoef
                | Section::HbondBcoef
                | Section::HBCut
                | Section::AmberAtomType
                | Section::TreeChainClassification
                | Section::JoinArray
                | Section::Irotat
                | Section::SolventPointers
                | Section::AtomsPerMolecule
                | Section::BoxDimensions
                | Section::CapInfo
                | Section::CapInfo2
                | Section::RadiusSet
                | Section::Radii
                | Section::Ipol
                | Section::Screen
                | Section::Polarizability => {
                    end = Self::skip_section(&mut line, &mut reader, &mut already_read_flag_header);
                }
            }
        }

        for i in 0..pointers.NATOM as usize {
            let name = atom_names
                .as_ref()
                .and_then(|v| v.get(i))
                .map(|s| s.trim().to_string())
                .unwrap_or_default();

            let symbol = if let Some(ref nums) = atomic_numbers {
                symbol_from_atomic_number(nums[i]).to_string()
            } else if let Some(ref m) = masses {
                symbol_from_mass(m[i]).to_string()
            } else {
                name.clone()
            };

            let mut atom = Atom::with_symbol(&name, &symbol);
            if let Some(ref m) = masses {
                atom.mass = m[i];
            }
            if let Some(charges) = charges.as_ref() {
                atom.charge = charges[i] * 0.054_877_8;
            }
            topology.add_atom(atom);
        }

        if let Some(residue_labels) = residue_labels
            && let Some(residue_pointers) = residue_pointers
        {
            // AMBER residue pointers are 1-based. The end of residue i is the
            // start of residue i+1; for the last residue it is NATOM+1.
            for (i, label) in residue_labels.iter().enumerate() {
                let start = residue_pointers[i];
                let end = residue_pointers
                    .get(i + 1)
                    .copied()
                    .unwrap_or(pointers.NATOM + 1);
                let mut residue = Residue::new_from_name(label);
                residue.id = Some(i as i64 + 1);
                for id in start..end {
                    residue.add_atom((id - 1) as usize);
                }
                topology.add_residue(residue)?;
            }
        }

        let all_bonds: Vec<usize> = bonds_h.into_iter().chain(bonds_noh).flatten().collect();
        for chunk in all_bonds.chunks(3) {
            // AMBER encodes atom indices as atom_index * 3 + 1 (FORTRAN coordinate array offset).
            // as indices are 1-based, we omit the +1 -1.
            let from = chunk[0] / 3;
            let to = chunk[1] / 3;
            topology.add_bond(from, to, BondOrder::Unknown)?;
        }

        Ok(Some(topology))
    }

    fn skip_section(
        line_buffer: &mut String,
        file_reader: &mut BufReader<File>,
        already_read_flag: &mut bool,
    ) -> bool {
        while file_reader.read_line(line_buffer).unwrap_or(0) > 0 {
            if line_buffer.contains("%FLAG") {
                *already_read_flag = true;
                return false;
            }
            line_buffer.clear();
        }

        // we reached EOF
        true
    }

    pub fn len(&self) -> Result<usize, CError> {
        match self.convention {
            Convention::Restart => Ok(1),
            Convention::Amber => self
                .reader
                .as_ref()
                .expect("we should have init reader")
                .data_set()
                .num_records()
                .ok_or(CError::GenericError(
                    "expected unlimited 'frame' dimension to be defined".to_string(),
                )),
        }
    }

    pub fn is_empty(&self) -> Result<bool, CError> {
        self.len().map(|n| n == 0)
    }

    pub fn write(&mut self, frame: &Frame) -> Result<(), CError> {
        if !self.initialized {
            self.n_atoms = frame.size();
            self.has_velocities = frame.velocities().is_some();
            if let Some(Property::String(name)) = frame.properties.get("name") {
                self.file_title.clone_from(name);
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
            .and_then(super::super::property::Property::as_double)
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
            Some(p) => p.clone(),
            None => return Ok(()),
        };

        if self.buffered_frames.is_empty() {
            return Ok(());
        }

        let existing_count = self.index - self.buffered_frames.len();
        if self.mode == FileMode::Append && existing_count > 0 {
            self.finish_append(&path)
        } else {
            self.finish_create(&path)
        }
    }

    fn finish_create(&mut self, path: &Path) -> Result<(), CError> {
        let has_time = self.buffered_frames.iter().any(|f| f.time.is_some());
        let total_frames = self.buffered_frames.len();

        let data_set = self.build_data_set(total_frames, has_time, self.has_velocities)?;

        let mut writer = FileWriter::open(path)?;
        writer.set_def(&data_set, Version::Offset64Bit, 0)?;

        writer.write_var_u8("spatial", b"xyz")?;
        writer.write_var_u8("cell_spatial", b"abc")?;
        writer.write_var_u8("cell_angular", b"alphabeta gamma")?;

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

    fn finish_append(&mut self, path: &Path) -> Result<(), CError> {
        let new_count = self.buffered_frames.len();
        let existing_count = self.index - new_count;
        let total_count = self.index;

        let has_velocities = self.variables.velocities.is_some();
        let has_time = self.variables.time.is_some();
        let version = self.version.clone().unwrap_or(Version::Offset64Bit);

        let data_set = self.build_data_set(total_count, has_time, has_velocities)?;

        let output_file = OpenOptions::new().write(true).open(path)?;
        let mut writer = FileWriter::open_seek_write(
            path.as_os_str().to_str().expect("valid unicode path"),
            Box::new(output_file),
        )?;

        // Rewrite the header at offset 0 with the updated record count; the data
        // area already contains the existing records and is left intact.
        writer.set_def(&data_set, version, 0)?;

        // Write only the new records at their correct positions after the existing ones.
        for (i, frame) in self.buffered_frames.iter().enumerate() {
            let record_index = existing_count + i;
            writer.write_record_f32("coordinates", record_index, &frame.positions)?;
            writer.write_record_f32("cell_lengths", record_index, &frame.cell_lengths)?;
            writer.write_record_f32("cell_angles", record_index, &frame.cell_angles)?;
            if let Some(ref velocities) = frame.velocities {
                writer.write_record_f32("velocities", record_index, velocities)?;
            }
            if let Some(time) = frame.time {
                writer.write_record_f32("time", record_index, &[time])?;
            }
        }
        // Drop writer without calling close(): close() would fill unwritten record
        // slots (0..existing_count) with NC_FILL, overwriting the existing data.

        self.buffered_frames.clear();
        Ok(())
    }

    fn build_data_set(
        &self,
        total_frames: usize,
        has_time: bool,
        has_velocities: bool,
    ) -> Result<DataSet, CError> {
        let mut data_set = DataSet::new();
        data_set.add_global_attr_string("Conventions", "AMBER")?;
        data_set.add_global_attr_string("ConventionVersion", "1.0")?;
        data_set.add_global_attr_string("program", "molio")?;
        data_set.add_global_attr_string("programVersion", env!("CARGO_PKG_VERSION"))?;

        if !self.file_title.is_empty() {
            data_set.add_global_attr_string("title", &self.file_title)?;
        }

        data_set.set_unlimited_dim("frame", total_frames)?;
        data_set.add_fixed_dim("spatial", 3)?;
        data_set.add_fixed_dim("atom", self.n_atoms)?;
        data_set.add_fixed_dim("cell_spatial", 3)?;
        data_set.add_fixed_dim("cell_angular", 3)?;
        data_set.add_fixed_dim("label", 5)?;

        data_set.add_var_u8("spatial", &["spatial"])?;
        data_set.add_var_u8("cell_spatial", &["cell_spatial"])?;
        data_set.add_var_u8("cell_angular", &["cell_angular", "label"])?;

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

        if has_velocities {
            data_set.add_var_f32("velocities", &["frame", "atom", "spatial"])?;
            data_set.add_var_attr_string("velocities", "units", "angstrom/picosecond")?;
        }

        Ok(data_set)
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use assert_approx_eq::assert_approx_eq;
    use tempfile::Builder;

    use crate::{
        bond::Bond,
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
        assert_approx_eq!(positions[0][0], 0.4172_191, 1e-4);
        assert_approx_eq!(positions[0][1], 8.303_366, 1e-4);
        assert_approx_eq!(positions[0][2], 11.737_17, 1e-4);

        assert_approx_eq!(positions[296][0], 6.664_049, 1e-4);
        assert_approx_eq!(positions[296][1], 11.614_18, 1e-4);
        assert_approx_eq!(positions[296][2], 12.961_49, 1e-4);

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
        assert_approx_eq!(positions[0][0], 0.299_095_2, 1e-4);
        assert_approx_eq!(positions[0][1], 8.31003, 1e-4);
        assert_approx_eq!(positions[0][2], 11.72146, 1e-4);

        assert_approx_eq!(positions[296][0], 6.797_599, 1e-4);
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
        assert_approx_eq!(positions[0][0], 0.318_558_6, 1e-4);
        assert_approx_eq!(positions[0][1], 8.776_042, 1e-4);
        assert_approx_eq!(positions[0][2], 11.8927, 1e-4);

        assert_approx_eq!(positions[296][0], 7.089_802, 1e-4);
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
    fn read_topology() {
        let path = Path::new("./src/tests-data/amber/7OAP_BA4_dry.nc");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.len(), 500);
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.len(), 4931);

        let topology = frame.topology();
        assert_eq!(topology.len(), 4931);
        assert_eq!(topology[0].name, "N");
        assert_eq!(topology[0].symbol, "N");
        assert_eq!(topology[1].name, "H1");
        assert_eq!(topology[1].symbol, "H");
        assert_eq!(topology[2].name, "H2");
        assert_eq!(topology[2].symbol, "H");
        assert_eq!(topology[39].name, "CB");
        assert_eq!(topology[39].symbol, "C");
        assert_eq!(topology.atoms.last().unwrap().name, "O");
        assert_eq!(topology.atoms.last().unwrap().symbol, "O");

        // read the charges
        assert_approx_eq!(topology[0].charge, 0.1493, 1e-5);
        assert_approx_eq!(topology[10].charge, 0.0331, 1e-5);
        assert_approx_eq!(topology[49].charge, 0.4251, 1e-5);
        assert_approx_eq!(topology.atoms.last().unwrap().charge, -0.5679, 1e-5);

        assert_eq!(topology.residues.len(), 320);

        let residue = &topology.residues[0];
        assert_eq!(residue.len(), 19);
        assert!(residue.contains(0));
        assert!(residue.contains(10));
        assert!(residue.contains(12));

        assert_eq!(residue.name, "GLN");
        assert_eq!(topology.residues[35].name, "TRP");
        assert_eq!(topology.residues.last().unwrap().name, "SER");

        assert_eq!(residue.id.unwrap(), 1);
        assert_eq!(topology.residues.last().unwrap().id.unwrap(), 320);

        let bonds = topology.bonds();
        // NBONH=2401 (bonds with H) + MBONA=2606 (bonds without H)
        assert_eq!(bonds.len(), 5007);
        // N-terminus: N(0) bonded to H1(1) and H2(2)
        assert!(bonds.contains(&Bond::new(0, 1)));
        assert!(bonds.contains(&Bond::new(0, 2)));
        // randomly assorted
        assert!(bonds.contains(&Bond::new(1661, 1662)));
        assert!(bonds.contains(&Bond::new(2595, 2594)));

        // Check that there are NO bonds between the following atoms
        assert!(!bonds.contains(&Bond::new(3436, 3437)));
        assert!(!bonds.contains(&Bond::new(1588, 1627)));
        assert!(!bonds.contains(&Bond::new(2607, 3417)));
    }

    #[test]
    fn read_topology_no_atomic_number() {
        // the file uses Hydrogen Mass Repartioning so it's difficult to directly estimate the masses
        // therefore, in the following, we only check the name
        let path = Path::new("./src/tests-data/amber/7OAP_BA4_dry_no_atomic_number.nc");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.len(), 500);
        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.len(), 4931);

        let topology = frame.topology();
        assert_eq!(topology.len(), 4931);
        assert_eq!(topology[0].name, "N");
        assert_eq!(topology[0].symbol, "");
        assert_eq!(topology[1].name, "H1");
        assert_eq!(topology[1].symbol, "");
        assert_eq!(topology[2].name, "H2");
        assert_eq!(topology[2].symbol, "");
        assert_eq!(topology[39].name, "CB");
        assert_eq!(topology[39].symbol, "");

        // the last atom actually has the correct weight
        assert_eq!(topology.atoms.last().unwrap().name, "O");
        assert_eq!(topology.atoms.last().unwrap().symbol, "O");
    }

    #[test]
    fn scale_factor() {
        let path = Path::new("./src/tests-data/netcdf/scaled_traj.nc");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.len(), 26);

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
        assert_approx_eq!(velocities[1400][0], 0.6854_072 * -0.856, 1e-4);
        assert_approx_eq!(velocities[1400][1], 0.091_960_11 * -0.856, 1e-4);
        assert_approx_eq!(velocities[1400][2], 2.260_214 * -0.856, 1e-4);

        assert_approx_eq!(velocities[1600][0], -0.3342_645 * -0.856, 1e-4);
        assert_approx_eq!(velocities[1600][1], 0.322_594 * -0.856, 1e-4);
        assert_approx_eq!(velocities[1600][2], -2.446_901 * -0.856, 1e-4);
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

    fn generate_frame() -> Frame {
        let mut frame = Frame::from_unitcell(
            UnitCell::new_from_lengths_angles([2.0, 3.0, 4.0].into(), [80.0, 90.0, 120.0].into())
                .unwrap(),
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
                [1.0 * f64::from(i), 2.0 * f64::from(i), 3.0 * f64::from(i)],
                [-3.0, -2.0, -1.0],
            );
        }

        frame
    }

    #[test]
    fn write_files_in_netcdf_format() {
        let frame = generate_frame();
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

    #[test]
    fn append_to_existing_file() {
        let frame = generate_frame();
        let named_tmpfile = Builder::new()
            .prefix("netcdf-test-append")
            .suffix(".nc")
            .tempfile()
            .unwrap();

        {
            let mut trajectory = Trajectory::create(named_tmpfile.path()).unwrap();
            trajectory.write(&frame).unwrap();
        }

        {
            let mut trajectory = Trajectory::append(named_tmpfile.path()).unwrap();
            trajectory.write(&frame).unwrap();
        }

        let mut trajectory = Trajectory::open(named_tmpfile.path()).unwrap();
        assert_eq!(trajectory.len(), 2);
        check_frame(&trajectory.read().unwrap().unwrap());
        check_frame(&trajectory.read().unwrap().unwrap());
    }

    #[test]
    fn append_to_new_file() {
        let frame = generate_frame();
        let named_tmpfile = Builder::new()
            .prefix("netcdf-test-append")
            .suffix(".nc")
            .tempfile()
            .unwrap();

        {
            let mut trajectory = Trajectory::append(named_tmpfile.path()).unwrap();
            trajectory.write(&frame).unwrap();
            trajectory.write(&frame).unwrap();
        }

        let mut trajectory = Trajectory::open(named_tmpfile.path()).unwrap();
        assert_eq!(trajectory.len(), 2);
        check_frame(&trajectory.read().unwrap().unwrap());
        check_frame(&trajectory.read().unwrap().unwrap());
    }

    // AMBER restart format

    #[test]
    fn rst_water() {
        let path = Path::new("./src/tests-data/netcdf/water.ncrst");
        let mut trajectory = Trajectory::open(path).unwrap();
        assert_eq!(trajectory.len(), 1);

        let frame = trajectory.read().unwrap().unwrap();
        assert_eq!(frame.size(), 297);
        assert_eq!(
            frame.properties.get("name").unwrap().as_string().unwrap(),
            "Cpptraj Generated Restart"
        );
    }
}
