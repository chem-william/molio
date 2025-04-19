use crate::atom::Atom;
use crate::property::Properties;
use crate::unit_cell::UnitCell;

#[derive(Debug, Default)]
pub struct Frame {
    pub atoms: Vec<Atom>,
    pub unit_cell: UnitCell,
    pub properties: Properties,
}

impl Frame {
    pub fn new() -> Self {
        Frame {
            atoms: vec![],
            unit_cell: UnitCell::new(),
            properties: Properties::new(),
        }
    }
    pub fn size(&self) -> usize {
        self.atoms.len()
    }

    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.atoms.iter().map(|a| [a.x, a.y, a.z]).collect()
    }
    pub fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom)
    }
}
