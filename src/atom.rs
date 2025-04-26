use crate::property::Properties;

#[derive(Debug, Default, Clone)]
pub struct Atom {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub symbol: String,
    pub name: String,
    // pub mass: f64,
    // pub charge: f64,
    pub properties: Properties,
}

impl Atom {
    pub fn new(x: f64, y: f64, z: f64, symbol: String, name: String) -> Self {
        Self {
            x,
            y,
            z,
            symbol,
            name,
            properties: Properties::new(),
        }
    }
}
