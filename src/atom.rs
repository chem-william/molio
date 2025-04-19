use crate::property::Properties;

#[derive(Debug)]
pub struct Atom {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub symbol: String,
    // pub mass: f64,
    // pub charge: f64,
    pub properties: Properties,
}
