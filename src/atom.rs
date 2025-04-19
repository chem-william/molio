use crate::property::Properties;

#[derive(Debug)]
pub struct Atom {
    pub id: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub properties: Properties,
}
