use crate::property::Properties;

#[derive(Debug, Default, Clone)]
pub struct Atom {
    pub symbol: String,
    pub name: String,
    // pub mass: f64,
    // pub charge: f64,
    pub properties: Properties,
}

impl Atom {
    pub fn new(symbol: String, name: String) -> Self {
        Self {
            symbol,
            name,
            properties: Properties::new(),
        }
    }
}
