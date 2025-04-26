use crate::property::Properties;

#[derive(Debug, Default, Clone)]
pub struct Atom {
    pub name: String,
    pub symbol: String,
    // pub mass: f64,
    // pub charge: f64,
    pub properties: Properties,
}

impl Atom {
    pub fn new(name: String) -> Self {
        Self {
            name: name.clone(),
            symbol: name,
            properties: Properties::new(),
        }
    }

    pub fn with_symbol(name: String, symbol: String) -> Self {
        Self {
            name,
            symbol,
            properties: Properties::new(),
        }
    }
}
