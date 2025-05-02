// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2025 William Bro-Jørgensen
// Copyright (c) 2020 Guillaume Fraux and contributors
//
// See LICENSE at the project root for full text.

use crate::property::Properties;

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Atom {
    pub name: String,
    pub symbol: String,
    // pub mass: f64,
    pub charge: f64,
    pub properties: Properties,
}

impl Atom {
    pub fn new(name: String) -> Self {
        Self {
            name: name.clone(),
            symbol: name,
            charge: 0.0,
            properties: Properties::new(),
        }
    }

    pub fn with_symbol(name: String, symbol: String) -> Self {
        Self {
            name,
            symbol,
            charge: 0.0,
            properties: Properties::new(),
        }
    }
}
