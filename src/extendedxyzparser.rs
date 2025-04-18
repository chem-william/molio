use std::{collections::HashMap, str::Chars};

/// The two kinds of values we can see in extended‑XYZ comments:
/// - `Bool(true)` for flags with no “=value` (ASE treats them as true)
/// - `Str(...)` for quoted or unquoted string values
#[derive(Debug, Clone, PartialEq)]
pub enum Property {
    Bool(bool),
    Str(String),
}

pub type ExtendedXyzProperties = HashMap<String, Property>;

pub struct ExtendedXyzParser<'a> {
    line: &'a str,
    pos: usize,
    done: bool,
}

impl<'a> ExtendedXyzParser<'a> {
    pub fn new(line: &'a str) -> Self {
        ExtendedXyzParser {
            line,
            pos: 0,
            done: false,
        }
    }

    pub fn parse(&mut self) -> ExtendedXyzProperties {
        let mut props = HashMap::new();

        while !self.done {
            let name = self.next_substring();
            let value = self.next_substring();
            props.insert(name, Property::Str(value.to_string()));
            self.done = true;
        }
        props
    }

    /// Grabs either a quoted substring (allowing spaces) or an unquoted
    /// one (stopping at `=` or whitespace).  Doesn't include the quotes.
    pub fn next_substring(&mut self) -> String {
        let line = self.line.trim();

        let mut char_indices = line[self.pos..].char_indices().peekable();

        let result = if line[self.pos..].starts_with('"') {
            // Quoted string: skip the opening quote
            char_indices.next();
            let quoted_start = self.pos + 1;

            // Only a quote at the end
            if self.pos >= line.len() {
                return "".to_string();
            }

            let pos_quoted_end = line[quoted_start..].chars().position(|x| x == '"');

            // Update position to after the closing quote (if any)
            if let Some(pos) = pos_quoted_end {
                self.pos += pos;
                line[quoted_start..(quoted_start + pos)].to_string()
            } else {
                self.pos = line.len();
                line[quoted_start..line.len()].to_string()
            }
        } else {
            // Unquoted string: ends at whitespace or '='
            let start_pos = self.pos;
            let pos_unquoted_end = line[self.pos..]
                .chars()
                .position(|x| x.is_whitespace() || x == '=');

            if let Some(pos) = pos_unquoted_end {
                self.pos += pos + 1;
                line[start_pos..pos].to_string()
            } else {
                self.pos = line.len();
                line[start_pos..line.len()].to_string()
            }
        };

        result
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn simple_flags_and_values() {
//         let mut p =
//             ExtendedXyzParser::new(r#"Properties=species:S:1:pos:R:3 name='test file' debug"#);
//         let m = p.parse();
//         assert_eq!(m["Properties"], Property::Str("species:S:1:pos:R:3".into()));
//         assert_eq!(m["name"], Property::Str("test file".into()));
//         assert_eq!(m["debug"], Property::Bool(true));
//     }
// }
