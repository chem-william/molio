use crate::property::Property;
use std::collections::HashMap;

pub type ExtendedXyzProperties = HashMap<String, Property>;

pub struct ExtendedXyzParser<'a> {
    line: &'a str,
    pos: usize,
    done: bool,
}
#[derive(Debug, PartialEq)]
enum State {
    /// We're not in a token right now
    Whitespace,
    /// We're building a token, unquoted
    Unquoted,
    /// We're inside a double‐quoted token
    DoubleQuoted,
    /// We're inside a single‐quoted token
    SingleQuoted,
}

impl<'a> ExtendedXyzParser<'a> {
    pub fn new(line: &'a str) -> Self {
        ExtendedXyzParser {
            line,
            pos: 0,
            done: false,
        }
    }

    /// Splits an input string into tokens, honoring single and double quotes.
    fn split_tokens(input: &str) -> Vec<String> {
        let mut tokens = Vec::new();
        let mut current = String::new();
        let mut state = State::Whitespace;

        for ch in input.chars() {
            match state {
                State::Whitespace => {
                    if ch.is_whitespace() {
                        // skip
                    } else if ch == '"' {
                        state = State::DoubleQuoted;
                    } else if ch == '\'' {
                        state = State::SingleQuoted;
                    } else {
                        state = State::Unquoted;
                        current.push(ch);
                    }
                }
                State::Unquoted => {
                    if ch.is_whitespace() {
                        // end token
                        tokens.push(current.clone());
                        current.clear();
                        state = State::Whitespace;
                    } else if ch == '"' {
                        state = State::DoubleQuoted;
                    } else if ch == '\'' {
                        state = State::SingleQuoted;
                    } else {
                        current.push(ch);
                    }
                }
                State::DoubleQuoted => {
                    if ch == '"' {
                        // end of double‐quoted section
                        state = State::Unquoted;
                    } else {
                        current.push(ch);
                    }
                }
                State::SingleQuoted => {
                    if ch == '\'' {
                        // end of single‐quoted section
                        state = State::Unquoted;
                    } else {
                        current.push(ch);
                    }
                }
            }
        }

        // if we ended while building a token, push it
        if state != State::Whitespace {
            tokens.push(current);
        }

        tokens
    }

    /// Parses a single line of KEY=VALUE and standalone flags into a HashMap.
    pub fn parse(&self) -> ExtendedXyzProperties {
        let mut map = HashMap::new();
        for tok in ExtendedXyzParser::split_tokens(self.line) {
            if let Some(idx) = tok.find('=') {
                let key = &tok[..idx];
                let mut val = tok[idx + 1..].to_string();
                // strip surrounding quotes if present
                if (val.starts_with('"') && val.ends_with('"'))
                    || (val.starts_with('\'') && val.ends_with('\''))
                {
                    val = val[1..val.len() - 1].to_string();
                }
                map.insert(key.to_string(), Property::Str(val));
            } else {
                // no '=', treat as boolean flag
                map.insert(tok, Property::Bool(true));
            }
        }
        map
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_flags_and_values() {
        let p = ExtendedXyzParser::new(r#"Properties=species:S:1:pos:R:3 name='test file' debug"#);
        let m = p.parse();
        assert_eq!(m["Properties"], Property::Str("species:S:1:pos:R:3".into()));
        assert_eq!(m["name"], Property::Str("test file".into()));
        assert_eq!(m["debug"], Property::Bool(true));
    }
}
