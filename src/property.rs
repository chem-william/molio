use std::collections::HashMap;
use std::ops::{Deref, DerefMut};

use nalgebra::Matrix3;

use crate::error::CError;

const EPSILON: f64 = 1e-12;

#[derive(PartialEq, Clone, Debug)]
pub enum PropertyKind {
    Bool,
    Double,
    String,
    Vector3D,
    Matrix3x3,
    VectorXD,
}

#[derive(Debug, Clone)]
pub enum Property {
    Bool(bool),
    Double(f64),
    String(String),
    Vector3D([f64; 3]),
    Matrix3x3(Matrix3<f64>),
    VectorXD(Vec<f64>),
}

impl Default for Property {
    fn default() -> Self {
        Property::Bool(false)
    }
}

/// Returns `true` if `a` and `b` are both finite and within `epsilon` of each other.
/// Any `NaN` or infinite value always compares as `false`.
fn almost_eq(a: f64, b: f64, epsilon: f64) -> bool {
    // Reject NaN outright
    if a.is_nan() || b.is_nan() {
        return false;
    }
    // Reject infinities outright
    if a.is_infinite() || b.is_infinite() {
        return false;
    }
    // Both are finite: compare absolute difference
    (a - b).abs() <= epsilon
}

impl PartialEq for Property {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Property::Bool(a), Property::Bool(b)) => a == b,

            (Property::Double(a), Property::Double(b)) => almost_eq(*a, *b, EPSILON),

            (Property::String(a), Property::String(b)) => a == b,

            (Property::Vector3D(a), Property::Vector3D(b)) => a
                .iter()
                .zip(b.iter())
                .all(|(x, y)| almost_eq(*x, *y, EPSILON)),

            // nalgebra’s Matrix3<f64>: compare each entry’s bits
            (Property::Matrix3x3(a), Property::Matrix3x3(b)) => {
                for i in 0..3 {
                    for j in 0..3 {
                        if !almost_eq(a[(i, j)], b[(i, j)], EPSILON) {
                            return false;
                        }
                    }
                }
                true
            }

            // variable‐length vector of f64: same pattern
            (Property::VectorXD(a), Property::VectorXD(b)) => {
                if a.len() != b.len() {
                    return false;
                }
                a.iter()
                    .zip(b.iter())
                    .all(|(x, y)| almost_eq(*x, *y, EPSILON))
            }

            // different variants are never equal
            _ => false,
        }
    }
}

impl Eq for Property {}

#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct Properties(HashMap<String, Property>);

impl Properties {
    pub fn new() -> Self {
        Properties(HashMap::new())
    }
}

impl Deref for Properties {
    type Target = HashMap<String, Property>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Properties {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl IntoIterator for Properties {
    type Item = (String, Property);
    type IntoIter = <HashMap<String, Property> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl Property {
    #[must_use]
    pub fn as_bool(&self) -> Option<bool> {
        if let Property::Bool(b) = *self {
            Some(b)
        } else {
            None
        }
    }

    pub fn expect_bool(&self) -> bool {
        match *self {
            Property::Bool(b) => b,
            ref other => panic!("expected Bool, found {other:?}"),
        }
    }
    pub fn as_double(&self) -> Option<f64> {
        if let Property::Double(x) = *self {
            Some(x)
        } else {
            None
        }
    }
    pub fn expect_double(&self) -> f64 {
        match *self {
            Property::Double(d) => d,
            ref other => panic!("expected Double, found {other:?}"),
        }
    }
    pub fn as_string(&self) -> Option<&str> {
        if let Property::String(ref s) = *self {
            Some(s)
        } else {
            None
        }
    }
    pub fn expect_string(&self) -> &str {
        match *self {
            Property::String(ref s) => s,
            ref other => panic!("expected String, found {other:?}"),
        }
    }
    pub fn as_vector3d(&self) -> Option<[f64; 3]> {
        if let Property::Vector3D(v) = *self {
            Some(v)
        } else {
            None
        }
    }
    pub fn expect_vector3d(&self) -> [f64; 3] {
        match *self {
            Property::Vector3D(v) => v,
            ref other => panic!("expected Vector3D, found {other:?}"),
        }
    }
    pub fn as_matrix3x3(&self) -> Option<[f64; 9]> {
        if let Property::Matrix3x3(m) = *self {
            let mut array = [0.0; 9];
            array.copy_from_slice(m.as_slice());
            Some(array)
        } else {
            None
        }
    }
    pub fn expect_matrix3x3(&self) -> [f64; 9] {
        match *self {
            Property::Matrix3x3(m) => {
                let mut array = [0.0; 9];
                array.copy_from_slice(m.as_slice());
                array
            }
            ref other => panic!("expected Matrix3x3, found {other:?}"),
        }
    }

    pub fn parse_value(value: &str, kind: &PropertyKind) -> Result<Property, CError> {
        match kind {
            PropertyKind::String => StringParser::parse(value),
            PropertyKind::Bool => BoolParser::parse(value),
            PropertyKind::Double => DoubleParser::parse(value),
            PropertyKind::Vector3D => Vector3DParser::parse(value),
            PropertyKind::Matrix3x3 => Matrix3x3Parser::parse(value),
            PropertyKind::VectorXD => VectorXDParser::parse(value),
        }
    }
}
/// Helper trait for parsing values into Property
pub trait ValueParser {
    fn parse(value: &str) -> Result<Property, CError>;
}

pub struct StringParser;
pub struct BoolParser;
pub struct DoubleParser;
pub struct Vector3DParser;
pub struct Matrix3x3Parser;
pub struct VectorXDParser;

impl ValueParser for StringParser {
    fn parse(value: &str) -> Result<Property, CError> {
        Ok(Property::String(value.to_string()))
    }
}

impl ValueParser for BoolParser {
    fn parse(value: &str) -> Result<Property, CError> {
        match value.to_lowercase().as_str() {
            "t" | "true" => Ok(Property::Bool(true)),
            "f" | "false" => Ok(Property::Bool(false)),
            _ => Err(CError::GenericError(format!(
                "Invalid boolean value: {value}"
            ))),
        }
    }
}

impl ValueParser for DoubleParser {
    fn parse(value: &str) -> Result<Property, CError> {
        value
            .trim()
            .parse::<f64>()
            .map(Property::Double)
            .map_err(|e| CError::GenericError(format!("Failed to parse number: {e}")))
    }
}

impl ValueParser for Vector3DParser {
    fn parse(value: &str) -> Result<Property, CError> {
        let parts: Vec<&str> = value.split_whitespace().collect();
        if parts.len() != 3 {
            return Err(CError::GenericError(format!(
                "Vector3D requires exactly 3 components, got {}: {parts:?}",
                parts.len()
            )));
        }
        let x = parts[0]
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("Failed to parse x component: {e}")))?;
        let y = parts[1]
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("Failed to parse y component: {e}")))?;
        let z = parts[2]
            .parse::<f64>()
            .map_err(|e| CError::GenericError(format!("Failed to parse z component: {e}")))?;
        Ok(Property::Vector3D([x, y, z]))
    }
}

impl ValueParser for Matrix3x3Parser {
    fn parse(value: &str) -> Result<Property, CError> {
        let parts: Vec<&str> = value.split_whitespace().collect();
        if parts.len() != 9 {
            return Err(CError::GenericError(format!(
                "Matrix3x3 requires exactly 9 components, got {}: {parts:?}",
                parts.len()
            )));
        }
        let nums: Result<Vec<f64>, _> = parts.iter().map(|p| p.parse::<f64>()).collect();
        nums.map(|n| Property::Matrix3x3(Matrix3::from_iterator(n)))
            .map_err(|e| CError::GenericError(format!("Failed to parse matrix components: {e}")))
    }
}

impl ValueParser for VectorXDParser {
    fn parse(value: &str) -> Result<Property, CError> {
        let parts: Vec<&str> = value.split_whitespace().collect();
        let nums: Result<Vec<f64>, _> = parts.iter().map(|p| p.parse::<f64>()).collect();
        nums.map(Property::VectorXD)
            .map_err(|e| CError::GenericError(format!("Failed to parse vector components: {e}")))
    }
}

#[derive(Debug)]
pub struct ExtendedProperty {
    pub name: String,
    pub kind: PropertyKind,
}

#[cfg(test)]
mod tests {
    use core::f64;

    use assert_approx_eq::assert_approx_eq;

    use super::*;

    #[test]
    fn test_bool_property() {
        let prop = Property::Bool(true);
        assert_eq!(prop.as_bool(), Some(true));
        assert!(prop.expect_bool());
        assert_eq!(prop.as_double(), None);
    }

    #[test]
    fn test_double_property() {
        let prop = Property::Double(f64::consts::PI);
        assert_eq!(prop.as_double(), Some(f64::consts::PI));
        assert_approx_eq!(prop.expect_double(), f64::consts::PI);
        assert_eq!(prop.as_bool(), None);
    }

    #[test]
    fn test_string_property() {
        let prop = Property::String("test".to_string());
        assert_eq!(prop.as_string(), Some("test"));
        assert_eq!(prop.expect_string(), "test");
        assert_eq!(prop.as_double(), None);
    }

    #[test]
    fn test_vector3d_property() {
        let vec = [1.0, 2.0, 3.0];
        let prop = Property::Vector3D(vec);
        assert_eq!(prop.as_vector3d(), Some(vec));
        assert_eq!(prop.expect_vector3d(), vec);
        assert_eq!(prop.as_double(), None);
    }

    #[test]
    fn test_matrix3x3_property() {
        let matrix = Matrix3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        let prop = Property::Matrix3x3(matrix);
        let expected = [1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0];

        assert_eq!(
            prop.as_matrix3x3(),
            Some([1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0])
        );

        for (v, e) in expected.iter().zip(prop.expect_matrix3x3()) {
            assert_approx_eq!(v, e);
        }
        assert_eq!(prop.as_double(), None);
    }

    #[test]
    #[should_panic(expected = "expected Bool")]
    fn test_expect_bool_panic() {
        let prop = Property::Double(1.0);
        prop.expect_bool();
    }

    #[test]
    #[should_panic(expected = "expected Double")]
    fn test_expect_double_panic() {
        let prop = Property::String("test".to_string());
        prop.expect_double();
    }

    #[test]
    #[should_panic(expected = "expected String")]
    fn test_expect_string_panic() {
        let prop = Property::Double(1.0);
        prop.expect_string();
    }

    #[test]
    #[should_panic(expected = "expected Vector3D")]
    fn test_expect_vector3d_panic() {
        let prop = Property::Double(1.0);
        prop.expect_vector3d();
    }

    #[test]
    #[should_panic(expected = "expected Matrix3x3")]
    fn test_expect_matrix3x3_panic() {
        let prop = Property::Double(1.0);
        prop.expect_matrix3x3();
    }
}
