use nalgebra::Matrix3;
use std::collections::HashMap;

#[derive(PartialEq, Clone, Debug)]
pub enum PropertyKind {
    Bool,
    Double,
    String,
    Vector3D,
    Matrix3x3,
    VectorXD,
}

#[derive(Debug)]
pub enum Property {
    Bool(bool),
    Double(f64),
    String(String),
    Vector3D([f64; 3]),
    Matrix3x3(Matrix3<f64>),
    VectorXD(Vec<f64>),
}

pub type Properties = HashMap<String, Property>;

impl Property {
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
            ref other => panic!("expected Bool, found {:?}", other),
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
            ref other => panic!("expected Double, found {:?}", other),
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
            ref other => panic!("expected String, found {:?}", other),
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
            ref other => panic!("expected Vector3D, found {:?}", other),
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
            ref other => panic!("expected Matrix3x3, found {:?}", other),
        }
    }
}

#[derive(Debug)]
pub struct ExtendedProperty {
    pub name: String,
    pub kind: PropertyKind,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bool_property() {
        let prop = Property::Bool(true);
        assert_eq!(prop.as_bool(), Some(true));
        assert_eq!(prop.expect_bool(), true);
        assert_eq!(prop.as_double(), None);
    }

    #[test]
    fn test_double_property() {
        let prop = Property::Double(3.14);
        assert_eq!(prop.as_double(), Some(3.14));
        assert_eq!(prop.expect_double(), 3.14);
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
        let matrix = Matrix3::new(
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0
        );
        let prop = Property::Matrix3x3(matrix);
        assert_eq!(prop.as_matrix3x3(), Some([1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0]));
        assert_eq!(prop.expect_matrix3x3(), [1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0]);
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
