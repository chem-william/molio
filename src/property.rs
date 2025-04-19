use nalgebra::{Matrix3, VectorN};

// #[derive(Debug, Clone, PartialEq)]
// pub enum PropertyKind {
//     Bool(bool),
//     Double(f64),
//     Str(String),
//     Vector3D(f64),
// }
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
pub enum AtomProperty {
    Bool(bool),
    Double(f64),
    String(String),
    Vector3D([f64; 3]),
    Matrix3x3(Matrix3<f64>),
    VectorXD(Vec<f64>),
}
impl AtomProperty {
    pub fn as_bool(&self) -> Option<bool> {
        if let AtomProperty::Bool(b) = *self {
            Some(b)
        } else {
            None
        }
    }
    pub fn expect_bool(&self) -> bool {
        match *self {
            AtomProperty::Bool(b) => b,
            ref other => panic!("expected Bool, found {:?}", other),
        }
    }
    pub fn as_double(&self) -> Option<f64> {
        if let AtomProperty::Double(x) = *self {
            Some(x)
        } else {
            None
        }
    }
    pub fn expect_double(&self) -> f64 {
        match *self {
            AtomProperty::Double(d) => d,
            ref other => panic!("expected Double, found {:?}", other),
        }
    }
    pub fn as_string(&self) -> Option<&str> {
        if let AtomProperty::String(ref s) = *self {
            Some(s)
        } else {
            None
        }
    }
    pub fn expect_string(&self) -> &str {
        match *self {
            AtomProperty::String(ref s) => s,
            ref other => panic!("expected String, found {:?}", other),
        }
    }
    pub fn as_vector3d(&self) -> Option<[f64; 3]> {
        if let AtomProperty::Vector3D(v) = *self {
            Some(v)
        } else {
            None
        }
    }
    pub fn expect_vector3d(&self) -> [f64; 3] {
        match *self {
            AtomProperty::Vector3D(v) => v,
            ref other => panic!("expected Vector3D, found {:?}", other),
        }
    }
    pub fn as_matrix3x3(&self) -> Option<[f64; 9]> {
        if let AtomProperty::Matrix3x3(m) = *self {
            let mut array = [0.0; 9];
            array.copy_from_slice(m.as_slice());
            Some(array)
        } else {
            None
        }
    }
    pub fn expect_matrix3x3(&self) -> [f64; 9] {
        match *self {
            AtomProperty::Matrix3x3(m) => {
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
