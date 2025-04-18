// #[derive(Debug, Clone, PartialEq)]
// pub enum PropertyKind {
//     Bool(bool),
//     Double(f64),
//     Str(String),
//     Vector3D(f64),
// }
#[derive(PartialEq, Clone, Copy, Debug)]
pub enum PropertyKind {
    Bool,
    Double,
    String,
    Vector3D,
}
#[derive(Debug)]
pub struct ExtendedProperty {
    pub name: String,
    pub kind: PropertyKind,
}
