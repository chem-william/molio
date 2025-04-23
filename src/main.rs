use std::{hint::black_box, path::Path};

use molio::trajectory::Trajectory;

fn main() {
    let path = Path::new("./helium.xyz");
    let mut trajectory = Trajectory::new(path).unwrap();
    let mut dummy = 0;
    while let Some(next_frame) = trajectory.read().unwrap() {
        dummy += next_frame.size();
    }
    black_box(dummy);
}
