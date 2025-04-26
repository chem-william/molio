use std::{hint::black_box, path::Path};

use molio::trajectory::Trajectory;

fn main() {
    // let path = Path::new("./helium.xyz");
    let path = Path::new("./src/tests-data/xyz/helium.xyz");
    // let path = Path::new("./src/tests-data/pdb/water.pdb");
    // let path = Path::new("./water.pdb");
    for _ in 0..32 {
        let mut trajectory = Trajectory::open(path).unwrap();
        let mut dummy = 0;
        while let Some(next_frame) = trajectory.read().unwrap() {
            dummy += next_frame.size();
        }
        black_box(dummy);
    }
}
