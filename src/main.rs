use std::{hint::black_box, path::Path};

use molio::{frame::Frame, trajectory::Trajectory};

fn main() {
    let path = Path::new("./helium.xyz");
    let mut trajectory = Trajectory::new(path).unwrap();
    let mut frame = Frame::new();
    let mut dummy = 0;
    while let Some(next_frame) = trajectory.read().unwrap() {
        frame = next_frame;
        dummy += frame.size();
    }
    black_box(dummy);
}
