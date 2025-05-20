use std::{env, path::Path};

fn main() {
    let args: Vec<String> = env::args().collect();

    // Check for profiling mode flag
    let mut profiling_mode = false;
    let mut file_path = "./water.pdb";

    // Parse command line arguments
    for arg in &args[1..] {
        if arg == "--profile" || arg == "-p" {
            profiling_mode = true;
        } else if !arg.starts_with('-') {
            // Assume this is a file path
            file_path = arg;
        }
    }

    let path = Path::new(file_path);

    if profiling_mode {
        // Profiling mode: run multiple iterations for flamegraph
        println!("Running in profiling mode (32 iterations)");
        for _ in 0..32 {
            molio::read_trajectory(path);
        }
    } else {
        // Normal benchmarking mode: run once
        molio::read_trajectory(path);
    }
}
