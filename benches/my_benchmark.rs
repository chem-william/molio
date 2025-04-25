use criterion::{Criterion, criterion_group, criterion_main};
use molio::trajectory::Trajectory;
use std::hint::black_box;
use std::path::Path;
use std::time::Duration;

fn load_trajectory(path: &Path) -> usize {
    let mut trajectory = Trajectory::new(path).unwrap();
    let mut dummy = 0;
    while let Some(next_frame) = trajectory.read().unwrap() {
        dummy += next_frame.size();
    }
    black_box(dummy)
}

fn criterion_benchmark(c: &mut Criterion) {
    let path = Path::new("./src/tests-data/xyz/helium.xyz");
    let mut group = c.benchmark_group("my_group");
    group.measurement_time(Duration::from_secs(6));
    group.bench_function("read large xyz", |b| b.iter(|| load_trajectory(path)));

    let path = Path::new("./src/tests-data/pdb/water.pdb");
    group.bench_function("read large pdb", |b| b.iter(|| load_trajectory(path)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
