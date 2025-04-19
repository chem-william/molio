use criterion::{Criterion, criterion_group, criterion_main};
use molio::trajectory::Trajectory;
use std::path::Path;
use std::time::Duration;

fn load_trajectory(path: &Path) -> usize {
    let trajectory = Trajectory::new(path).unwrap();
    trajectory.size
}

fn criterion_benchmark(c: &mut Criterion) {
    let path = Path::new("./src/tests-data/xyz/helium.xyz");
    let mut group = c.benchmark_group("my_group");
    group.measurement_time(Duration::from_secs(6));
    group.bench_function("read large xyz", |b| b.iter(|| load_trajectory(path)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
