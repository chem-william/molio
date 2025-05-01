use criterion::{criterion_group, criterion_main, Criterion};
use std::path::Path;
use std::time::Duration;

// Import the shared function from the main module
fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("format_reading");
    group.measurement_time(Duration::from_secs(6));

    // Benchmark XYZ format
    let xyz_path = Path::new("./src/tests-data/xyz/helium.xyz");
    group.bench_function("read helium.xyz", |b| {
        b.iter(|| molio::read_trajectory(xyz_path))
    });

    // Benchmark PDB format
    let pdb_path = Path::new("./src/tests-data/pdb/water.pdb");
    group.bench_function("read water.pdb", |b| {
        b.iter(|| molio::read_trajectory(pdb_path))
    });

    // Benchmark SDF format
    let sdf_path = Path::new("./src/tests-data/sdf/kinases.sdf");
    group.bench_function("read kinases.sdf", |b| {
        b.iter(|| molio::read_trajectory(sdf_path))
    });

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
