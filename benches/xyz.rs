use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::Path;

fn criterion_benchmark(c: &mut Criterion) {
    let xyz_path = Path::new("./src/tests-data/xyz/helium.xyz");
    c.bench_function("read helium.xyz", |b| {
        b.iter(|| black_box(molio::read_trajectory(xyz_path)));
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
