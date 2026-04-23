use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::Path;

fn criterion_benchmark(c: &mut Criterion) {
    // Benchmark TNG format (large file: 38376 atoms × 6 frames)
    let tng_path = Path::new("./src/tests-data/tng/1aki.tng");
    c.bench_function("read 1aki.tng", |b| {
        b.iter(|| black_box(molio::read_trajectory(tng_path)))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
