use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::Path;

fn criterion_benchmark(c: &mut Criterion) {
    let sdf_path = Path::new("./src/tests-data/sdf/kinases.sdf");
    c.bench_function("read kinases.sdf", |b| {
        b.iter(|| black_box(molio::read_trajectory(sdf_path)));
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
