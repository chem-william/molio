use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::Path;

fn criterion_benchmark(c: &mut Criterion) {
    let amber_path = Path::new("./src/tests-data/amber/7OAP_BA4_dry.nc");
    c.bench_function("read 7OAP_BA4_dry.nc", |b| {
        b.iter(|| black_box(molio::read_trajectory(amber_path)));
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
