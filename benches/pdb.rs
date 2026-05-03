use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::Path;

fn criterion_benchmark(c: &mut Criterion) {
    let pdb_path = Path::new("./src/tests-data/pdb/water.pdb");
    c.bench_function("read water.pdb", |b| {
        b.iter(|| black_box(molio::read_trajectory(pdb_path)));
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
