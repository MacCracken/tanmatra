use criterion::{Criterion, black_box, criterion_group, criterion_main};

fn binding_energy_1000(c: &mut Criterion) {
    c.bench_function("nucleus/binding_energy_1000", |b| {
        b.iter(|| {
            for z in 1u16..=100 {
                for a_offset in 0u16..10 {
                    let a = z + a_offset;
                    let _ = black_box(tanmatra::nucleus::binding_energy_mev(z, a));
                }
            }
        });
    });
}

fn spectral_line_1000(c: &mut Criterion) {
    c.bench_function("atomic/spectral_line_1000", |b| {
        b.iter(|| {
            for z in 1u16..=10 {
                for n_upper in 2u16..=101 {
                    let _ = black_box(tanmatra::atomic::spectral_line_nm(z, n_upper, 1));
                }
            }
        });
    });
}

fn electron_config_all_36(c: &mut Criterion) {
    c.bench_function("atomic/electron_config_all_36", |b| {
        b.iter(|| {
            for z in 1u16..=36 {
                let _ = black_box(tanmatra::atomic::electron_configuration(z));
            }
        });
    });
}

fn decay_chain_u238_10steps(c: &mut Criterion) {
    let u238 = tanmatra::decay::uranium238();
    c.bench_function("decay/decay_chain_u238_10steps", |b| {
        b.iter(|| {
            let _ = black_box(tanmatra::decay::decay_chain(&u238, 10));
        });
    });
}

criterion_group!(
    benches,
    binding_energy_1000,
    spectral_line_1000,
    electron_config_all_36,
    decay_chain_u238_10steps,
);
criterion_main!(benches);
