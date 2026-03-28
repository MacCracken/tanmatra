use criterion::{Criterion, black_box, criterion_group, criterion_main};
use tanmatra::prelude::*;

fn binding_energy_1000(c: &mut Criterion) {
    c.bench_function("nucleus/binding_energy_1000", |b| {
        let nuclei: Vec<Nucleus> = (1..=100)
            .flat_map(|z| {
                let a_min = z;
                let a_max = z * 3;
                (a_min..=a_max.min(a_min + 9)).filter_map(move |a| Nucleus::new(z, a).ok())
            })
            .collect();
        b.iter(|| {
            for n in &nuclei {
                black_box(n.binding_energy());
            }
        });
    });
}

fn spectral_line_1000(c: &mut Criterion) {
    c.bench_function("atomic/spectral_line_1000", |b| {
        b.iter(|| {
            for z in 1..=10 {
                for n1 in 1..=10 {
                    for n2 in (n1 + 1)..=10u32.min(n1 + 10) {
                        black_box(spectral_line_nm(z, n1, n2).ok());
                    }
                }
            }
        });
    });
}

fn electron_config_36(c: &mut Criterion) {
    c.bench_function("atomic/electron_config_36", |b| {
        b.iter(|| {
            for z in 1..=36 {
                black_box(electron_configuration(z).ok());
            }
        });
    });
}

fn decay_chain_10(c: &mut Criterion) {
    c.bench_function("decay/decay_chain_10", |b| {
        let u238 = Nucleus::uranium_238();
        b.iter(|| {
            black_box(decay_chain(&u238, 10));
        });
    });
}

criterion_group!(
    benches,
    binding_energy_1000,
    spectral_line_1000,
    electron_config_36,
    decay_chain_10,
);
criterion_main!(benches);
