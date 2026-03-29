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

fn shell_occupation_126(c: &mut Criterion) {
    c.bench_function("nucleus/shell_occupation_126", |b| {
        b.iter(|| {
            for n in 1..=126 {
                black_box(tanmatra::nucleus::shell_occupation(n));
            }
        });
    });
}

fn ionization_energy_118(c: &mut Criterion) {
    c.bench_function("atomic/ionization_energy_118", |b| {
        b.iter(|| {
            for z in 1..=118 {
                black_box(ionization_energy_ev(z).ok());
            }
        });
    });
}

fn lorentz_gamma_1000(c: &mut Criterion) {
    c.bench_function("relativity/lorentz_gamma_1000", |b| {
        b.iter(|| {
            for i in 1..=1000 {
                let beta = i as f64 / 1001.0;
                black_box(tanmatra::relativity::lorentz_gamma(beta));
            }
        });
    });
}

fn rutherford_scattering_1000(c: &mut Criterion) {
    c.bench_function("scattering/rutherford_1000", |b| {
        b.iter(|| {
            for i in 1..=1000 {
                let theta = i as f64 * core::f64::consts::PI / 1001.0;
                black_box(tanmatra::scattering::rutherford_differential(
                    2, 79, 5.0, theta,
                ));
            }
        });
    });
}

fn bateman_chain_3(c: &mut Criterion) {
    c.bench_function("decay/bateman_chain_3", |b| {
        let lambdas = [0.1, 0.05, 0.0];
        b.iter(|| {
            black_box(tanmatra::decay::bateman_chain(&lambdas, 1e6, 100.0));
        });
    });
}

fn radial_wavefunction_100(c: &mut Criterion) {
    c.bench_function("atomic/radial_wavefunction_100", |b| {
        b.iter(|| {
            for i in 0..100 {
                let r = i as f64 * 0.1;
                black_box(tanmatra::atomic::radial_wavefunction(1, 3, 2, r).ok());
            }
        });
    });
}

fn known_isotopes_alloc(c: &mut Criterion) {
    c.bench_function("decay/known_isotopes_alloc", |b| {
        b.iter(|| {
            black_box(tanmatra::decay::known_isotopes());
        });
    });
}

criterion_group!(
    benches,
    binding_energy_1000,
    spectral_line_1000,
    electron_config_36,
    decay_chain_10,
    shell_occupation_126,
    ionization_energy_118,
    lorentz_gamma_1000,
    rutherford_scattering_1000,
    bateman_chain_3,
    radial_wavefunction_100,
    known_isotopes_alloc,
);
criterion_main!(benches);
