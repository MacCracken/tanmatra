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

fn binding_energy_shell_corrected_1000(c: &mut Criterion) {
    c.bench_function("nucleus/binding_energy_shell_corrected_1000", |b| {
        let nuclei: Vec<Nucleus> = (1..=100)
            .flat_map(|z| (z..=z + 9).filter_map(move |a| Nucleus::new(z, a).ok()))
            .collect();
        b.iter(|| {
            for n in &nuclei {
                black_box(n.binding_energy_shell_corrected());
            }
        });
    });
}

fn spin_parity_odd_odd_100(c: &mut Criterion) {
    c.bench_function("nucleus/ground_state_spin_parity_100", |b| {
        let nuclei: Vec<Nucleus> = (1..=99)
            .step_by(2)
            .filter_map(|z| Nucleus::new(z, 2 * z + 1).ok())
            .chain(
                (1..=99)
                    .step_by(2)
                    .filter_map(|z| Nucleus::new(z, 2 * z).ok()),
            )
            .collect();
        b.iter(|| {
            for n in &nuclei {
                black_box(ground_state_spin_parity(n));
            }
        });
    });
}

fn bateman_u238_series(c: &mut Criterion) {
    c.bench_function("decay/bateman_u238_series_15", |b| {
        let y = 365.2422 * 86400.0;
        let hl = [
            4.463e9 * y,
            24.107 * 86400.0,
            69.54,
            2.455e5 * y,
            7.54e4 * y,
            1600.0 * y,
            3.8215 * 86400.0,
            185.88,
            1623.6,
            1194.0,
            163.47e-6,
            22.2 * y,
            5.012 * 86400.0,
            138.376 * 86400.0,
        ];
        let mut lam: Vec<f64> = hl.iter().map(|h| core::f64::consts::LN_2 / h).collect();
        lam.push(0.0);
        b.iter(|| black_box(bateman_chain(black_box(&lam), 1.0, black_box(y))));
    });
}

fn einstein_a_n10(c: &mut Criterion) {
    c.bench_function("atomic/einstein_a_to_n2_upto_n10", |b| {
        b.iter(|| {
            for n in 3..=10 {
                black_box(einstein_a_coefficient(1, n, 1, 2, 0).ok());
                black_box(einstein_a_coefficient(1, n, 2, 2, 1).ok());
            }
        });
    });
}

fn radial_wavefunction_n10(c: &mut Criterion) {
    c.bench_function("atomic/radial_wavefunction_n10_100", |b| {
        b.iter(|| {
            for i in 0..100 {
                black_box(radial_wavefunction(1, 10, 3, f64::from(i) * 0.5).ok());
            }
        });
    });
}

fn lamb_shift_levels(c: &mut Criterion) {
    c.bench_function("atomic/lamb_shift_nlj_all_n4", |b| {
        b.iter(|| {
            for n in 1..=4 {
                for l in 0..n {
                    black_box(lamb_shift_nlj_ev(1, n, l, 2 * l + 1).ok());
                }
            }
        });
    });
}

fn pair_production_sweep(c: &mut Criterion) {
    c.bench_function("scattering/pair_production_1000", |b| {
        b.iter(|| {
            for i in 0..1000 {
                black_box(pair_production_cross_section(82, 1.1 + f64::from(i)));
            }
        });
    });
}

fn atomic_instant_add(c: &mut Criterion) {
    c.bench_function("timekeeping/atomic_instant_add_seconds_1000", |b| {
        let t = AtomicInstant::new(2_100_000_000, 123_456_789);
        b.iter(|| {
            let mut x = t;
            for _ in 0..1000 {
                x = x.add_seconds(black_box(1.000_000_001));
            }
            black_box(x)
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
    binding_energy_shell_corrected_1000,
    spin_parity_odd_odd_100,
    bateman_u238_series,
    einstein_a_n10,
    radial_wavefunction_n10,
    lamb_shift_levels,
    pair_production_sweep,
    atomic_instant_add,
);
criterion_main!(benches);
