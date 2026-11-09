//! `relativity.tsv` — port plan M1 grid.
//!
//! - β ∈ {0, ±0.5, ±(1 − 1e-12), ±1, ±1.5}, widened with more points near 0 and 1, ±∞
//!   and NaN.
//! - Four-momenta from 1 keV to 10 GeV for every mass in [`masses`], along z and along an
//!   oblique direction, plus massless, spacelike, negative-energy and NaN edges.
//!
//! Classes: ±, ×, ÷ and `sqrt` only → B. Anything through `momentum_magnitude`, which is
//! `libm::hypot`, → C. Constructors that copy their inputs → A.

use crate::grids::log_grid;
use crate::tsv::{Class, Module, Tok, f};
use tanmatra::constants::{AMU_MEV, C, ELECTRON_MASS_MEV, NEUTRON_MASS_MEV, PROTON_MASS_MEV};
use tanmatra::particle::{Boson, Lepton};
use tanmatra::relativity::{
    FourMomentum, beta_to_velocity, de_broglie_wavelength_fm, gamma_to_beta,
    invariant_mass_two_body, lorentz_gamma, relativistic_energy, relativistic_momentum,
    velocity_addition, velocity_to_beta,
};

fn betas() -> Vec<f64> {
    let mut v = vec![0.0, -0.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN];
    for x in [
        1e-12,
        1e-8,
        0.1,
        0.5,
        0.9,
        0.99,
        0.999_999,
        1.0 - 1e-12,
        1.0,
        1.5,
    ] {
        v.push(x);
        v.push(-x);
    }
    v
}

fn gammas() -> Vec<f64> {
    vec![
        -1.0,
        0.0,
        0.5,
        1.0,
        1.0 + 1e-12,
        1.0 + 1e-8,
        1.0001,
        1.5,
        2.0,
        10.0,
        1e3,
        1e6,
        1e15,
        1e300,
        f64::INFINITY,
        f64::NAN,
    ]
}

/// Rest masses (MeV): massless, the tabulated particles, and a negative edge.
fn masses() -> Vec<f64> {
    vec![
        0.0,
        ELECTRON_MASS_MEV,
        Lepton::Muon.mass_mev(),
        PROTON_MASS_MEV,
        NEUTRON_MASS_MEV,
        AMU_MEV,
        Lepton::Tau.mass_mev(),
        Boson::Higgs.mass_mev(),
        -1.0,
    ]
}

/// Momenta / kinetic energies (MeV): 1 keV to 10 GeV at 4 per decade, plus edges.
fn momenta() -> Vec<f64> {
    let mut v = vec![0.0, 1e-9];
    v.extend(log_grid(1e-3, 7, 4));
    v.extend([1e8, -1.0, f64::INFINITY, f64::NAN]);
    v
}

fn fm(p: &FourMomentum) -> [Tok; 4] {
    [f(p.energy), f(p.px), f(p.py), f(p.pz)]
}

/// The four-momenta the methods are probed on.
fn four_momenta() -> Vec<FourMomentum> {
    let mut v = Vec::new();
    for m in masses() {
        for p in momenta() {
            v.push(FourMomentum::from_mass_and_momentum(m, p));
            // Oblique: (3, 4, 12) / 13 is a unit vector with exact rational components.
            let e = relativistic_energy(m, p);
            v.push(FourMomentum::new(
                e,
                3.0 * p / 13.0,
                4.0 * p / 13.0,
                12.0 * p / 13.0,
            ));
        }
    }
    v.extend([
        FourMomentum::new(0.0, 0.0, 0.0, 0.0),
        FourMomentum::new(5.0, 5.0, 0.0, 0.0), // massless along x
        FourMomentum::new(5.0, 0.0, 0.0, -5.0), // massless along −z
        FourMomentum::new(1.0, 2.0, 0.0, 0.0), // spacelike
        FourMomentum::new(-1.0, 0.0, 0.0, 0.0), // negative energy
        FourMomentum::new(f64::NAN, 0.0, 0.0, 0.0),
        FourMomentum::new(1.0, f64::NAN, 0.0, 0.0),
        FourMomentum::new(f64::INFINITY, 1.0, 0.0, 0.0),
    ]);
    v
}

pub fn emit() -> Module {
    let mut m = Module::new("relativity");

    for &beta in &betas() {
        m.row(
            "lorentz_gamma",
            Class::B,
            &[f(beta)],
            &[f(lorentz_gamma(beta))],
        );
        m.row(
            "beta_to_velocity",
            Class::B,
            &[f(beta)],
            &[f(beta_to_velocity(beta))],
        );
        let v = beta * C;
        m.row(
            "velocity_to_beta",
            Class::B,
            &[f(v)],
            &[f(velocity_to_beta(v))],
        );
    }
    for v in [1.0, 1e3, 3e7] {
        m.row(
            "velocity_to_beta",
            Class::B,
            &[f(v)],
            &[f(velocity_to_beta(v))],
        );
    }
    for &g in &gammas() {
        m.row("gamma_to_beta", Class::B, &[f(g)], &[f(gamma_to_beta(g))]);
    }
    for &b1 in &betas() {
        for &b2 in &betas() {
            let out = velocity_addition(b1, b2);
            m.row("velocity_addition", Class::B, &[f(b1), f(b2)], &[f(out)]);
        }
    }
    for &mass in &masses() {
        for &p in &momenta() {
            let e = relativistic_energy(mass, p);
            m.row("relativistic_energy", Class::B, &[f(mass), f(p)], &[f(e)]);
            let q = relativistic_momentum(mass, p);
            m.row("relativistic_momentum", Class::B, &[f(mass), f(p)], &[f(q)]);
            let fmp = FourMomentum::from_mass_and_momentum(mass, p);
            m.row(
                "FourMomentum::from_mass_and_momentum",
                Class::B,
                &[f(mass), f(p)],
                &fm(&fmp),
            );
        }
        let rest = FourMomentum::at_rest(mass);
        m.row("FourMomentum::at_rest", Class::A, &[f(mass)], &fm(&rest));
    }
    for &p in &momenta() {
        m.row(
            "de_broglie_wavelength_fm",
            Class::B,
            &[f(p)],
            &[f(de_broglie_wavelength_fm(p))],
        );
    }

    let probes = four_momenta();
    for p in &probes {
        let ins = fm(p);
        let made = FourMomentum::new(p.energy, p.px, p.py, p.pz);
        m.row("FourMomentum::new", Class::A, &ins, &fm(&made));
        m.row(
            "FourMomentum::momentum_magnitude",
            Class::C,
            &ins,
            &[f(p.momentum_magnitude())],
        );
        m.row(
            "FourMomentum::invariant_mass",
            Class::C,
            &ins,
            &[f(p.invariant_mass())],
        );
        m.row("FourMomentum::beta", Class::C, &ins, &[f(p.beta())]);
        m.row("FourMomentum::gamma", Class::C, &ins, &[f(p.gamma())]);
        m.row(
            "FourMomentum::kinetic_energy",
            Class::C,
            &ins,
            &[f(p.kinetic_energy())],
        );
    }

    // Pairs: every 29th probe against every 23rd (22 × 28 = 616 pairs) samples every mass,
    // momentum and edge without the full 638² product.
    for p1 in probes.iter().step_by(29) {
        for p2 in probes.iter().step_by(23) {
            let mut ins = fm(p1).to_vec();
            ins.extend(fm(p2));
            m.row("FourMomentum::add", Class::B, &ins, &fm(&(*p1 + *p2)));
            m.row(
                "invariant_mass_two_body",
                Class::C,
                &ins,
                &[f(invariant_mass_two_body(p1, p2))],
            );
        }
    }

    m
}
