//! `particle.tsv` — every property of every Standard Model variant (port plan M1).
//!
//! Classes: tabulated values are A. `Lepton::decay_width_mev` and
//! `Boson::lifetime_seconds` compute ħ / τ and ħ / Γ, so they are B.

use crate::tsv::{Class, Module, Tok, b, e, f, i};
use tanmatra::particle::{Boson, FundamentalForce, Lepton, Quark};

pub const QUARKS: [Quark; 6] = [
    Quark::Up,
    Quark::Down,
    Quark::Charm,
    Quark::Strange,
    Quark::Top,
    Quark::Bottom,
];

pub const LEPTONS: [Lepton; 6] = [
    Lepton::Electron,
    Lepton::ElectronNeutrino,
    Lepton::Muon,
    Lepton::MuonNeutrino,
    Lepton::Tau,
    Lepton::TauNeutrino,
];

pub const BOSONS: [Boson; 6] = [
    Boson::Photon,
    Boson::Gluon,
    Boson::WPlus,
    Boson::WMinus,
    Boson::Z,
    Boson::Higgs,
];

pub const FORCES: [FundamentalForce; 4] = [
    FundamentalForce::Strong,
    FundamentalForce::Electromagnetic,
    FundamentalForce::Weak,
    FundamentalForce::Gravity,
];

pub fn emit() -> Module {
    let mut m = Module::new("particle");

    for q in QUARKS {
        let v = [e("Quark", q)];
        m.row("Quark::mass_mev", Class::A, &v, &[f(q.mass_mev())]);
        m.row("Quark::charge", Class::A, &v, &[f(q.charge())]);
        m.row("Quark::generation", Class::A, &v, &[i(q.generation())]);
    }

    for l in LEPTONS {
        let v = [e("Lepton", l)];
        m.row("Lepton::mass_mev", Class::A, &v, &[f(l.mass_mev())]);
        m.row("Lepton::charge", Class::A, &v, &[f(l.charge())]);
        m.row("Lepton::generation", Class::A, &v, &[i(l.generation())]);
        m.row(
            "Lepton::lifetime_seconds",
            Class::A,
            &v,
            &[f(l.lifetime_seconds())],
        );
        m.row(
            "Lepton::decay_width_mev",
            Class::B,
            &v,
            &[f(l.decay_width_mev())],
        );
        m.row("Lepton::is_stable", Class::A, &v, &[b(l.is_stable())]);
    }

    for bo in BOSONS {
        let v = [e("Boson", bo)];
        m.row("Boson::mass_mev", Class::A, &v, &[f(bo.mass_mev())]);
        m.row("Boson::spin", Class::A, &v, &[i(bo.spin())]);
        m.row("Boson::charge", Class::A, &v, &[f(bo.charge())]);
        m.row(
            "Boson::decay_width_mev",
            Class::A,
            &v,
            &[f(bo.decay_width_mev())],
        );
        m.row(
            "Boson::lifetime_seconds",
            Class::B,
            &v,
            &[f(bo.lifetime_seconds())],
        );
        m.row("Boson::is_stable", Class::A, &v, &[b(bo.is_stable())]);
    }

    for force in FORCES {
        let v = [e("FundamentalForce", force)];
        m.row(
            "FundamentalForce::relative_strength",
            Class::A,
            &v,
            &[f(force.relative_strength())],
        );
        m.row(
            "FundamentalForce::range_meters",
            Class::A,
            &v,
            &[f(force.range_meters())],
        );
        let med = force.mediator();
        let mut out = vec![Tok::V(med.len())];
        out.extend(med.iter().map(|&x| e("Boson", x)));
        m.row("FundamentalForce::mediator", Class::A, &v, &out);
    }

    m
}
