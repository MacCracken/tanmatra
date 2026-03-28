use tanmatra::atomic::{OrbitalType, balmer_series, electron_configuration, ionization_energy_ev};
use tanmatra::constants::{FINE_STRUCTURE, PROTON_MASS_MEV};
use tanmatra::decay::{
    DecayMode, alpha_decay, beta_minus_decay, carbon14, decay_constant, uranium238,
};
use tanmatra::nucleus::{Nucleus, binding_energy_mev, binding_energy_per_nucleon, is_magic_number};
use tanmatra::particle::Quark;
use tanmatra::reaction::dt_fusion;

#[test]
fn fe56_binding_energy_per_nucleon() {
    let bea = binding_energy_per_nucleon(26, 56).unwrap();
    assert!(
        bea > 8.5 && bea < 9.1,
        "Fe-56 B/A = {bea}, expected ~8.8 MeV"
    );
}

#[test]
fn he4_binding_energy() {
    let b = binding_energy_mev(2, 4).unwrap();
    // Semi-empirical formula is approximate; real value is 28.3 MeV
    assert!(b > 24.0 && b < 32.0, "He-4 B = {b}, expected ~28.3 MeV");
}

#[test]
fn hydrogen_balmer_h_alpha() {
    let lambda = balmer_series(3).unwrap();
    assert!(
        (lambda - 656.3).abs() < 0.5,
        "H-alpha = {lambda} nm, expected ~656.3 nm"
    );
}

#[test]
fn c14_decay_constant() {
    let c14 = carbon14();
    let lambda = decay_constant(c14.half_life_s).unwrap();
    let expected = 3.836e-12;
    let rel_err = (lambda - expected).abs() / expected;
    assert!(
        rel_err < 0.01,
        "C-14 lambda = {lambda}, expected ~{expected}"
    );
}

#[test]
fn alpha_decay_u238_to_th234() {
    let daughter = alpha_decay(92, 238).unwrap();
    assert_eq!(daughter.z, 90, "Thorium Z");
    assert_eq!(daughter.a, 234, "Thorium A");
}

#[test]
fn beta_minus_c14_to_n14() {
    let daughter = beta_minus_decay(6, 14).unwrap();
    assert_eq!(daughter.z, 7, "Nitrogen Z");
    assert_eq!(daughter.a, 14, "Nitrogen A");
}

#[test]
fn electron_config_iron() {
    // Fe (Z=26): [Ar] 3d6 4s2
    let config = electron_configuration(26).unwrap();
    let d3 = config
        .iter()
        .find(|&&(n, ot, _)| n == 3 && ot == OrbitalType::D);
    let s4 = config
        .iter()
        .find(|&&(n, ot, _)| n == 4 && ot == OrbitalType::S);
    assert_eq!(d3, Some(&(3, OrbitalType::D, 6)));
    assert_eq!(s4, Some(&(4, OrbitalType::S, 2)));
}

#[test]
fn electron_config_chromium_exception() {
    // Cr (Z=24): [Ar] 3d5 4s1
    let config = electron_configuration(24).unwrap();
    let d3 = config
        .iter()
        .find(|&&(n, ot, _)| n == 3 && ot == OrbitalType::D);
    let s4 = config
        .iter()
        .find(|&&(n, ot, _)| n == 4 && ot == OrbitalType::S);
    assert_eq!(d3, Some(&(3, OrbitalType::D, 5)));
    assert_eq!(s4, Some(&(4, OrbitalType::S, 1)));
}

#[test]
fn dt_fusion_q_value() {
    let r = dt_fusion();
    assert!(
        (r.q_value_mev - 17.6).abs() < 0.1,
        "D-T Q = {}, expected ~17.6 MeV",
        r.q_value_mev
    );
}

#[test]
fn proton_mass() {
    assert!(
        (PROTON_MASS_MEV - 938.272).abs() < 0.001,
        "proton mass = {PROTON_MASS_MEV}"
    );
}

#[test]
fn fine_structure_value() {
    let expected = 1.0 / 137.036;
    assert!(
        (FINE_STRUCTURE - expected).abs() < 1e-6,
        "alpha = {FINE_STRUCTURE}"
    );
}

#[test]
fn magic_numbers_recognized() {
    for n in [2, 8, 20, 28, 50, 82, 126] {
        assert!(is_magic_number(n), "{n} should be a magic number");
    }
}

#[test]
fn serde_roundtrip_nucleus() {
    let n = Nucleus { z: 26, a: 56 };
    let json = serde_json::to_string(&n).unwrap();
    let back: Nucleus = serde_json::from_str(&json).unwrap();
    assert_eq!(n, back);
}

#[test]
fn serde_roundtrip_quark() {
    let q = Quark::Top;
    let json = serde_json::to_string(&q).unwrap();
    let back: Quark = serde_json::from_str(&json).unwrap();
    assert_eq!(q, back);
}

#[test]
fn serde_roundtrip_decay_mode() {
    let dm = DecayMode::Alpha;
    let json = serde_json::to_string(&dm).unwrap();
    let back: DecayMode = serde_json::from_str(&json).unwrap();
    assert_eq!(dm, back);
}

#[test]
fn serde_roundtrip_quantum_numbers() {
    let qn = tanmatra::atomic::QuantumNumbers {
        n: 3,
        l: 2,
        ml: -1,
        ms: 0.5,
    };
    let json = serde_json::to_string(&qn).unwrap();
    let back: tanmatra::atomic::QuantumNumbers = serde_json::from_str(&json).unwrap();
    assert_eq!(qn, back);
}

#[test]
fn u238_decay_chain_10_steps() {
    let u238 = uranium238();
    let chain = tanmatra::decay::decay_chain(&u238, 10).unwrap();
    assert_eq!(chain.len(), 10);
    // After 10 alpha decays: Z = 92 - 20 = 72, A = 238 - 40 = 198
    assert_eq!(chain[9].z, 72);
    assert_eq!(chain[9].a, 198);
}

#[test]
fn ionization_energies_first_36() {
    // Spot-check a few NIST values
    assert!((ionization_energy_ev(1).unwrap() - 13.598).abs() < 0.001);
    assert!((ionization_energy_ev(2).unwrap() - 24.587).abs() < 0.001);
    assert!((ionization_energy_ev(10).unwrap() - 21.565).abs() < 0.001);
    assert!((ionization_energy_ev(26).unwrap() - 7.902).abs() < 0.001);
    assert!((ionization_energy_ev(36).unwrap() - 14.000).abs() < 0.001);
}
