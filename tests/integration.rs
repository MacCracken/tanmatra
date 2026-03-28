use tanmatra::prelude::*;

#[test]
fn fe56_binding_energy_per_nucleon() {
    let fe56 = Nucleus::iron_56();
    let bea = fe56.binding_energy_per_nucleon();
    assert!(
        bea > 8.4 && bea < 9.2,
        "Fe-56 B/A = {bea}, expected ~8.8 MeV"
    );
}

#[test]
fn he4_binding_energy() {
    let he4 = Nucleus::helium_4();
    let b = he4.binding_energy();
    // Semi-empirical formula is approximate; real value is 28.3 MeV
    assert!(b > 15.0 && b < 35.0, "He-4 B = {b}, expected ~28.3 MeV");
}

#[test]
fn hydrogen_h_alpha() {
    let lambda = spectral_line_nm(1, 2, 3).unwrap();
    assert!(
        (lambda - 656.3).abs() < 1.0,
        "H-alpha = {lambda} nm, expected ~656.3 nm"
    );
}

#[test]
fn c14_half_life_5730_years() {
    let isotopes = known_isotopes();
    let c14 = isotopes.iter().find(|i| i.name == "C-14").unwrap();
    let years = c14.half_life_seconds / (365.25 * 24.0 * 3600.0);
    assert!(
        (years - 5730.0).abs() < 1.0,
        "C-14 half-life = {years} years"
    );
}

#[test]
fn c14_decay_constant_value() {
    let t_half = 5730.0 * 365.25 * 24.0 * 3600.0;
    let lambda = decay_constant(t_half).unwrap();
    let expected = 3.836e-12;
    let rel_err = (lambda - expected).abs() / expected;
    assert!(
        rel_err < 0.01,
        "C-14 lambda = {lambda}, expected ~{expected}"
    );
}

#[test]
fn alpha_decay_u238_to_th234() {
    let u238 = Nucleus::uranium_238();
    let daughter = alpha_decay(&u238).unwrap();
    assert_eq!(daughter.z(), 90, "Thorium Z");
    assert_eq!(daughter.a(), 234, "Thorium A");
}

#[test]
fn beta_minus_c14_to_n14() {
    let c14 = Nucleus::new(6, 14).unwrap();
    let daughter = beta_minus_decay(&c14).unwrap();
    assert_eq!(daughter.z(), 7, "Nitrogen Z");
    assert_eq!(daughter.a(), 14, "Nitrogen A");
}

#[test]
fn electron_config_iron() {
    // Fe (Z=26): [Ar] 3d6 4s2
    let config = electron_configuration(26).unwrap();
    let full = tanmatra::atomic::format_configuration(&config);
    assert_eq!(full, "1s2 2s2 2p6 3s2 3p6 4s2 3d6");
}

#[test]
fn electron_config_iron_short() {
    let config = electron_configuration(26).unwrap();
    let short = format_configuration_short(&config, 26);
    assert_eq!(short, "[Ar] 4s2 3d6");
}

#[test]
fn electron_config_chromium_exception() {
    // Cr (Z=24): [Ar] 3d5 4s1
    let config = electron_configuration(24).unwrap();
    let short = format_configuration_short(&config, 24);
    assert_eq!(short, "[Ar] 4s1 3d5");
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
    let n = Nucleus::iron_56();
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
    let qn = QuantumNumbers::new(3, 2, -1, 1).unwrap();
    let json = serde_json::to_string(&qn).unwrap();
    let back: QuantumNumbers = serde_json::from_str(&json).unwrap();
    assert_eq!(qn, back);
}

#[test]
fn u238_decay_chain_starts_alpha() {
    let u238 = Nucleus::uranium_238();
    let chain = decay_chain(&u238, 10);
    assert!(!chain.is_empty());
    // First decay should be alpha
    assert_eq!(chain[0].1, DecayMode::Alpha);
}

#[test]
fn ionization_energies_first_36() {
    // Spot-check a few NIST values
    assert!((ionization_energy_ev(1).unwrap() - 13.598).abs() < 0.01);
    assert!((ionization_energy_ev(2).unwrap() - 24.587).abs() < 0.01);
    assert!((ionization_energy_ev(10).unwrap() - 21.565).abs() < 0.01);
    assert!((ionization_energy_ev(26).unwrap() - 7.902).abs() < 0.01);
    assert!((ionization_energy_ev(36).unwrap() - 14.000).abs() < 0.01);
}
