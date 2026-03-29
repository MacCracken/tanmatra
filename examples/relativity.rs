//! Relativistic kinematics and scattering.

use tanmatra::prelude::*;

fn main() {
    // --- Lorentz factor ---
    println!("=== Relativistic Kinematics ===");
    for beta in [0.1, 0.5, 0.9, 0.99, 0.999] {
        let gamma = lorentz_gamma(beta);
        println!("  beta={beta:.3}: gamma={gamma:.4}");
    }
    println!();

    // --- Velocity addition ---
    let v_total = velocity_addition(0.9, 0.9);
    println!("Velocity addition: 0.9c + 0.9c = {v_total:.6}c");
    println!();

    // --- Four-momentum ---
    let proton = FourMomentum::from_mass_and_momentum(PROTON_MASS_MEV, 1000.0);
    println!("Proton at 1 GeV/c:");
    println!("  Energy: {:.1} MeV", proton.energy);
    println!("  Mass: {:.3} MeV", proton.invariant_mass());
    println!("  beta: {:.4}", proton.beta());
    println!("  gamma: {:.3}", proton.gamma());
    println!("  KE: {:.1} MeV", proton.kinetic_energy());
    println!();

    // --- de Broglie ---
    let lambda = de_broglie_wavelength_fm(1000.0);
    println!("de Broglie wavelength at 1 GeV/c: {lambda:.3} fm");
    println!();

    // --- Invariant mass ---
    let e_plus = FourMomentum::from_mass_and_momentum(ELECTRON_MASS_MEV, 500.0);
    let e_minus = FourMomentum::new(
        relativistic_energy(ELECTRON_MASS_MEV, 500.0),
        0.0,
        0.0,
        -500.0,
    );
    let m_inv = invariant_mass_two_body(&e_plus, &e_minus);
    println!("e+e- invariant mass (head-on, 500 MeV/c each): {m_inv:.1} MeV");
    println!();

    // --- Rutherford scattering ---
    println!("=== Rutherford Scattering ===");
    println!("5 MeV alpha on Au-197:");
    for angle_deg in [10, 30, 60, 90, 120, 150] {
        let theta = angle_deg as f64 * core::f64::consts::PI / 180.0;
        let ds = rutherford_differential(2, 79, 5.0, theta);
        println!("  theta={angle_deg}: dsigma/dOmega = {ds:.1} fm^2/sr");
    }
    let dist = distance_of_closest_approach(2, 79, 5.0);
    println!("  Distance of closest approach: {dist:.1} fm");
}
