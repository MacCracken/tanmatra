//! Spectral series, fine structure, and Zeeman/Stark effects.

use tanmatra::prelude::*;

fn main() {
    // --- Named spectral series ---
    println!("=== Hydrogen Spectral Series ===");
    println!();

    let lyman = lyman_series(1, 6).expect("valid series");
    println!("Lyman series (UV, n -> 1):");
    for (n, wavelength) in &lyman {
        println!("  n={n} -> 1: {wavelength:.2} nm");
    }
    println!();

    let balmer = balmer_series(1, 8).expect("valid series");
    println!("Balmer series (visible, n -> 2):");
    for (n, wavelength) in &balmer {
        println!("  n={n} -> 2: {wavelength:.2} nm");
    }
    println!();

    let paschen = paschen_series(1, 7).expect("valid series");
    println!("Paschen series (near-IR, n -> 3):");
    for (n, wavelength) in &paschen {
        println!("  n={n} -> 3: {wavelength:.1} nm");
    }
    println!();

    // --- Fine structure ---
    println!("=== Fine Structure ===");
    let e_2s = hydrogen_level_energy_ev(1, 2, 1).expect("valid");
    let e_2p32 = hydrogen_level_energy_ev(1, 2, 3).expect("valid");
    println!("H n=2 levels:");
    println!("  2S1/2: {e_2s:.6} eV");
    println!("  2P3/2: {e_2p32:.6} eV");
    println!("  Splitting: {:.3} meV", (e_2s - e_2p32).abs() * 1000.0);
    println!();

    // --- Lamb shift ---
    println!("=== QED Corrections ===");
    let lamb = lamb_shift_ev(1, 2, 0);
    println!(
        "H 2S Lamb shift: {:.3e} eV ({:.3} MHz)",
        lamb,
        lamb / 4.136e-9
    );
    println!();

    // --- Zeeman effect ---
    println!("=== Zeeman Effect (B = 1 T) ===");
    let g = lande_g_factor(0, 1); // 1s1/2
    println!("1s1/2 Lande g-factor: {g:.4}");
    for two_mj in [-1, 1] {
        let de = zeeman_splitting_ev(0, 1, two_mj, 1.0);
        println!("  mj={}/2: dE = {:.3e} eV", two_mj, de);
    }
    println!();

    // --- Stark effect ---
    println!("=== Stark Effect (E = 1e6 V/m) ===");
    let de = stark_shift_hydrogen_ev(2, 1, 1e6);
    println!("H n=2, k=1: dE = {:.3e} eV", de);
}
