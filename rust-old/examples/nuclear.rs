//! Nuclear physics: binding energy, shell model, decay chains, and fission.

use tanmatra::prelude::*;

fn main() {
    // --- Binding energy ---
    let fe56 = Nucleus::iron_56();
    println!("Fe-56 binding energy: {:.2} MeV", fe56.binding_energy());
    println!(
        "Fe-56 B/A: {:.3} MeV/nucleon",
        fe56.binding_energy_per_nucleon()
    );
    println!(
        "Fe-56 B/A (shell-corrected): {:.3} MeV/nucleon",
        fe56.binding_energy_shell_corrected() / 56.0
    );
    println!();

    // --- Shell model ---
    let (two_j, parity) = ground_state_spin_parity(&fe56);
    let parity_str = if parity > 0 { "+" } else { "-" };
    println!("Fe-56 spin-parity: {}/{}{parity_str}", two_j, 2);

    let o17 = Nucleus::new(8, 17).expect("valid nucleus");
    let (two_j, parity) = ground_state_spin_parity(&o17);
    let parity_str = if parity > 0 { "+" } else { "-" };
    println!("O-17 spin-parity: {}/2{parity_str}", two_j);
    println!();

    // --- Decay chain ---
    let u238 = Nucleus::uranium_238();
    let chain = decay_chain(&u238, 20);
    println!("U-238 decay chain ({} steps):", chain.len());
    for (nucleus, mode) in &chain {
        println!("  Z={}, A={} -> {:?}", nucleus.z(), nucleus.a(), mode);
    }
    println!();

    // --- Bateman equations ---
    // Simple 3-species chain: A -> B -> C (stable)
    let lambdas = [0.1, 0.05, 0.0];
    let pops = bateman_chain(&lambdas, 1e6, 50.0);
    println!("Bateman chain at t=50s:");
    for (i, n) in pops.iter().enumerate() {
        println!("  Species {}: {:.0} atoms", i + 1, n);
    }
}
