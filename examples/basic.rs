//! Basic tanmatra usage: binding energy, electron configuration, spectral lines.

use tanmatra::atomic::{OrbitalType, balmer_series, electron_configuration};
use tanmatra::nucleus::{IRON56, binding_energy_per_nucleon};

fn main() {
    // Iron-56: peak of the binding energy curve
    let bea = binding_energy_per_nucleon(IRON56.z, IRON56.a).unwrap();
    println!("Fe-56 binding energy per nucleon: {bea:.3} MeV");

    // Electron configuration of iron (Z=26)
    let config = electron_configuration(26).unwrap();
    print!("Fe electron configuration: ");
    for (n, ot, count) in &config {
        let label = match ot {
            OrbitalType::S => "s",
            OrbitalType::P => "p",
            OrbitalType::D => "d",
            OrbitalType::F => "f",
            _ => "?",
        };
        print!("{n}{label}{count} ");
    }
    println!();

    // Hydrogen H-alpha line (Balmer series, n=3 -> n=2)
    let h_alpha = balmer_series(3).unwrap();
    println!("Hydrogen H-alpha wavelength: {h_alpha:.1} nm");
}
