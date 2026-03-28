//! Basic usage of tanmatra: binding energy, electron configuration, spectral lines.

use tanmatra::prelude::*;

fn main() {
    // --- Iron-56 binding energy ---
    let fe56 = Nucleus::iron_56();
    let be = fe56.binding_energy();
    let bea = fe56.binding_energy_per_nucleon();
    println!("Fe-56 binding energy: {be:.2} MeV");
    println!("Fe-56 B/A: {bea:.2} MeV/nucleon");
    println!("Fe-56 nuclear radius: {:.2} fm", fe56.nuclear_radius());
    println!("Fe-56 is magic: {}", fe56.is_magic());
    println!();

    // --- Electron configuration of iron ---
    let config = electron_configuration(26).unwrap();
    let full = format_configuration(&config);
    let short = format_configuration_short(&config, 26);
    println!("Fe electron config (full): {full}");
    println!("Fe electron config (short): {short}");
    println!(
        "Fe ionization energy: {:.3} eV",
        ionization_energy_ev(26).unwrap()
    );
    println!();

    // --- Chromium exception ---
    let cr_config = electron_configuration(24).unwrap();
    let cr_short = format_configuration_short(&cr_config, 24);
    println!("Cr electron config: {cr_short} (half-filled d shell)");
    println!();

    // --- H-alpha spectral line ---
    let h_alpha = spectral_line_nm(1, 2, 3).unwrap();
    println!("H-alpha wavelength: {h_alpha:.1} nm");

    let lyman_alpha = spectral_line_nm(1, 1, 2).unwrap();
    println!("Lyman-alpha wavelength: {lyman_alpha:.1} nm");
    println!();

    // --- Radioactive decay ---
    let isotopes = known_isotopes();
    let c14 = isotopes.iter().find(|i| i.name == "C-14").unwrap();
    let years = c14.half_life_seconds / (365.25 * 24.0 * 3600.0);
    println!("C-14 half-life: {years:.0} years");

    let frac = remaining_fraction(c14.half_life_seconds, 2.0 * c14.half_life_seconds);
    println!("C-14 remaining after 2 half-lives: {:.1}%", frac * 100.0);
    println!();

    // --- Nuclear reactions ---
    let dt = dt_fusion();
    println!("D-T fusion Q-value: {:.1} MeV", dt.q_value_mev);

    let d = Nucleus::new(1, 2).unwrap();
    let t = Nucleus::new(1, 3).unwrap();
    let barrier = coulomb_barrier(&d, &t).unwrap();
    println!("D-T Coulomb barrier: {barrier:.3} MeV");
    println!();

    // --- Standard Model ---
    println!("Top quark mass: {:.0} MeV", Quark::Top.mass_mev());
    println!("Higgs boson mass: {:.0} MeV", Boson::Higgs.mass_mev());
    println!(
        "Strong force relative strength: {:.0}",
        FundamentalForce::Strong.relative_strength()
    );
    println!(
        "Gravity relative strength: {:.1e}",
        FundamentalForce::Gravity.relative_strength()
    );
}
