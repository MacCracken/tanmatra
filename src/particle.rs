//! Standard Model particles and fundamental forces.
//!
//! Contains the three families of fundamental particles (quarks, leptons, bosons)
//! and the four fundamental forces with their relative strengths and ranges.
//!
//! All masses are from the Particle Data Group (PDG) 2024 Review.

use serde::{Deserialize, Serialize};

/// Quark flavors of the Standard Model.
///
/// Masses are current quark masses (MSbar scheme) from PDG 2024,
/// except top which is the pole mass.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Quark {
    /// Up quark (charge +2/3, mass ~2.16 MeV).
    Up,
    /// Down quark (charge -1/3, mass ~4.67 MeV).
    Down,
    /// Charm quark (charge +2/3, mass ~1270 MeV).
    Charm,
    /// Strange quark (charge -1/3, mass ~93.4 MeV).
    Strange,
    /// Top quark (charge +2/3, mass ~172570 MeV).
    Top,
    /// Bottom quark (charge -1/3, mass ~4180 MeV).
    Bottom,
}

impl Quark {
    /// Returns the quark mass in MeV/c^2 (PDG 2024).
    #[must_use]
    pub const fn mass_mev(self) -> f64 {
        match self {
            Self::Up => 2.16,
            Self::Down => 4.67,
            Self::Charm => 1_270.0,
            Self::Strange => 93.4,
            Self::Top => 172_570.0,
            Self::Bottom => 4_180.0,
        }
    }

    /// Returns the electric charge in units of elementary charge.
    #[must_use]
    pub const fn charge(self) -> f64 {
        match self {
            Self::Up | Self::Charm | Self::Top => 2.0 / 3.0,
            Self::Down | Self::Strange | Self::Bottom => -1.0 / 3.0,
        }
    }

    /// Returns the generation (1, 2, or 3).
    #[must_use]
    pub const fn generation(self) -> u8 {
        match self {
            Self::Up | Self::Down => 1,
            Self::Charm | Self::Strange => 2,
            Self::Top | Self::Bottom => 3,
        }
    }
}

/// Lepton flavors of the Standard Model.
///
/// Masses from CODATA 2022 / PDG 2024.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Lepton {
    /// Electron (charge -1, mass 0.51099895 MeV).
    Electron,
    /// Electron neutrino (charge 0, mass < 0.8 eV).
    ElectronNeutrino,
    /// Muon (charge -1, mass 105.6583755 MeV).
    Muon,
    /// Muon neutrino (charge 0, mass < 0.19 MeV).
    MuonNeutrino,
    /// Tau (charge -1, mass 1776.93 MeV).
    Tau,
    /// Tau neutrino (charge 0, mass < 18.2 MeV).
    TauNeutrino,
}

impl Lepton {
    /// Returns the lepton mass in MeV/c^2.
    ///
    /// Neutrino masses are upper bounds; we return 0.0 for computations.
    #[must_use]
    pub const fn mass_mev(self) -> f64 {
        match self {
            Self::Electron => 0.510_998_950,
            Self::ElectronNeutrino => 0.0,
            Self::Muon => 105.658_375_5,
            Self::MuonNeutrino => 0.0,
            Self::Tau => 1_776.93,
            Self::TauNeutrino => 0.0,
        }
    }

    /// Returns the electric charge in units of elementary charge.
    #[must_use]
    pub const fn charge(self) -> f64 {
        match self {
            Self::Electron | Self::Muon | Self::Tau => -1.0,
            Self::ElectronNeutrino | Self::MuonNeutrino | Self::TauNeutrino => 0.0,
        }
    }

    /// Returns the generation (1, 2, or 3).
    #[must_use]
    pub const fn generation(self) -> u8 {
        match self {
            Self::Electron | Self::ElectronNeutrino => 1,
            Self::Muon | Self::MuonNeutrino => 2,
            Self::Tau | Self::TauNeutrino => 3,
        }
    }

    /// Returns the mean lifetime in seconds (PDG 2024).
    ///
    /// Returns [`f64::INFINITY`] for stable particles (electron, neutrinos).
    #[must_use]
    #[inline]
    pub const fn lifetime_seconds(self) -> f64 {
        match self {
            Self::Muon => 2.196_981_1e-6, // PDG 2024
            Self::Tau => 2.903e-13,       // PDG 2024
            Self::Electron | Self::ElectronNeutrino | Self::MuonNeutrino | Self::TauNeutrino => {
                f64::INFINITY
            }
        }
    }

    /// Returns the total decay width in MeV (PDG 2024).
    ///
    /// Derived from Γ = ħ / τ. Returns 0.0 for stable particles.
    #[must_use]
    #[inline]
    pub fn decay_width_mev(self) -> f64 {
        let tau = self.lifetime_seconds();
        if tau.is_infinite() {
            return 0.0;
        }
        crate::constants::HBAR_MEV_S / tau
    }

    /// Returns `true` if this lepton is stable (infinite lifetime).
    #[must_use]
    pub const fn is_stable(self) -> bool {
        matches!(
            self,
            Self::Electron | Self::ElectronNeutrino | Self::MuonNeutrino | Self::TauNeutrino
        )
    }
}

/// Gauge bosons and the Higgs boson of the Standard Model.
///
/// Masses from PDG 2024.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Boson {
    /// Photon (massless, mediates electromagnetic force).
    Photon,
    /// Gluon (massless, mediates strong force).
    Gluon,
    /// W+ boson (mass 80369 MeV).
    WPlus,
    /// W- boson (mass 80369 MeV).
    WMinus,
    /// Z boson (mass 91188 MeV).
    Z,
    /// Higgs boson (mass 125250 MeV).
    Higgs,
}

impl Boson {
    /// Returns the boson mass in MeV/c^2 (PDG 2024).
    #[must_use]
    pub const fn mass_mev(self) -> f64 {
        match self {
            Self::Photon | Self::Gluon => 0.0,
            Self::WPlus | Self::WMinus => 80_369.0,
            Self::Z => 91_188.0,
            Self::Higgs => 125_250.0,
        }
    }

    /// Returns the spin quantum number.
    #[must_use]
    pub const fn spin(self) -> u8 {
        match self {
            Self::Higgs => 0,
            Self::Photon | Self::Gluon | Self::WPlus | Self::WMinus | Self::Z => 1,
        }
    }

    /// Returns the electric charge in units of elementary charge.
    #[must_use]
    pub const fn charge(self) -> f64 {
        match self {
            Self::Photon | Self::Gluon | Self::Z | Self::Higgs => 0.0,
            Self::WPlus => 1.0,
            Self::WMinus => -1.0,
        }
    }

    /// Returns the total decay width in MeV (PDG 2024).
    ///
    /// Returns 0.0 for stable particles (photon) and confined particles (gluon).
    #[must_use]
    #[inline]
    pub const fn decay_width_mev(self) -> f64 {
        match self {
            Self::Photon | Self::Gluon => 0.0,
            Self::WPlus | Self::WMinus => 2_085.0, // 2.085 GeV
            Self::Z => 2_495.5,                    // 2.4955 GeV
            Self::Higgs => 3.7,                    // 3.7 MeV (CMS off-shell)
        }
    }

    /// Returns the mean lifetime in seconds (PDG 2024).
    ///
    /// Returns [`f64::INFINITY`] for stable (photon) or confined (gluon) particles.
    #[must_use]
    #[inline]
    pub fn lifetime_seconds(self) -> f64 {
        let width = self.decay_width_mev();
        if width <= 0.0 {
            return f64::INFINITY;
        }
        crate::constants::HBAR_MEV_S / width
    }

    /// Returns `true` if this boson is stable (infinite lifetime).
    #[must_use]
    pub const fn is_stable(self) -> bool {
        matches!(self, Self::Photon | Self::Gluon)
    }
}

/// The four fundamental forces of nature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum FundamentalForce {
    /// Strong nuclear force (relative strength ~1, range ~1e-15 m).
    Strong,
    /// Electromagnetic force (relative strength ~1/137, infinite range).
    Electromagnetic,
    /// Weak nuclear force (relative strength ~1e-6, range ~1e-18 m).
    Weak,
    /// Gravitational force (relative strength ~6e-39, infinite range).
    Gravity,
}

impl FundamentalForce {
    /// Returns the approximate relative coupling strength compared to the
    /// strong force (which is normalized to 1).
    #[must_use]
    pub const fn relative_strength(self) -> f64 {
        match self {
            Self::Strong => 1.0,
            Self::Electromagnetic => 7.297_353e-3, // ~1/137
            Self::Weak => 1.0e-6,
            Self::Gravity => 6.0e-39,
        }
    }

    /// Returns the approximate range of the force in meters.
    ///
    /// Returns [`f64::INFINITY`] for infinite-range forces.
    #[must_use]
    pub const fn range_meters(self) -> f64 {
        match self {
            Self::Strong => 1.0e-15,
            Self::Electromagnetic => f64::INFINITY,
            Self::Weak => 1.0e-18,
            Self::Gravity => f64::INFINITY,
        }
    }

    /// Returns the mediating boson(s) for this force.
    #[must_use]
    pub fn mediator(self) -> &'static [Boson] {
        match self {
            Self::Strong => &[Boson::Gluon],
            Self::Electromagnetic => &[Boson::Photon],
            Self::Weak => &[Boson::WPlus, Boson::WMinus, Boson::Z],
            Self::Gravity => &[], // graviton not yet observed
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quark_masses_ordered() {
        assert!(Quark::Up.mass_mev() < Quark::Down.mass_mev());
        assert!(Quark::Down.mass_mev() < Quark::Strange.mass_mev());
        assert!(Quark::Strange.mass_mev() < Quark::Charm.mass_mev());
        assert!(Quark::Charm.mass_mev() < Quark::Bottom.mass_mev());
        assert!(Quark::Bottom.mass_mev() < Quark::Top.mass_mev());
    }

    #[test]
    fn lepton_masses_ordered() {
        assert!(Lepton::Electron.mass_mev() < Lepton::Muon.mass_mev());
        assert!(Lepton::Muon.mass_mev() < Lepton::Tau.mass_mev());
    }

    #[test]
    fn force_strengths_ordered() {
        assert!(
            FundamentalForce::Strong.relative_strength()
                > FundamentalForce::Electromagnetic.relative_strength()
        );
        assert!(
            FundamentalForce::Electromagnetic.relative_strength()
                > FundamentalForce::Weak.relative_strength()
        );
        assert!(
            FundamentalForce::Weak.relative_strength()
                > FundamentalForce::Gravity.relative_strength()
        );
    }

    #[test]
    fn quark_charges() {
        let eps = 1e-10;
        assert!((Quark::Up.charge() - 2.0 / 3.0).abs() < eps);
        assert!((Quark::Down.charge() - (-1.0 / 3.0)).abs() < eps);
    }

    #[test]
    fn boson_higgs_mass() {
        assert!((Boson::Higgs.mass_mev() - 125_250.0).abs() < 1.0);
    }

    #[test]
    fn serde_roundtrip_quark() {
        let q = Quark::Charm;
        let json = serde_json::to_string(&q).unwrap();
        let back: Quark = serde_json::from_str(&json).unwrap();
        assert_eq!(q, back);
    }

    #[test]
    fn serde_roundtrip_lepton() {
        let l = Lepton::Tau;
        let json = serde_json::to_string(&l).unwrap();
        let back: Lepton = serde_json::from_str(&json).unwrap();
        assert_eq!(l, back);
    }

    #[test]
    fn serde_roundtrip_boson() {
        let b = Boson::Z;
        let json = serde_json::to_string(&b).unwrap();
        let back: Boson = serde_json::from_str(&json).unwrap();
        assert_eq!(b, back);
    }

    #[test]
    fn serde_roundtrip_force() {
        let f = FundamentalForce::Strong;
        let json = serde_json::to_string(&f).unwrap();
        let back: FundamentalForce = serde_json::from_str(&json).unwrap();
        assert_eq!(f, back);
    }

    #[test]
    fn em_mediator_is_photon() {
        assert_eq!(
            FundamentalForce::Electromagnetic.mediator(),
            &[Boson::Photon]
        );
    }

    #[test]
    fn muon_lifetime() {
        let tau = Lepton::Muon.lifetime_seconds();
        assert!((tau - 2.197e-6).abs() < 1e-8, "muon lifetime={tau}");
    }

    #[test]
    fn tau_lifetime() {
        let tau = Lepton::Tau.lifetime_seconds();
        assert!((tau - 2.903e-13).abs() < 1e-15, "tau lifetime={tau}");
    }

    #[test]
    fn electron_is_stable() {
        assert!(Lepton::Electron.is_stable());
        assert!(Lepton::Electron.lifetime_seconds().is_infinite());
        assert!((Lepton::Electron.decay_width_mev()).abs() < 1e-30);
    }

    #[test]
    fn neutrinos_stable() {
        assert!(Lepton::ElectronNeutrino.is_stable());
        assert!(Lepton::MuonNeutrino.is_stable());
        assert!(Lepton::TauNeutrino.is_stable());
    }

    #[test]
    fn w_boson_decay_width() {
        let gamma = Boson::WPlus.decay_width_mev();
        // W width ≈ 2085 MeV (2.085 GeV)
        assert!((gamma - 2085.0).abs() < 50.0, "W width={gamma} MeV");
    }

    #[test]
    fn z_boson_decay_width() {
        let gamma = Boson::Z.decay_width_mev();
        // Z width ≈ 2495.5 MeV (2.4955 GeV)
        assert!((gamma - 2495.5).abs() < 5.0, "Z width={gamma} MeV");
    }

    #[test]
    fn photon_is_stable() {
        assert!(Boson::Photon.is_stable());
        assert!(Boson::Photon.lifetime_seconds().is_infinite());
    }

    #[test]
    fn width_lifetime_consistency() {
        // Verify Γ·τ = ħ for the muon
        let gamma = Lepton::Muon.decay_width_mev();
        let tau = Lepton::Muon.lifetime_seconds();
        let hbar = crate::constants::HBAR_MEV_S;
        let product = gamma * tau;
        let rel_err = (product - hbar).abs() / hbar;
        assert!(rel_err < 1e-6, "Γ·τ={product}, ħ={hbar}");
    }

    // --- Coverage: exercise all enum variants ---

    #[test]
    fn quark_generations() {
        assert_eq!(Quark::Up.generation(), 1);
        assert_eq!(Quark::Down.generation(), 1);
        assert_eq!(Quark::Charm.generation(), 2);
        assert_eq!(Quark::Strange.generation(), 2);
        assert_eq!(Quark::Top.generation(), 3);
        assert_eq!(Quark::Bottom.generation(), 3);
    }

    #[test]
    fn lepton_all_masses() {
        assert!(Lepton::Electron.mass_mev() > 0.0);
        assert!((Lepton::ElectronNeutrino.mass_mev()).abs() < 1e-10);
        assert!(Lepton::Muon.mass_mev() > 100.0);
        assert!((Lepton::MuonNeutrino.mass_mev()).abs() < 1e-10);
        assert!(Lepton::Tau.mass_mev() > 1700.0);
        assert!((Lepton::TauNeutrino.mass_mev()).abs() < 1e-10);
    }

    #[test]
    fn lepton_all_charges() {
        assert!((Lepton::Electron.charge() - (-1.0)).abs() < 1e-10);
        assert!((Lepton::Muon.charge() - (-1.0)).abs() < 1e-10);
        assert!((Lepton::Tau.charge() - (-1.0)).abs() < 1e-10);
        assert!((Lepton::ElectronNeutrino.charge()).abs() < 1e-10);
        assert!((Lepton::MuonNeutrino.charge()).abs() < 1e-10);
        assert!((Lepton::TauNeutrino.charge()).abs() < 1e-10);
    }

    #[test]
    fn lepton_all_generations() {
        assert_eq!(Lepton::Electron.generation(), 1);
        assert_eq!(Lepton::ElectronNeutrino.generation(), 1);
        assert_eq!(Lepton::Muon.generation(), 2);
        assert_eq!(Lepton::MuonNeutrino.generation(), 2);
        assert_eq!(Lepton::Tau.generation(), 3);
        assert_eq!(Lepton::TauNeutrino.generation(), 3);
    }

    #[test]
    fn boson_all_properties() {
        assert!((Boson::Photon.mass_mev()).abs() < 1e-10);
        assert!((Boson::Gluon.mass_mev()).abs() < 1e-10);
        assert!(Boson::WPlus.mass_mev() > 80_000.0);
        assert!(Boson::WMinus.mass_mev() > 80_000.0);
        assert!(Boson::Z.mass_mev() > 91_000.0);
        assert_eq!(Boson::Higgs.spin(), 0);
        assert_eq!(Boson::Photon.spin(), 1);
        assert!((Boson::WPlus.charge() - 1.0).abs() < 1e-10);
        assert!((Boson::WMinus.charge() - (-1.0)).abs() < 1e-10);
        assert!((Boson::Z.charge()).abs() < 1e-10);
        assert!((Boson::Higgs.charge()).abs() < 1e-10);
    }

    #[test]
    fn force_all_ranges() {
        assert!(FundamentalForce::Strong.range_meters() < 1e-14);
        assert!(
            FundamentalForce::Electromagnetic
                .range_meters()
                .is_infinite()
        );
        assert!(FundamentalForce::Weak.range_meters() < 1e-17);
        assert!(FundamentalForce::Gravity.range_meters().is_infinite());
    }

    #[test]
    fn force_all_mediators() {
        assert!(!FundamentalForce::Strong.mediator().is_empty());
        assert!(!FundamentalForce::Electromagnetic.mediator().is_empty());
        assert!(!FundamentalForce::Weak.mediator().is_empty());
        assert!(FundamentalForce::Gravity.mediator().is_empty());
    }

    #[test]
    fn boson_gluon_is_stable() {
        assert!(Boson::Gluon.is_stable());
    }

    #[test]
    fn tau_not_stable() {
        assert!(!Lepton::Tau.is_stable());
        assert!(!Lepton::Muon.is_stable());
    }

    #[test]
    fn higgs_decay_width() {
        assert!((Boson::Higgs.decay_width_mev() - 3.7).abs() < 0.1);
    }
}
