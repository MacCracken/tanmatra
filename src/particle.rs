//! Standard Model particles: quarks, leptons, bosons, and fundamental forces.

extern crate alloc;
use alloc::borrow::Cow;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Quarks
// ---------------------------------------------------------------------------

/// The six quark flavors of the Standard Model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Quark {
    /// Up quark.
    Up,
    /// Down quark.
    Down,
    /// Charm quark.
    Charm,
    /// Strange quark.
    Strange,
    /// Top quark.
    Top,
    /// Bottom quark.
    Bottom,
}

/// Physical properties of a quark.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct QuarkProperties {
    /// Electric charge in thirds of e (+2 means +2/3 e, -1 means -1/3 e).
    pub charge_thirds: i8,
    /// Rest mass in MeV/c^2.
    pub mass_mev: f64,
    /// Whether the quark carries color charge.
    pub color: bool,
}

/// Returns the physical properties of a quark flavor.
///
/// Mass values from PDG 2024 review (MS-bar scheme at 2 `GeV` for light quarks).
#[must_use]
#[inline]
pub fn quark_properties(quark: Quark) -> QuarkProperties {
    match quark {
        Quark::Up => QuarkProperties {
            charge_thirds: 2,
            mass_mev: 2.16,
            color: true,
        },
        Quark::Down => QuarkProperties {
            charge_thirds: -1,
            mass_mev: 4.67,
            color: true,
        },
        Quark::Charm => QuarkProperties {
            charge_thirds: 2,
            mass_mev: 1270.0,
            color: true,
        },
        Quark::Strange => QuarkProperties {
            charge_thirds: -1,
            mass_mev: 93.4,
            color: true,
        },
        Quark::Top => QuarkProperties {
            charge_thirds: 2,
            mass_mev: 172_760.0,
            color: true,
        },
        Quark::Bottom => QuarkProperties {
            charge_thirds: -1,
            mass_mev: 4180.0,
            color: true,
        },
    }
}

// ---------------------------------------------------------------------------
// Leptons
// ---------------------------------------------------------------------------

/// The six leptons of the Standard Model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Lepton {
    /// Electron (e-).
    Electron,
    /// Muon (mu-).
    Muon,
    /// Tau (tau-).
    Tau,
    /// Electron neutrino.
    ElectronNeutrino,
    /// Muon neutrino.
    MuonNeutrino,
    /// Tau neutrino.
    TauNeutrino,
}

/// Physical properties of a lepton.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LeptonProperties {
    /// Electric charge in units of e (-1, 0).
    pub charge: i8,
    /// Rest mass in MeV/c^2.
    pub mass_mev: f64,
}

/// Returns the physical properties of a lepton.
///
/// Neutrino masses are listed as upper bounds; effectively zero for computation.
#[must_use]
#[inline]
pub fn lepton_properties(lepton: Lepton) -> LeptonProperties {
    match lepton {
        Lepton::Electron => LeptonProperties {
            charge: -1,
            mass_mev: 0.510_998_95,
        },
        Lepton::Muon => LeptonProperties {
            charge: -1,
            mass_mev: 105.658_375_5,
        },
        Lepton::Tau => LeptonProperties {
            charge: -1,
            mass_mev: 1776.86,
        },
        Lepton::ElectronNeutrino | Lepton::MuonNeutrino | Lepton::TauNeutrino => LeptonProperties {
            charge: 0,
            mass_mev: 0.0,
        },
    }
}

// ---------------------------------------------------------------------------
// Bosons
// ---------------------------------------------------------------------------

/// Gauge bosons and the Higgs boson.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum Boson {
    /// Photon -- electromagnetic force carrier.
    Photon,
    /// Gluon -- strong force carrier.
    Gluon,
    /// W+ boson -- weak force carrier (positive).
    WPlus,
    /// W- boson -- weak force carrier (negative).
    WMinus,
    /// Z boson -- weak force carrier (neutral).
    Z,
    /// Higgs boson -- mass generation.
    Higgs,
}

/// Physical properties of a boson.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BosonProperties {
    /// Rest mass in MeV/c^2 (0 for massless bosons).
    pub mass_mev: f64,
    /// Spin quantum number.
    pub spin: f64,
    /// Electric charge in units of e.
    pub charge: i8,
}

/// Returns the physical properties of a boson.
#[must_use]
#[inline]
pub fn boson_properties(boson: Boson) -> BosonProperties {
    match boson {
        Boson::Photon | Boson::Gluon => BosonProperties {
            mass_mev: 0.0,
            spin: 1.0,
            charge: 0,
        },
        Boson::WPlus => BosonProperties {
            mass_mev: 80_377.0,
            spin: 1.0,
            charge: 1,
        },
        Boson::WMinus => BosonProperties {
            mass_mev: 80_377.0,
            spin: 1.0,
            charge: -1,
        },
        Boson::Z => BosonProperties {
            mass_mev: 91_187.6,
            spin: 1.0,
            charge: 0,
        },
        Boson::Higgs => BosonProperties {
            mass_mev: 125_250.0,
            spin: 0.0,
            charge: 0,
        },
    }
}

// ---------------------------------------------------------------------------
// Fundamental Forces
// ---------------------------------------------------------------------------

/// The four fundamental forces of nature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum FundamentalForce {
    /// Strong nuclear force.
    Strong,
    /// Electromagnetic force.
    Electromagnetic,
    /// Weak nuclear force.
    Weak,
    /// Gravitational force.
    Gravitational,
}

/// Properties of a fundamental force.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ForceProperties<'a> {
    /// Relative strength (normalized to strong = 1).
    pub relative_strength: f64,
    /// Description of the force's range.
    pub range_description: Cow<'a, str>,
    /// Name of the mediator particle.
    pub mediator: Cow<'a, str>,
}

/// Returns the properties of a fundamental force.
#[must_use]
pub fn force_properties(force: FundamentalForce) -> ForceProperties<'static> {
    match force {
        FundamentalForce::Strong => ForceProperties {
            relative_strength: 1.0,
            range_description: Cow::Borrowed("~1 fm (confined to nucleus)"),
            mediator: Cow::Borrowed("gluon"),
        },
        FundamentalForce::Electromagnetic => ForceProperties {
            relative_strength: 1.0 / 137.0,
            range_description: Cow::Borrowed("infinite (1/r^2)"),
            mediator: Cow::Borrowed("photon"),
        },
        FundamentalForce::Weak => ForceProperties {
            relative_strength: 1e-6,
            range_description: Cow::Borrowed("~0.001 fm (W/Z boson range)"),
            mediator: Cow::Borrowed("W+, W-, Z bosons"),
        },
        FundamentalForce::Gravitational => ForceProperties {
            relative_strength: 6e-39,
            range_description: Cow::Borrowed("infinite (1/r^2)"),
            mediator: Cow::Borrowed("graviton (hypothetical)"),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quark_charges_sum_to_proton() {
        let u = quark_properties(Quark::Up);
        let d = quark_properties(Quark::Down);
        let total = u.charge_thirds + u.charge_thirds + d.charge_thirds;
        assert_eq!(total, 3);
    }

    #[test]
    fn quark_charges_sum_to_neutron() {
        let u = quark_properties(Quark::Up);
        let d = quark_properties(Quark::Down);
        let total = u.charge_thirds + d.charge_thirds + d.charge_thirds;
        assert_eq!(total, 0);
    }

    #[test]
    fn all_quarks_have_color() {
        let quarks = [
            Quark::Up,
            Quark::Down,
            Quark::Charm,
            Quark::Strange,
            Quark::Top,
            Quark::Bottom,
        ];
        for q in &quarks {
            assert!(quark_properties(*q).color);
        }
    }

    #[test]
    fn top_quark_heaviest() {
        let top = quark_properties(Quark::Top);
        let quarks = [
            Quark::Up,
            Quark::Down,
            Quark::Charm,
            Quark::Strange,
            Quark::Bottom,
        ];
        for q in &quarks {
            assert!(top.mass_mev > quark_properties(*q).mass_mev);
        }
    }

    #[test]
    fn electron_mass_matches_constant() {
        let e = lepton_properties(Lepton::Electron);
        let diff = libm::fabs(e.mass_mev - crate::constants::ELECTRON_MASS_MEV);
        assert!(diff < 1e-6);
    }

    #[test]
    fn neutrinos_neutral() {
        let neutrinos = [
            Lepton::ElectronNeutrino,
            Lepton::MuonNeutrino,
            Lepton::TauNeutrino,
        ];
        for n in &neutrinos {
            assert_eq!(lepton_properties(*n).charge, 0);
        }
    }

    #[test]
    fn photon_massless() {
        let p = boson_properties(Boson::Photon);
        assert!((p.mass_mev - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn w_boson_pair_charges() {
        let wp = boson_properties(Boson::WPlus);
        let wm = boson_properties(Boson::WMinus);
        assert_eq!(wp.charge, 1);
        assert_eq!(wm.charge, -1);
    }

    #[test]
    fn strong_force_strongest() {
        let forces = [
            FundamentalForce::Strong,
            FundamentalForce::Electromagnetic,
            FundamentalForce::Weak,
            FundamentalForce::Gravitational,
        ];
        let strong = force_properties(FundamentalForce::Strong);
        for f in &forces[1..] {
            assert!(strong.relative_strength > force_properties(*f).relative_strength);
        }
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
        let l = Lepton::Muon;
        let json = serde_json::to_string(&l).unwrap();
        let back: Lepton = serde_json::from_str(&json).unwrap();
        assert_eq!(l, back);
    }

    #[test]
    fn serde_roundtrip_boson() {
        let b = Boson::Higgs;
        let json = serde_json::to_string(&b).unwrap();
        let back: Boson = serde_json::from_str(&json).unwrap();
        assert_eq!(b, back);
    }

    #[test]
    fn serde_roundtrip_force() {
        let f = FundamentalForce::Weak;
        let json = serde_json::to_string(&f).unwrap();
        let back: FundamentalForce = serde_json::from_str(&json).unwrap();
        assert_eq!(f, back);
    }
}
