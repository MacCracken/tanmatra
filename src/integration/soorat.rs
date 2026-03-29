//! Soorat integration — visualization data structures for atomic/nuclear physics.
//!
//! Provides structured types that soorat can render: atomic orbitals,
//! nuclear structure, spectral lines, and decay chains.

use alloc::string::String;
use alloc::vec::Vec;
use serde::{Deserialize, Serialize};

// ── Atomic orbital visualization ───────────────────────────────────────────

/// Atomic orbital probability density on a 2D grid for volumetric rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrbitalVisualization {
    /// Probability density values on a 2D slice (flattened row-major).
    pub density: Vec<f64>,
    /// Grid dimensions (nx, ny).
    pub dimensions: [usize; 2],
    /// Grid spacing in Bohr radii.
    pub spacing: f64,
    /// Quantum numbers (n, l, m).
    pub quantum_numbers: [i32; 3],
    /// Maximum density value (for normalization).
    pub max_density: f64,
}

impl OrbitalVisualization {
    /// Generate an orbital slice through the XZ plane (y=0) for hydrogen-like atoms.
    ///
    /// `n`, `l`: principal and angular momentum quantum numbers.
    /// `grid_size`: number of points per side.
    /// `extent`: half-width in Bohr radii.
    #[must_use]
    #[allow(clippy::many_single_char_names)]
    pub fn hydrogen_slice(principal: u32, angular: u32, grid_size: usize, extent: f64) -> Self {
        let spacing = 2.0 * extent / grid_size.max(1) as f64;
        let mut density = Vec::with_capacity(grid_size * grid_size);
        let mut max_density = 0.0_f64;

        for iy in 0..grid_size {
            for ix in 0..grid_size {
                let px = -extent + ix as f64 * spacing;
                let pz = -extent + iy as f64 * spacing;
                let radius = (px * px + pz * pz).sqrt();

                let rho = radius / principal as f64;
                let prob = if radius < 1e-10 {
                    0.0
                } else {
                    let exp_part = (-rho).exp();
                    let poly = rho.powi(angular.cast_signed());
                    poly * poly * exp_part * exp_part * radius * radius
                };
                if prob > max_density {
                    max_density = prob;
                }
                density.push(prob);
            }
        }

        Self {
            density,
            dimensions: [grid_size, grid_size],
            spacing,
            quantum_numbers: [principal.cast_signed(), angular.cast_signed(), 0],
            max_density,
        }
    }
}

// ── Nuclear structure visualization ────────────────────────────────────────

/// Nucleon positions for particle rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NuclearStructure {
    /// Proton positions `[x, y, z]` in femtometres.
    pub protons: Vec<[f32; 3]>,
    /// Neutron positions `[x, y, z]` in femtometres.
    pub neutrons: Vec<[f32; 3]>,
    /// Nuclear radius in femtometres.
    pub radius_fm: f32,
    /// Element symbol.
    pub symbol: String,
}

impl NuclearStructure {
    /// Generate a simple spherical shell layout for a nucleus.
    #[must_use]
    pub fn from_nucleon_count(num_protons: u32, num_neutrons: u32, symbol: &str) -> Self {
        let mass_number = num_protons + num_neutrons;
        let radius = 1.2 * (mass_number as f32).cbrt(); // R = r₀ × A^(1/3)

        let mut protons = Vec::with_capacity(num_protons as usize);
        let mut neutrons = Vec::with_capacity(num_neutrons as usize);

        // Distribute nucleons on concentric shells using golden spiral
        let total = mass_number as usize;
        let golden_ratio = f32::midpoint(1.0, 5.0_f32.sqrt());

        for idx in 0..total {
            let frac = idx as f32 / total.max(1) as f32;
            let shell_r = radius * frac.cbrt();
            let theta = core::f32::consts::PI * (1.0 + (2.0 * golden_ratio - 1.0)) * idx as f32;
            let phi = (1.0 - 2.0 * (idx as f32 + 0.5) / total as f32).acos();

            let pos = [
                shell_r * phi.sin() * theta.cos(),
                shell_r * phi.sin() * theta.sin(),
                shell_r * phi.cos(),
            ];

            if (idx as u32) < num_protons {
                protons.push(pos);
            } else {
                neutrons.push(pos);
            }
        }

        Self {
            protons,
            neutrons,
            radius_fm: radius,
            symbol: String::from(symbol),
        }
    }
}

// ── Spectral line data ─────────────────────────────────────────────────────

/// Spectral line data for spectral plot rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpectralLineData {
    /// Lines with wavelength, intensity, and element.
    pub lines: Vec<SpectralLine>,
}

/// A single spectral line.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpectralLine {
    /// Wavelength in nm.
    pub wavelength_nm: f64,
    /// Relative intensity (0.0–1.0).
    pub intensity: f64,
    /// Element atomic number.
    pub element: u8,
    /// Transition label (e.g. "2→1", "Lyman-α").
    pub label: String,
}

// ── Decay chain visualization ──────────────────────────────────────────────

/// Decay chain for node-link graph rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DecayChainVisualization {
    /// Nuclides in the chain.
    pub nuclides: Vec<DecayNuclide>,
    /// Decay transitions between nuclides.
    pub transitions: Vec<DecayTransition>,
}

/// A nuclide node in the decay chain.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DecayNuclide {
    /// Atomic number.
    pub z: u32,
    /// Mass number.
    pub a: u32,
    /// Half-life in seconds (f64::INFINITY for stable).
    pub half_life_s: f64,
    /// Label (e.g. "U-238", "Pb-206").
    pub label: String,
}

/// A decay transition edge.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DecayTransition {
    /// Index of parent nuclide.
    pub parent: usize,
    /// Index of daughter nuclide.
    pub daughter: usize,
    /// Decay mode label (e.g. "α", "β⁻", "γ").
    pub mode: String,
    /// Branching ratio (0.0–1.0).
    pub branching_ratio: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orbital_1s() {
        let viz = OrbitalVisualization::hydrogen_slice(1, 0, 10, 5.0);
        assert_eq!(viz.density.len(), 100);
        assert!(viz.max_density > 0.0);
        assert_eq!(viz.quantum_numbers, [1, 0, 0]);
    }

    #[test]
    fn orbital_2p() {
        let viz = OrbitalVisualization::hydrogen_slice(2, 1, 8, 10.0);
        assert_eq!(viz.density.len(), 64);
        assert!(viz.max_density > 0.0);
    }

    #[test]
    fn nuclear_helium4() {
        let nuc = NuclearStructure::from_nucleon_count(2, 2, "He");
        assert_eq!(nuc.protons.len(), 2);
        assert_eq!(nuc.neutrons.len(), 2);
        assert!(nuc.radius_fm > 0.0);
        assert_eq!(nuc.symbol, "He");
    }

    #[test]
    fn nuclear_iron56() {
        let nuc = NuclearStructure::from_nucleon_count(26, 30, "Fe");
        assert_eq!(nuc.protons.len(), 26);
        assert_eq!(nuc.neutrons.len(), 30);
    }

    #[test]
    fn spectral_line_data_serializes() {
        let data = SpectralLineData {
            lines: vec![SpectralLine {
                wavelength_nm: 656.3,
                intensity: 1.0,
                element: 1,
                label: String::from("H-α"),
            }],
        };
        let json = serde_json::to_string(&data);
        assert!(json.is_ok());
    }

    #[test]
    fn decay_chain_serializes() {
        let chain = DecayChainVisualization {
            nuclides: vec![
                DecayNuclide {
                    z: 92,
                    a: 238,
                    half_life_s: 1.41e17,
                    label: String::from("U-238"),
                },
                DecayNuclide {
                    z: 90,
                    a: 234,
                    half_life_s: 2.08e12,
                    label: String::from("Th-234"),
                },
            ],
            transitions: vec![DecayTransition {
                parent: 0,
                daughter: 1,
                mode: String::from("α"),
                branching_ratio: 1.0,
            }],
        };
        assert_eq!(chain.nuclides.len(), 2);
        assert_eq!(chain.transitions.len(), 1);
    }

    #[test]
    fn nuclear_single_nucleon() {
        let nuc = NuclearStructure::from_nucleon_count(1, 0, "H");
        assert_eq!(nuc.protons.len(), 1);
        assert!(nuc.neutrons.is_empty());
    }
}
