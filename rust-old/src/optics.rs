//! Integration with the prakash optics crate.
//!
//! Provides conversion from tanmatra spectral line data to prakash
//! spectral power distributions (`Spd`), enabling visualization and
//! optical simulation of atomic emission and absorption spectra.
//!
//! Requires the `optics` feature flag.

use alloc::vec::Vec;
use prakash::spectral::Spd;

/// A spectral emission or absorption line.
///
/// Represents a discrete spectral feature at a specific wavelength with
/// an associated relative intensity. Used as an intermediate format for
/// building spectral power distributions.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SpectralLine {
    /// Wavelength of the line in nanometers.
    pub wavelength_nm: f64,
    /// Relative intensity (arbitrary units, >= 0).
    pub intensity: f64,
}

/// Builds a prakash `Spd` (spectral power distribution) from a set of
/// spectral lines using Gaussian line profiles.
///
/// Each line is represented as a Gaussian peak centered at its wavelength
/// with the specified full-width at half-maximum (FWHM).
///
/// Parameters:
/// - `lines`: spectral lines with wavelengths and intensities
/// - `fwhm_nm`: line width (FWHM) in nm (e.g., 0.1 for narrow lines)
/// - `start_nm`: start of the SPD wavelength range
/// - `end_nm`: end of the SPD wavelength range
/// - `step_nm`: wavelength step size for sampling
///
/// Returns a prakash `Spd` sampled over the specified range. Invalid ranges
/// (non-finite bounds, `end_nm < start_nm`, `step_nm <= 0`, or more than
/// 10⁷ samples) return an empty distribution; a non-positive FWHM returns
/// zero power at every sample.
#[must_use]
pub fn lines_to_spd(
    lines: &[SpectralLine],
    fwhm_nm: f64,
    start_nm: f64,
    end_nm: f64,
    step_nm: f64,
) -> Spd {
    if !(start_nm.is_finite() && end_nm.is_finite() && step_nm > 0.0 && step_nm.is_finite())
        || end_nm < start_nm
    {
        return Spd::new(
            if start_nm.is_finite() { start_nm } else { 0.0 },
            1.0,
            Vec::new(),
        );
    }
    let span_steps = ((end_nm - start_nm) / step_nm).ceil();
    if span_steps.is_nan() || span_steps > 1.0e7 {
        return Spd::new(start_nm, step_nm, Vec::new());
    }
    #[allow(clippy::cast_sign_loss)]
    let num_samples = span_steps as usize + 1;
    if fwhm_nm.is_nan() || fwhm_nm <= 0.0 {
        return Spd::new(start_nm, step_nm, alloc::vec![0.0; num_samples]);
    }
    let sigma = fwhm_nm / (2.0 * libm::sqrt(2.0 * libm::log(2.0)));
    let inv_2sigma2 = 1.0 / (2.0 * sigma * sigma);
    let mut values = Vec::with_capacity(num_samples);

    for i in 0..num_samples {
        let wavelength = start_nm + i as f64 * step_nm;
        let mut power = 0.0;

        for line in lines {
            let dw = wavelength - line.wavelength_nm;
            power += line.intensity * libm::exp(-dw * dw * inv_2sigma2);
        }

        values.push(power);
    }

    Spd::new(start_nm, step_nm, values)
}

/// Generates the hydrogen Balmer series emission lines (visible spectrum).
///
/// Returns spectral lines for transitions n=3..=8 -> n=2 (Hα through H8) at
/// their vacuum wavelengths including the proton reduced mass (Hα 656.470 nm).
/// Intensities are Case B recombination emissivities relative to Hα for
/// T_e = 10⁴ K and n_e = 10² cm⁻³ (Storey & Hummer, MNRAS 272, 41 (1995),
/// VizieR VI/64): Hα/Hβ = 2.8632, Hγ/Hβ = 0.4683, Hδ/Hβ = 0.2589,
/// Hε/Hβ = 0.1590, H8/Hβ = 0.1050.
#[must_use]
pub fn hydrogen_balmer_series() -> Vec<SpectralLine> {
    let mut lines = Vec::new();
    // (n_upper, emissivity relative to Hα)
    let names_and_intensities = [
        (3, 1.0),
        (4, 1.0 / 2.8632),
        (5, 0.4683 / 2.8632),
        (6, 0.2589 / 2.8632),
        (7, 0.1590 / 2.8632),
        (8, 0.1050 / 2.8632),
    ];

    for &(n_upper, intensity) in &names_and_intensities {
        if let Ok(wavelength) = crate::atomic::spectral_line_vacuum_nm(1, 1, 2, n_upper) {
            lines.push(SpectralLine {
                wavelength_nm: wavelength,
                intensity,
            });
        }
    }

    lines
}

/// Generates emission lines for a hydrogen-like atom for a given spectral series.
///
/// Parameters:
/// - `z`: atomic number
/// - `n_lower`: lower level of the series (1=Lyman, 2=Balmer, 3=Paschen, etc.)
/// - `n_upper_max`: maximum upper level to include
///
/// Wavelengths assume an infinitely heavy nucleus. The relative intensity of
/// each line is its total spontaneous emission rate Σ A(n_upper l′ → n_lower l)
/// summed over all dipole-allowed l, l′ pairs and weighted by (2l′ + 1), i.e.
/// a statistically populated upper level; values are normalized to the first
/// line. (For recombination-dominated nebulae use tabulated Case B intensities.)
#[must_use]
pub fn spectral_series_lines(z: u32, n_lower: u32, n_upper_max: u32) -> Vec<SpectralLine> {
    let mut lines = Vec::new();
    let mut first: Option<f64> = None;

    for n_upper in (n_lower + 1)..=n_upper_max {
        if let Ok(wavelength) = crate::atomic::spectral_line_nm(z, n_lower, n_upper) {
            let mut emission = 0.0;
            for lu in 0..n_upper {
                for ll in 0..n_lower {
                    if let Ok(a) =
                        crate::atomic::einstein_a_coefficient(z, n_upper, lu, n_lower, ll)
                    {
                        emission += f64::from(2 * lu + 1) * a;
                    }
                }
            }
            // Photon energy weighting is not applied: intensities are photon rates.
            let norm = *first.get_or_insert(emission);
            let intensity = if norm > 0.0 { emission / norm } else { 0.0 };
            lines.push(SpectralLine {
                wavelength_nm: wavelength,
                intensity,
            });
        }
    }

    lines
}

/// Converts a spectral line wavelength to a prakash RGB color.
///
/// Uses prakash's `wavelength_to_rgb` for accurate visible-spectrum colors.
///
/// # Errors
///
/// Returns [`prakash::PrakashError`] if the wavelength is out of visible range.
pub fn line_to_rgb(line: &SpectralLine) -> Result<prakash::spectral::Rgb, prakash::PrakashError> {
    prakash::spectral::wavelength_to_rgb(line.wavelength_nm)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn balmer_series_has_6_lines() {
        let lines = hydrogen_balmer_series();
        assert_eq!(lines.len(), 6);
        // H-alpha should be ~656 nm
        assert!((lines[0].wavelength_nm - 656.3).abs() < 1.0);
    }

    #[test]
    fn lines_to_spd_produces_valid_spd() {
        let lines = hydrogen_balmer_series();
        let spd = lines_to_spd(&lines, 2.0, 380.0, 780.0, 1.0);
        // SPD should have ~401 samples for 380-780nm at 1nm step
        let val = spd.at(656.0);
        assert!(val > 0.0, "SPD should have power at H-alpha");
    }

    #[test]
    fn spectral_series_lyman() {
        let lines = spectral_series_lines(1, 1, 6);
        assert_eq!(lines.len(), 5); // n=2..=6 -> 5 lines
        // Lyman-alpha should be ~121.6 nm
        assert!((lines[0].wavelength_nm - 121.6).abs() < 1.0);
    }

    #[test]
    fn h_alpha_visible_color() {
        let line = SpectralLine {
            wavelength_nm: 656.3,
            intensity: 1.0,
        };
        let rgb = line_to_rgb(&line).unwrap();
        // H-alpha is red
        assert!(rgb.r > rgb.b, "H-alpha should be reddish");
    }

    #[test]
    fn serde_roundtrip_spectral_line() {
        let line = SpectralLine {
            wavelength_nm: 486.1,
            intensity: 0.35,
        };
        let json = serde_json::to_string(&line).unwrap();
        let back: SpectralLine = serde_json::from_str(&json).unwrap();
        assert_eq!(line, back);
    }
}
