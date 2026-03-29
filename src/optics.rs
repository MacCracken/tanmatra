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
/// Returns a prakash `Spd` sampled over the specified range.
#[must_use]
pub fn lines_to_spd(
    lines: &[SpectralLine],
    fwhm_nm: f64,
    start_nm: f64,
    end_nm: f64,
    step_nm: f64,
) -> Spd {
    let sigma = fwhm_nm / (2.0 * libm::sqrt(2.0 * libm::log(2.0)));
    let inv_2sigma2 = 1.0 / (2.0 * sigma * sigma);

    #[allow(clippy::cast_sign_loss)]
    let num_samples = ((end_nm - start_nm) / step_nm).ceil().max(0.0) as usize + 1;
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
/// Returns spectral lines for transitions n=3..=8 -> n=2 (H-alpha through H-zeta).
#[must_use]
pub fn hydrogen_balmer_series() -> Vec<SpectralLine> {
    let mut lines = Vec::new();
    // Balmer series: n_upper -> n_lower=2
    let names_and_intensities = [
        (3, 1.0),  // H-alpha: ~656.3 nm
        (4, 0.35), // H-beta: ~486.1 nm
        (5, 0.15), // H-gamma: ~434.0 nm
        (6, 0.08), // H-delta: ~410.2 nm
        (7, 0.04), // H-epsilon: ~397.0 nm
        (8, 0.02), // H-zeta: ~388.9 nm
    ];

    for &(n_upper, intensity) in &names_and_intensities {
        if let Ok(wavelength) = crate::atomic::spectral_line_nm(1, 2, n_upper) {
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
/// Returns spectral lines with relative intensities decreasing with n.
#[must_use]
pub fn spectral_series_lines(z: u32, n_lower: u32, n_upper_max: u32) -> Vec<SpectralLine> {
    let mut lines = Vec::new();

    for n_upper in (n_lower + 1)..=n_upper_max {
        if let Ok(wavelength) = crate::atomic::spectral_line_nm(z, n_lower, n_upper) {
            // Intensity roughly scales as 1/n³ (from oscillator strengths)
            let intensity = 1.0 / ((n_upper as f64) * (n_upper as f64) * (n_upper as f64));
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
