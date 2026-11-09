//! `src/optics_data.cyr` — the Storey & Hummer Case B Balmer intensities (port plan M6).
//!
//! The intensities come through the oracle's public `hydrogen_balmer_series`, each the f64
//! the oracle computed (Hβ/Hα = 1/2.8632, …), so none is typed here. The upper levels n = 3..=8
//! are proven: line k's wavelength must be `spectral_line_vacuum_nm(1, 1, 2, n_k)`, bit for bit.

use crate::cyr::{Elem, header, table};
use std::fmt::Write as _;
use tanmatra::atomic::spectral_line_vacuum_nm;
use tanmatra::optics::hydrogen_balmer_series;

const N_UPPER: [u8; 6] = [3, 4, 5, 6, 7, 8];
const NAMES: [&str; 6] = [
    "Hα (n=3)",
    "Hβ (n=4)",
    "Hγ (n=5)",
    "Hδ (n=6)",
    "Hε (n=7)",
    "H8 (n=8)",
];

pub fn emit(provenance: &str) -> Result<String, String> {
    let lines = hydrogen_balmer_series();
    if lines.len() != N_UPPER.len() {
        return Err(format!(
            "hydrogen_balmer_series: {} lines, expected 6",
            lines.len()
        ));
    }
    for (l, &n) in lines.iter().zip(&N_UPPER) {
        let want = spectral_line_vacuum_nm(1, 1, 2, u32::from(n)).map_err(|e| format!("{e:?}"))?;
        if want.to_bits() != l.wavelength_nm.to_bits() {
            return Err(format!(
                "Balmer line n={n}: the wavelength does not match spectral_line_vacuum_nm"
            ));
        }
    }
    let mut out = header(
        "optics data: the Balmer Case B intensities",
        provenance,
        &[
            "Through the oracle's hydrogen_balmer_series; the upper levels are proven against its wavelengths\n\
           (scripts/oracle/src/t_optics.rs). tests/golden/optics.tsv pins every use.",
        ],
    );
    let lab: Vec<String> = NAMES.iter().map(|s| (*s).to_string()).collect();
    let _ = writeln!(
        out,
        "# The six Balmer lines, Hα through H8.\nconst _TM_BALMER_COUNT = {};\n",
        N_UPPER.len()
    );
    out.push_str(&table(
        "_TM_BALMER_N_UPPER",
        "Column: the upper level n (the lower is 2).",
        &Elem::U8(N_UPPER.to_vec()),
        &lab,
    )?);
    out.push_str(&table(
        "_TM_BALMER_INTENSITY",
        "Column: Case B emissivity relative to Hα, T_e = 10⁴ K, n_e = 10² cm⁻³ (Storey & Hummer, MNRAS\n\
         272, 41 (1995), VizieR VI/64): Hα/Hβ = 2.8632, Hγ/Hβ = 0.4683, Hδ/Hβ = 0.2589, Hε/Hβ = 0.1590,\n\
         H8/Hβ = 0.1050.",
        &Elem::F64(lines.iter().map(|l| l.intensity).collect()),
        &lab,
    )?);
    Ok(out.trim_end().to_string() + "\n")
}
