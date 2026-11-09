//! Shared input grids.
//!
//! Every value here must be reproducible bit for bit on any machine, so `--check` passes
//! anywhere. Use the pinned `libm` crate (pure Rust) for any transcendental — never the
//! `std` f64 methods, which call the platform's C library.

/// Special values every float function is probed with.
pub fn specials() -> Vec<f64> {
    vec![
        0.0,
        -0.0,
        1.0,
        -1.0,
        f64::MIN_POSITIVE,
        f64::from_bits(1), // smallest subnormal
        f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
    ]
}

/// `lo · 10^(k / per_decade)` for k = 0 ..= decades · per_decade. `lo` must be positive.
pub fn log_grid(lo: f64, decades: u32, per_decade: u32) -> Vec<f64> {
    (0..=decades * per_decade)
        .map(|k| lo * libm::pow(10.0, f64::from(k) / f64::from(per_decade)))
        .collect()
}

/// Deterministic xorshift64* — fixed seeds keep the sweeps identical across runs.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    /// Uniform in [0, 1): the top 53 bits as a fraction.
    pub fn unit(&mut self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let top = (self.next_u64() >> 11) as f64;
        top * (1.0 / 9_007_199_254_740_992.0) // 2^-53
    }

    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// `10^u` with u uniform in [lo_exp10, hi_exp10).
    pub fn log_uniform(&mut self, lo_exp10: f64, hi_exp10: f64) -> f64 {
        libm::pow(10.0, self.uniform(lo_exp10, hi_exp10))
    }

    /// As `log_uniform`, with a random sign.
    pub fn signed_log_uniform(&mut self, lo_exp10: f64, hi_exp10: f64) -> f64 {
        let x = self.log_uniform(lo_exp10, hi_exp10);
        if self.next_u64() & 1 == 1 { -x } else { x }
    }

    /// Uniform integer in [lo, hi].
    pub fn int(&mut self, lo: i64, hi: i64) -> i64 {
        #[allow(clippy::cast_possible_wrap, clippy::cast_sign_loss)]
        let span = (hi - lo + 1) as u64;
        #[allow(clippy::cast_possible_wrap)]
        let off = (self.next_u64() % span) as i64;
        lo + off
    }
}
