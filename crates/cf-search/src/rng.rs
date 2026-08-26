//! Deterministic PRNG for the search engine.
//!
//! Ticket P1-23.
//!
//! Same construction as `cf-constraints`' `RepairRng`, for the same reason:
//! every stochastic decision in a run has to be reproducible from the seed
//! recorded in the lockfile, which rules out a thread-local RNG, and none of
//! this justifies a dependency. `GATE-DET-001` eventually requires bit-identical
//! output across thread counts, so the streams must be assigned per block or
//! per individual rather than shared — that assignment is `P1-33`'s, and this
//! type is what it will hand out.

/// xorshift64*, seeded per stream.
#[derive(Clone, Debug)]
pub struct SearchRng(u64);

impl SearchRng {
    /// Seed a stream.
    ///
    /// A zero state is absorbing for xorshift, so it has to be excluded — but
    /// the obvious exclusion, `seed | 1`, silently maps every even seed onto
    /// its odd successor. Seeds `2026` and `2027` would then produce
    /// *identical runs*, and since the root seed comes from the lockfile, two
    /// jobs a user deliberately seeded differently would agree and look like a
    /// reproducibility success.
    ///
    /// The splitmix64 finalizer is a bijection on `u64`, so distinct seeds stay
    /// distinct; the initial XOR keeps zero from mapping to zero, and the
    /// remaining single preimage of zero is redirected.
    pub fn new(seed: u64) -> Self {
        let mut z = seed ^ 0x9E37_79B9_7F4A_7C15;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        Self(if z == 0 { 0x2545_F491_4F6C_DD1D } else { z })
    }

    /// A child stream for one block or individual.
    ///
    /// Derived rather than drawn, so the same `(seed, index)` gives the same
    /// stream regardless of how many other streams were created first — which
    /// is what makes thread-count independence achievable at all.
    pub fn stream(seed: u64, index: u64) -> Self {
        Self::new(seed ^ index.wrapping_mul(0x9E37_79B9_7F4A_7C15))
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    /// Uniform in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in `0..n`. Returns 0 for `n == 0`.
    pub fn below(&mut self, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        // Lemire's multiply-shift, without the rejection loop: the bias is
        // below 2^-64 for any n the search will ever use, and a rejection loop
        // would make the number of RNG draws data-dependent, which breaks
        // stream reproducibility across thread counts.
        ((self.next_u64() as u128 * n as u128) >> 64) as usize
    }

    /// A draw from `Geometric(p)` on `{1, 2, ...}`, by inversion.
    ///
    /// Clamped to `max` so a very small `p` cannot produce an unbounded block
    /// length; the clamp is a bound on work, not a distributional claim.
    pub fn geometric(&mut self, p: f64, max: usize) -> usize {
        if p >= 1.0 {
            return 1;
        }
        if p <= 0.0 {
            return max;
        }
        let u = self.next_f64().max(f64::MIN_POSITIVE);
        let draw = (u.ln() / (1.0 - p).ln()).floor() as usize + 1;
        draw.clamp(1, max.max(1))
    }

    /// Index sampled proportionally to `weights`.
    ///
    /// Returns `None` when every weight is zero, rather than falling back to a
    /// uniform draw: an all-zero weight vector means the caller has no
    /// admissible choice, and quietly inventing one would hide that.
    pub fn weighted(&mut self, weights: &[f64]) -> Option<usize> {
        let total: f64 = weights.iter().filter(|w| w.is_finite() && **w > 0.0).sum();
        // `total` cannot be NaN here — every summand was checked finite and
        // positive — so a plain comparison is sound as well as readable.
        if total <= 0.0 {
            return None;
        }
        let mut target = self.next_f64() * total;
        for (index, weight) in weights.iter().enumerate() {
            if !weight.is_finite() || *weight <= 0.0 {
                continue;
            }
            target -= weight;
            if target <= 0.0 {
                return Some(index);
            }
        }
        // Reachable only through floating-point drift at the very top of the
        // range; the last positive weight is the right answer there.
        weights.iter().rposition(|w| w.is_finite() && *w > 0.0)
    }

    pub fn coin(&mut self) -> bool {
        self.next_u64() & 1 == 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_streams_are_reproducible_and_distinct() {
        let a: Vec<u64> = (0..8).map(|_| SearchRng::stream(7, 3).next_u64()).collect();
        assert!(a.windows(2).all(|w| w[0] == w[1]), "same stream, same draw");

        let mut x = SearchRng::stream(7, 3);
        let mut y = SearchRng::stream(7, 4);
        assert_ne!(x.next_u64(), y.next_u64());

        // Order-independence: stream 4 is the same whether or not stream 3 was
        // created first. Without this, thread count changes the numbers.
        let mut z = SearchRng::stream(7, 4);
        assert_eq!(SearchRng::stream(7, 4).next_u64(), z.next_u64());
    }

    /// Adjacent seeds must not alias.
    ///
    /// `seed | 1` — the obvious way to keep xorshift out of its absorbing zero
    /// state — maps every even seed onto its odd successor, so half of all
    /// consecutive seed pairs produce identical runs. Caught by a determinism
    /// test that expected two seeds to diverge and found they did not.
    #[test]
    fn unit_distinct_seeds_give_distinct_streams() {
        let first: Vec<u64> = (0..64u64)
            .map(|seed| SearchRng::new(seed).next_u64())
            .collect();
        let unique: std::collections::BTreeSet<u64> = first.iter().copied().collect();
        assert_eq!(
            unique.len(),
            first.len(),
            "seeds 0..64 produced only {} distinct streams",
            unique.len()
        );

        // The pair that found this.
        assert_ne!(
            SearchRng::new(2026).next_u64(),
            SearchRng::new(2027).next_u64()
        );
        // And zero is a usable seed rather than an absorbing state.
        let mut zero = SearchRng::new(0);
        assert_ne!(zero.next_u64(), 0);
        assert_ne!(zero.next_u64(), 0);
    }

    #[test]
    fn unit_below_stays_in_range() {
        let mut rng = SearchRng::new(1);
        for n in [1usize, 2, 3, 7, 64, 1000] {
            for _ in 0..2000 {
                assert!(rng.below(n) < n);
            }
        }
        assert_eq!(rng.below(0), 0);
    }

    #[test]
    fn unit_weighted_reports_an_empty_choice_rather_than_inventing_one() {
        let mut rng = SearchRng::new(9);
        assert_eq!(rng.weighted(&[0.0, 0.0, 0.0]), None);
        assert_eq!(rng.weighted(&[]), None);
        assert_eq!(rng.weighted(&[0.0, 1.0, 0.0]), Some(1));

        let mut counts = [0usize; 3];
        for _ in 0..30_000 {
            counts[rng.weighted(&[1.0, 3.0, 0.0]).unwrap()] += 1;
        }
        assert_eq!(counts[2], 0, "a zero weight is never chosen");
        let ratio = counts[1] as f64 / counts[0] as f64;
        assert!((2.7..3.3).contains(&ratio), "ratio {ratio} is not near 3");
    }

    #[test]
    fn unit_geometric_is_bounded() {
        let mut rng = SearchRng::new(5);
        for _ in 0..10_000 {
            let draw = rng.geometric(1.0 / 30.0, 200);
            assert!((1..=200).contains(&draw));
        }
        assert_eq!(rng.geometric(1.0, 200), 1);
    }
}
