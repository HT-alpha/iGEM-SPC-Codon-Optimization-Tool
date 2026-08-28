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

/// xoshiro256++, seeded per stream.
///
/// The algorithm is named, not incidental. `landscape_lockfile.schema.json`'s
/// `rng.algorithm` is a closed enum — `pcg64`, `chacha8`, `xoshiro256pp` — and a
/// job has to record which one produced its numbers. This was plain xorshift64
/// (shifts 13/7/17), which is none of the three, so a lockfile would have had
/// to name an algorithm that did not run and a replay trusting it would have
/// been relying on a false record. `DEP-DET-001`: the implementation was
/// changed to match the contract rather than the enum widened to admit a weaker
/// generator, because `GATE-DET-001` and `P1-33` both rest on this layer.
#[derive(Clone, Debug)]
pub struct SearchRng([u64; 4]);

#[inline]
fn rotl(x: u64, k: u32) -> u64 {
    x.rotate_left(k)
}

impl SearchRng {
    /// Seed a stream.
    ///
    /// splitmix64 expands the seed into the 256-bit state, which is what
    /// xoshiro's own reference recommends and what keeps a poor seed — zero,
    /// or a small integer — from producing a poor state. It is also a
    /// bijection, so distinct seeds stay distinct.
    ///
    /// That last property is not decorative. An earlier version used
    /// `seed | 1` to keep xorshift out of its absorbing zero state, which maps
    /// every even seed onto its odd successor: seeds `2026` and `2027` produced
    /// identical runs. Since the root seed comes from the lockfile, two jobs a
    /// user seeded deliberately differently would have agreed and looked like a
    /// reproducibility success.
    pub fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut state = [0u64; 4];
        for slot in state.iter_mut() {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            *slot = x ^ (x >> 31);
        }
        // An all-zero state is the one xoshiro cannot leave. splitmix64 reaches
        // it from no input, but the state is not the caller's to choose and a
        // guard costs one comparison at construction.
        if state == [0, 0, 0, 0] {
            state = [
                0x2545_F491_4F6C_DD1D,
                0x9E37_79B9_7F4A_7C15,
                0xBF58_476D_1CE4_E5B9,
                0x94D0_49BB_1331_11EB,
            ];
        }
        Self(state)
    }

    /// The algorithm name this generator records in a lockfile.
    ///
    /// A constant rather than a string literal at the call site, so the record
    /// and the implementation cannot drift apart silently.
    pub const ALGORITHM: &'static str = "xoshiro256pp";

    /// A child stream for one block, generation or individual.
    ///
    /// Derived rather than drawn, so the same `(seed, index)` gives the same
    /// stream regardless of how many other streams were created first — which
    /// is what makes thread-count independence achievable at all.
    pub fn stream(seed: u64, index: u64) -> Self {
        Self::new(seed ^ index.wrapping_mul(0x9E37_79B9_7F4A_7C15))
    }

    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.0;
        let result = rotl(s[0].wrapping_add(s[3]), 23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        result
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
