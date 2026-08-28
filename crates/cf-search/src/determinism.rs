//! The determinism harness.
//!
//! Ticket P1-33. Gate: `GATE-DET-001` — bit-identical output across 1, 4 and 16
//! threads, and across two machines of the same architecture.
//!
//! ## What this measures, and what it could not measure before
//!
//! `GATE-DET-001` is about *parallel* execution. Until `P1-33` there was no
//! parallelism anywhere in the workspace — no rayon, no threads — so a harness
//! comparing "1, 4 and 16 threads" would have run the same sequential code
//! three times and passed unconditionally. It would have been the most
//! reassuring test in the repository and worth nothing.
//!
//! So the parallelism came first: tier 2's same-colour blocks are proposed
//! through `rayon`, which is what ADR-0001 named for this sweep and what §3.1's
//! conditional independence exists to make sound. This module is what checks
//! it, and its own negative control is a scheduling-dependent implementation:
//! [`compare_thread_counts`] must fail for one.
//!
//! ## The two claims, kept separate
//!
//! * **Thread-count independence.** Same input, same seed, different pool
//!   size: identical bytes. Rests on streams keyed by block index rather than
//!   by iteration order, and on a commit ordered by block index rather than by
//!   completion.
//! * **Split independence.** A run stopped and resumed equals an uninterrupted
//!   one. `P1-29` already established this for the segment runner, and
//!   [`compare_splits`] is that fixture generalised rather than a second copy
//!   of it — a harness that re-derived the property could pass while the
//!   thing it was meant to guard had regressed.

use std::fmt::Debug;

use cf_types::{code, CfError, Result};

/// Thread counts `GATE-DET-001` names.
pub const GATE_THREAD_COUNTS: [usize; 3] = [1, 4, 16];

/// Run `run` under a rayon pool of exactly `threads` threads.
///
/// A private pool rather than the global one, so a harness comparing counts
/// does not depend on process-wide state and several comparisons can run in
/// the same test binary without fighting over it.
pub fn with_threads<T: Send>(threads: usize, run: impl Fn() -> T + Sync + Send) -> Result<T> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| {
            CfError::new(
                code::SRCH_DETERMINISM,
                format!("could not build a {threads}-thread pool: {e}"),
            )
        })?;
    Ok(pool.install(run))
}

/// Compare a run's output across thread counts.
///
/// `observe` must return something that is `Eq`: comparing a summary statistic
/// rather than the output itself is how a determinism harness passes while the
/// output differs, so callers are expected to hand over the whole result.
pub fn compare_thread_counts<T: Eq + Debug + Send>(
    counts: &[usize],
    run: impl Fn() -> T + Sync + Send,
) -> Result<T> {
    let (first, rest) = counts
        .split_first()
        .ok_or_else(|| CfError::new(code::SRCH_DETERMINISM, "no thread counts to compare"))?;

    let reference = with_threads(*first, &run)?;
    for &threads in rest {
        let observed = with_threads(threads, &run)?;
        if observed != reference {
            return Err(CfError::new(
                code::SRCH_DETERMINISM,
                format!(
                    "output differs between {first} and {threads} threads; a parallel \
                     search whose answer depends on the pool size is not reproducible, \
                     and GATE-DET-001 blocks release on it"
                ),
            ));
        }
    }
    Ok(reference)
}

/// Compare a run advanced in one piece against the same run advanced in pieces.
///
/// `advance(state, n)` must move the run forward `n` steps. Every split of the
/// same total is checked against the whole, because a run that resumed
/// correctly at one boundary and not another would pass a single-split test —
/// which is the shape `P1-29`'s fixture already had, and the reason this
/// generalises it rather than replacing it.
pub fn compare_splits<S, O>(
    total: u32,
    splits: &[Vec<u32>],
    start: impl Fn() -> S,
    advance: impl Fn(S, u32) -> Result<S>,
    observe: impl Fn(&S) -> O,
) -> Result<()>
where
    O: Eq + Debug,
{
    for split in splits {
        let sum: u32 = split.iter().sum();
        if sum != total {
            return Err(CfError::new(
                code::SRCH_DETERMINISM,
                format!("split {split:?} sums to {sum}, not {total}; it is not the same run"),
            ));
        }
    }

    let mut whole = start();
    whole = advance(whole, total)?;
    let reference = observe(&whole);

    for split in splits {
        let mut state = start();
        for &piece in split {
            state = advance(state, piece)?;
        }
        if observe(&state) != reference {
            return Err(CfError::new(
                code::SRCH_DETERMINISM,
                format!(
                    "a run advanced in pieces {split:?} differs from one advanced in a \
                     single step of {total}; resuming is not transparent"
                ),
            ));
        }
    }
    Ok(())
}

/// The splits `P1-29` established, reused rather than re-derived.
pub fn standard_splits() -> Vec<Vec<u32>> {
    vec![
        vec![6, 6],
        vec![4, 4, 4],
        vec![3, 3, 3, 3],
        vec![1, 11],
        vec![11, 1],
        vec![1; 12],
        vec![5, 2, 5],
    ]
}

/// Sum in a fixed order, independent of how the values were produced.
///
/// Floating-point addition is not associative, so a reduction whose order
/// follows completion order gives a different answer at a different thread
/// count — the same class of failure as an unordered stream assignment, but
/// arriving through arithmetic rather than through control flow. Reductions
/// that cross a parallel boundary go through here.
pub fn ordered_sum(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |acc, value| acc + value)
}
