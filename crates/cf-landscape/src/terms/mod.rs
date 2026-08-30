//! The built-in terms.
//!
//! Ticket P1-21. Three, and only three:
//!
//! | term | capability | role cap |
//! | --- | --- | --- |
//! | [`CaiTerm`] (`usage.cai`) | L1 | Pareto axis |
//! | [`GcContentTerm`] (`usage.gc_content`) | L1 | diagnostic |
//! | [`ViolationCountTerm`] (`constraint.violation_count`) | L1 or L0, per set | Pareto axis |
//!
//! Every other term named in the planning package arrives later and with its
//! own ticket: the epistatic surface in `P3-05` (a synonymous-restricted Potts
//! model with `O(W·S)` exact deltas, specified in detail), windowed structure
//! accessibility in `P2-11`, the burden ledger in `P2-14`, folding kinetics in
//! `P2-17`, LM naturalness in `P3-14`.
//!
//! None of them is stubbed here. A stub that returned a plausible number would
//! be indistinguishable from the real term at every call site, and a search
//! optimizing against a placeholder epistatic surface would produce results
//! that look right. Asking for an unregistered term fails with `CF-LAND-001`
//! naming it, which is the loud failure the standing rules require.

pub mod constraints;
pub mod structure;
pub mod usage;

pub use constraints::ViolationCountTerm;
pub use usage::{CaiTerm, CodonUsageTable, GcContentTerm, UsageSubset};
