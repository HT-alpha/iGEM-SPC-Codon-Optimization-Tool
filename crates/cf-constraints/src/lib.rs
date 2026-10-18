//! `cf-constraints` — constraint-repair layer.
//!
//! Hard-constraint specification, incremental violation detection, the
//! minimal-change repair operator, the constraint-feasibility analyser, and the
//! tabu ledger that keeps repair from fighting the search loop.
//!
//! This crate is a pure function of `(sequence, ConstraintSet, optional delta
//! table)`. It has no knowledge of populations, generations or objectives, and
//! it must not acquire any: that is what keeps the module graph acyclic
//! (`planning/01-modules/02-module-dependency-graph.md` section 3.2) and what
//! makes repair testable against a hand-written delta table with no search
//! engine present.
//!
//! Declared dependencies: cf-types.
//! Contract: `planning/02-schemas/constraint_set.schema.json`.
//!
//! Tickets landed: P1-13 (spec and expansion), P1-14 (motif automaton),
//! P1-15 (k-mer index), P1-16 (homopolymer and GC), P1-19 (feasibility
//! analyser), P1-17 (repair operator), P1-18a (tabu ledger).
//!
//! P1-19 precedes P1-17: repair's termination is only principled in
//! combination with an analyser that rejects unsatisfiable constraint sets
//! once, at submission, rather than letting repair rediscover them on every
//! candidate. See planning/03-algorithms/06-constraint-repair.md section 5.1.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod detect;
pub mod feasibility;
pub mod kmer;
pub mod legacy;
pub mod motif;
pub mod repair;
pub mod seed;
pub mod spec;
pub mod spend;
pub mod tabu;
pub mod window;

pub use detect::{ConstraintChecker, ConstructView, KmerMasking, Violation};
pub use feasibility::{
    analyse, CheckId, DesignSpace, FeasibilityReport, FeasibilityWitness, Finding, Satisfiability,
    WitnessScope,
};
pub use kmer::{count_overlapping, KmerIndex};
pub use motif::MotifAutomaton;
pub use repair::{
    partition, repair, ActionableViolation, Edit, GiveUp, GiveUpReason, NoObjective,
    ObjectiveDelta, RepairOutcome, RepairRng, RepairTarget, WorkAccounting,
};
pub use seed::{
    explain_deviations, resolve_feasible_seed, resolve_feasible_seed_limited, BaselineDeviation,
    SeedReport, SeedResidual, UnsearchedRule,
};
pub use spec::{
    AssemblySpec, Constraint, ConstraintClass, ConstraintSet, ExpansionReport, RepairPolicy, Rule,
    Scope, Strand,
};
pub use spend::SpendLimit;
pub use tabu::TabuLedger;
pub use window::{find_runs, GcIndex, Run};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
