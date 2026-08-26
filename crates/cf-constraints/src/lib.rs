//! `cf-constraints` — constraint-repair layer.
//!
//! Hard-constraint specification, incremental violation detection, the
//! minimal-change repair operator, the constraint-feasibility analyser, and
//! the tabu ledger that keeps repair from fighting the search loop. Pure
//! function of (sequence, ConstraintSet, optional delta table) - it has no
//! knowledge of populations, generations or objectives.
//!
//! Declared dependencies: cf-types.
//! Contract: `planning/02-schemas/constraint_set.schema.json`.

#![forbid(unsafe_code)]
