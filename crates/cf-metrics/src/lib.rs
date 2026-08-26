//! `cf-metrics` — support layer.
//!
//! Pure multi-objective indicator math shared by the search engine and the
//! decision-support layer: non-dominated sorting, crowding distance,
//! Das-Dennis reference directions, WFG hypervolume, IGD+, epsilon indicator.
//! Deliberately a leaf so that `cf-search` and `cf-decision` never need to
//! depend on each other.
//!
//! Declared dependencies: none.
//! Contract: `planning/03-algorithms/03-nsga-loop-and-operators.md`.

#![forbid(unsafe_code)]
