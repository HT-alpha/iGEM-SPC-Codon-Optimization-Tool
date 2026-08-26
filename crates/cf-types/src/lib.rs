//! `cf-types` — foundation layer.
//!
//! Canonical codon-indexed sequence representation, genetic-code tables,
//! identifier newtypes, units, the objective vector, and the `CF-*` error
//! registry. Contains no biology model and no I/O.
//!
//! Declared dependencies: none. This is a rank-0 leaf in the module graph
//! (`planning/01-modules/modules.yaml`), and it must stay one: everything else
//! depends on it, so a dependency added here becomes a dependency of the whole
//! system.
//!
//! Contract: `planning/02-schemas/sequence.schema.json` and
//! `planning/02-schemas/common.schema.json`.
//!
//! Tickets: P1-04 (errors), P1-05 (codons and tables), P1-06 (sequence),
//! P1-07 (canonical JSON and digests), P1-08 (identifiers and units).

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod canonical;
pub mod codon;
pub mod error;
pub mod ids;
pub mod seqrec;
pub mod tables;
pub mod units;

pub use canonical::{canonicalize, Digest};
pub use codon::{Codon, BASES, CODON_COUNT};
pub use error::{code, CfCode, CfError, Result, REGISTRY};
pub use ids::{ArtifactId, HostId, JobId, SchemaVersion, SemVer, TenantId, TermId};
pub use seqrec::{CodonSequence, LockReason, LockedSpan};
pub use tables::GeneticCode;
pub use units::{Direction, ObjectiveVector, Unit};

/// Version of this crate, reported through the compiled core so that a running
/// system can be tied back to a build.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
