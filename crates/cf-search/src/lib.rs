//! `cf-search` — three-tier search engine.
//!
//! Ticket P1-23 delivers the genotype and the variation operators; tiers 1–3
//! and the joint-search driver follow in `P1-24` … `P1-26`.
//!
//! Talks only to `cf-landscape-api`. It never depends on `cf-landscape`, on a
//! predictor crate, or on `cf-epistasis` — the epistatic conditional two of the
//! mutation operators need arrives through [`context::ContextModel`], which is
//! what lets a predictor be swapped without the search engine noticing.
//!
//! ## The two properties this crate is built around
//!
//! **Every operator produces a valid individual.** The genotype is a product
//! space over synonymous choices, so a codon-boundary operator cannot violate
//! the protein constraint — the constraint that matters most, and the one most
//! sequence-design MOEAs spend their complexity budget on. `cf-search` gets it
//! for free by never exposing a mutable view of a genotype's codons: the single
//! mutation path is [`genotype::Genotype::apply`], which delegates to
//! `CodonSequence::substitute` and inherits both the synonymy check
//! (`CF-SEQ-004`) and the locked-span check (`CF-SEQ-006`) from `P1-06`.
//!
//! **Every operator move is a landscape delta query.** [`genotype::Move`]
//! carries `cf_landscape_api::Edit` values, the same type `Target::Delta`
//! holds, so a move reaches the exact-delta path by being handed over rather
//! than by being reconstructed from a diff. Reconstructing would cost `O(L)`
//! per move and would put the query back on the expensive path by construction.
//!
//! Declared dependencies: cf-types, cf-landscape-api, cf-constraints,
//! cf-metrics.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod baseline;
pub mod context;
pub mod determinism;
pub mod generation;
pub mod genotype;
pub mod job;
pub mod operators;
pub mod population;
pub mod rng;
pub mod rationing;
pub mod tier1;
pub mod tier2;
pub mod tier3;

pub use baseline::{baseline_backtracking, BaselineBudget, BaselineOutcome, BaselineReport};
pub use context::{CodonWeights, ContextModel};
pub use determinism::{compare_splits, compare_thread_counts, with_threads, GATE_THREAD_COUNTS};
pub use generation::{
    run_generations, GenerationConfig, GenerationRecord, GenerationReport, GenerationRequest,
};
pub use genotype::{Genotype, Move, MoveKind};
pub use job::{run_segment, JobCheckpoint, SegmentReport, SegmentRequest};
pub use operators::{
    BlockGibbs, BlockUniform, CrossoverOperator, MutationContext, MutationOperator, OnePoint,
    OperatorMix, PottsGuided, TwoPoint, UnavailableOperator, UniformSynonymous, ViolationDirected,
};
pub use population::{Offspring, Population, SeedReport};
pub use rng::SearchRng;
pub use tier1::{tier1_seed, Tier1Config, Tier1Report, Tier1Request};
pub use tier2::{
    partition, tier2_block_ca, Block, Colour, Tier2Config, Tier2Plan, Tier2Report, Tier2Request,
};
