//! The resumable segment runner.
//!
//! Ticket P1-29, Rust half.
//!
//! ## Why the run is chopped into segments
//!
//! The job lifecycle belongs to `cf-orchestrator`, which is Python. The search
//! belongs here. ADR-0001 caps boundary crossings at one per generation, so a
//! Python loop calling into Rust once per generation is exactly at the limit
//! and a Python loop that checked a cancellation flag mid-generation would be
//! over it.
//!
//! So Rust advances *whole segments* — a run of complete generations — and
//! returns a state the caller can persist. The job layer decides between
//! segments whether to continue, checkpoint or stop. One crossing per
//! checkpoint interval is coarser than the limit, and the state it returns is
//! the checkpoint: there is no second representation to keep in sync.
//!
//! ## Why cancellation and checkpointing cannot interleave badly
//!
//! Both happen only between segments, and a segment is a whole number of
//! generations. So the persisted state is always a completed generation and
//! there is no such thing as a checkpoint taken mid-cancellation. That is a
//! property of where the boundary sits rather than of any locking.

use cf_constraints::{ConstraintChecker, TabuLedger};
use cf_types::{code, CfError, Codon, CodonSequence, Result};
use serde_json::{json, Value};

use crate::context::CodonWeights;
use crate::generation::{run_generations, GenerationConfig, GenerationRecord};
use crate::genotype::Genotype;
use crate::operators::OperatorMix;

/// A resumable point in a run.
///
/// Serializable in full: resuming from a checkpoint must be indistinguishable
/// from never having stopped, so everything that influences the next generation
/// has to be in here. The tabu ledger is the easy thing to forget — it carries
/// state across generations by design, and a resume that dropped it would
/// re-propose moves repair had just reverted, which is the thrash `P1-18b`
/// exists to break, restarted by the checkpoint mechanism itself.
#[derive(Clone, Debug)]
pub struct JobCheckpoint {
    /// Generations already completed.
    pub generation: u32,
    pub incumbent: CodonSequence,
    /// Live tabu entries as `(region, codon_index, codon, remaining)`.
    pub tabu: Vec<(usize, usize, Codon, u32)>,
    pub repairs_attempted: u64,
    pub repairs_succeeded: u64,
    /// Every generation so far, oldest first.
    pub history: Vec<GenerationRecord>,
}

impl JobCheckpoint {
    pub fn start(incumbent: &Genotype) -> Self {
        Self {
            generation: 0,
            incumbent: incumbent.sequence().clone(),
            tabu: Vec::new(),
            repairs_attempted: 0,
            repairs_succeeded: 0,
            history: Vec::new(),
        }
    }

    pub fn genotype(&self) -> Genotype {
        Genotype::new(self.incumbent.clone())
    }

    pub fn to_json(&self) -> Result<Value> {
        Ok(json!({
            "generation": self.generation,
            "incumbent": self.incumbent.to_json()?,
            "tabu": self.tabu.iter().map(|(region, index, codon, remaining)| json!({
                "region": region,
                "codon_index": index,
                "codon": codon.to_string_upper(),
                "remaining": remaining,
            })).collect::<Vec<_>>(),
            "repairs_attempted": self.repairs_attempted,
            "repairs_succeeded": self.repairs_succeeded,
            "history": self.history.iter().map(|record| json!({
                "generation": record.generation,
                "proposed": record.proposed,
                "needed_repair": record.needed_repair,
                "repaired": record.repaired,
                "tabu_live": record.tabu_live,
                "feasible": record.feasible,
            })).collect::<Vec<_>>(),
        }))
    }

    pub fn from_json(value: &Value) -> Result<Self> {
        let bad = |what: &str| {
            CfError::new(
                code::ORCH_CHECKPOINT_INVALID,
                format!("checkpoint is missing or malformed: {what}"),
            )
        };
        let map = value.as_object().ok_or_else(|| bad("not an object"))?;

        let mut tabu = Vec::new();
        for entry in map
            .get("tabu")
            .and_then(Value::as_array)
            .ok_or_else(|| bad("tabu"))?
        {
            let entry = entry.as_object().ok_or_else(|| bad("tabu entry"))?;
            tabu.push((
                entry
                    .get("region")
                    .and_then(Value::as_u64)
                    .ok_or_else(|| bad("region"))? as usize,
                entry
                    .get("codon_index")
                    .and_then(Value::as_u64)
                    .ok_or_else(|| bad("codon_index"))? as usize,
                Codon::parse(
                    entry
                        .get("codon")
                        .and_then(Value::as_str)
                        .ok_or_else(|| bad("codon"))?,
                )?,
                entry
                    .get("remaining")
                    .and_then(Value::as_u64)
                    .ok_or_else(|| bad("remaining"))? as u32,
            ));
        }

        let mut history = Vec::new();
        for record in map
            .get("history")
            .and_then(Value::as_array)
            .ok_or_else(|| bad("history"))?
        {
            let record = record.as_object().ok_or_else(|| bad("history entry"))?;
            history.push(GenerationRecord {
                generation: record
                    .get("generation")
                    .and_then(Value::as_u64)
                    .ok_or_else(|| bad("generation"))? as u32,
                proposed: record
                    .get("proposed")
                    .and_then(Value::as_array)
                    .ok_or_else(|| bad("proposed"))?
                    .iter()
                    .filter_map(|v| v.as_u64().map(|n| n as usize))
                    .collect(),
                needed_repair: record
                    .get("needed_repair")
                    .and_then(Value::as_bool)
                    .unwrap_or(false),
                repaired: record
                    .get("repaired")
                    .and_then(Value::as_array)
                    .map(|a| {
                        a.iter()
                            .filter_map(|v| v.as_u64().map(|n| n as usize))
                            .collect()
                    })
                    .unwrap_or_default(),
                tabu_live: record.get("tabu_live").and_then(Value::as_u64).unwrap_or(0) as usize,
                feasible: record
                    .get("feasible")
                    .and_then(Value::as_bool)
                    .unwrap_or(true),
            });
        }

        Ok(Self {
            generation: map
                .get("generation")
                .and_then(Value::as_u64)
                .ok_or_else(|| bad("generation"))? as u32,
            incumbent: CodonSequence::from_json(
                map.get("incumbent").ok_or_else(|| bad("incumbent"))?,
            )?,
            tabu,
            repairs_attempted: map
                .get("repairs_attempted")
                .and_then(Value::as_u64)
                .unwrap_or(0),
            repairs_succeeded: map
                .get("repairs_succeeded")
                .and_then(Value::as_u64)
                .unwrap_or(0),
            history,
        })
    }
}

/// What one segment produced.
#[derive(Clone, Debug)]
pub struct SegmentReport {
    pub checkpoint: JobCheckpoint,
    /// Warnings from the operator mix. Collected here, delivered by the job
    /// layer into `job_status.warnings`.
    pub warnings: Vec<CfError>,
    pub generations_run: u32,
}

/// Advance `generations` complete generations from a checkpoint.
///
/// The generation index continues from the checkpoint, and each generation's
/// RNG stream is derived from `(root_seed, absolute generation index)` — not
/// from a counter that restarts. That is what makes a resumed run identical to
/// an uninterrupted one: generation `k` sees the same stream whether it is the
/// first generation of a fresh segment or the fiftieth of a long one.
#[derive(Clone, Copy, Debug)]
pub struct SegmentRequest<'a> {
    pub checker: &'a ConstraintChecker,
    pub mix: &'a OperatorMix,
    pub weights: &'a CodonWeights,
    pub leak: f64,
    pub region: usize,
    pub masking: cf_constraints::KmerMasking,
    pub root_seed: u64,
}

pub fn run_segment(
    checkpoint: &JobCheckpoint,
    request: &SegmentRequest<'_>,
    generations: u32,
) -> Result<SegmentReport> {
    let SegmentRequest {
        checker,
        mix,
        weights,
        leak,
        region,
        masking,
        root_seed,
    } = *request;

    if generations == 0 {
        return Ok(SegmentReport {
            checkpoint: checkpoint.clone(),
            warnings: mix.warnings(),
            generations_run: 0,
        });
    }

    let mut ledger = TabuLedger::new(leak);
    for (region, index, codon, remaining) in &checkpoint.tabu {
        ledger.record(*region, *index, *codon, *remaining);
    }

    let config = GenerationConfig {
        generations,
        tabu_leak: leak,
        region,
        masking,
        first_generation: checkpoint.generation,
    };
    let start = checkpoint.genotype();
    let report = run_generations(
        &crate::generation::GenerationRequest {
            start: &start,
            checker,
            mix,
            weights,
            config: &config,
            root_seed,
        },
        Some(ledger),
    )?;

    let mut history = checkpoint.history.clone();
    history.extend(report.history.iter().cloned());

    Ok(SegmentReport {
        checkpoint: JobCheckpoint {
            generation: checkpoint.generation + generations,
            incumbent: report.incumbent.sequence().clone(),
            tabu: report.tabu_snapshot,
            repairs_attempted: checkpoint.repairs_attempted + report.repairs_attempted,
            repairs_succeeded: checkpoint.repairs_succeeded + report.repairs_succeeded,
            history,
        },
        warnings: report.warnings,
        generations_run: generations,
    })
}
