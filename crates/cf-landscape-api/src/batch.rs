//! The batch-yield protocol: `step` / `absorb` / `snapshot` / `resume`.
//!
//! Ticket P2-05. Specification: `03-nsga-loop-and-operators.md` §6.
//!
//! ## What this is for
//!
//! An L0 term is one that cannot be scored incrementally — a language model, a
//! global fold — and §6 measures it as roughly 95 % of tier 3's wall time. It is
//! rationed: a bounded `E_max` individuals per generation get an exact score and
//! the rest are `UNKNOWN`. The scoring happens in **one** call per generation,
//! because ADR-0001 caps boundary crossings at one per generation and because a
//! batched model is the only affordable way to run one.
//!
//! So the search must be able to stop, hand out a batch, and continue when the
//! answers arrive. That is this state machine.
//!
//! ## It is synchronous, and that is not a limitation being papered over
//!
//! `step` yields a batch and `absorb` takes the results; nothing runs
//! concurrently with the search. §6's pseudocode is explicit —
//! `results <- WORKER.SCORE_BATCH(chosen)` blocks — and the batching exists to
//! amortize *model* cost over many sequences, not to overlap work. The protocol
//! therefore needs no worker pool, which is the capability `P1-29` declined to
//! build and `CF-CLI-003` still reports as absent.
//!
//! What the protocol *does* buy is that the evaluator may live anywhere: in
//! process, in a subprocess, on a GPU host. `DEP-LAND-004` records that the
//! out-of-process transport is deferred — the seam is shaped so it can be added
//! without changing the protocol, and adding a subprocess and an Arrow
//! dependency whose only consumer is a stub would be building the plumbing
//! before the thing that needs it exists.

use std::collections::BTreeSet;

use cf_types::{code, CfError, CodonSequence, Result, TermId};

/// A batch handed out for scoring.
#[derive(Clone, Debug)]
pub struct Batch {
    /// Which generation asked. Carried so a late or duplicated `absorb` can be
    /// recognised rather than silently applied to the wrong generation.
    pub generation: u32,
    pub term: TermId,
    /// The sequences to score, in a fixed order that `absorb` must match.
    pub items: Vec<CodonSequence>,
}

impl Batch {
    pub fn len(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
}

/// Scores returned for a batch, in the batch's order.
#[derive(Clone, Debug, PartialEq)]
pub struct BatchResults {
    pub generation: u32,
    /// One entry per batch item. `None` is `UNKNOWN` — a scorer that could not
    /// produce a value says so rather than returning a placeholder.
    pub values: Vec<Option<f64>>,
}

/// What the protocol wants next.
#[derive(Clone, Debug, PartialEq)]
pub enum Yield {
    /// Score these, then call `absorb`.
    Score(u32, usize),
    /// Nothing outstanding; the caller may advance.
    Ready,
}

/// Scores a whole batch in one call.
///
/// The single method is the point: a trait with a per-item `score` would invite
/// a caller to loop, which is the crossing pattern ADR-0001 forbids and the
/// cost §6 exists to avoid.
pub trait BatchScorer {
    fn term(&self) -> &TermId;

    /// Score every item. Must return one value per item, in order.
    ///
    /// `None` for an item the scorer could not evaluate. `TST-LM-001` requires
    /// this to agree with a per-item `score_one` for every batch size, which is
    /// what makes a batched implementation checkable against an obvious one.
    fn score_batch(&self, items: &[CodonSequence]) -> Result<Vec<Option<f64>>>;

    /// Score one item. The reference a batch is checked against.
    fn score_one(&self, item: &CodonSequence) -> Result<Option<f64>> {
        Ok(self
            .score_batch(std::slice::from_ref(item))?
            .pop()
            .flatten())
    }
}

/// The resumable protocol state.
///
/// Deliberately holds no scorer and no sequences beyond the outstanding batch:
/// everything it carries has to survive a `snapshot`/`resume` pair, and the
/// less that is, the less there is to lose. The rule that anything the reported
/// result depends on must travel in the checkpoint has been learned four times
/// in this repository; the cheapest way to obey it is to depend on little.
#[derive(Clone, Debug, PartialEq)]
pub struct BatchProtocol {
    generation: u32,
    /// `Some` while a batch is outstanding.
    outstanding: Option<usize>,
    /// Scores absorbed for the current generation, in batch order.
    absorbed: Vec<Option<f64>>,
    /// Exact evaluations spent, across every generation.
    evaluations: u64,
    /// Items deferred because the batch was truncated to `E_max`.
    ///
    /// §6's component (d): they go to the head of the next generation's queue
    /// rather than being forgotten, so a bounded budget delays an evaluation
    /// instead of dropping it.
    deferred: Vec<usize>,
    /// Every item that has ever been deferred, minus those later evaluated.
    ///
    /// The visible trace behind `CF-SRCH-016`. Deferring correctly is not the
    /// same as saying so: an item the budget never reached stays `UNKNOWN` to
    /// the end, which means the run searched less than its plan describes and
    /// nothing in the output admits it. That is structurally what `CF-ORCH-006`
    /// exists to surface for a skipped tier.
    ///
    /// Lives here rather than in `cf-search`'s rationer because this type owns
    /// the deferral queue; tracking it in both places would be two mechanisms
    /// for one rule (`DEP-SRCH-011`).
    stranded: BTreeSet<usize>,
    /// Batch order of the outstanding batch, so `absorb` knows which items its
    /// values belong to.
    in_flight: Vec<usize>,
}

impl Default for BatchProtocol {
    fn default() -> Self {
        Self::new()
    }
}

impl BatchProtocol {
    pub fn new() -> Self {
        Self {
            generation: 0,
            outstanding: None,
            absorbed: Vec::new(),
            evaluations: 0,
            deferred: Vec::new(),
            stranded: BTreeSet::new(),
            in_flight: Vec::new(),
        }
    }

    pub fn generation(&self) -> u32 {
        self.generation
    }

    pub fn evaluations(&self) -> u64 {
        self.evaluations
    }

    pub fn deferred(&self) -> &[usize] {
        &self.deferred
    }

    /// Items deferred at some point and never exactly evaluated since.
    ///
    /// The orchestrator turns a non-empty result into a `CF-SRCH-016` warning
    /// with stage and subject provenance.
    pub fn stranded(&self) -> Vec<usize> {
        self.stranded.iter().copied().collect()
    }

    pub fn is_waiting(&self) -> bool {
        self.outstanding.is_some()
    }

    /// Choose a batch and yield it.
    ///
    /// `candidates` are indices into the caller's population, already ranked by
    /// whatever selection §6 specifies; this takes the first `budget` of them
    /// with the previously-deferred ones first. Truncation defers rather than
    /// discards.
    pub fn step(&mut self, candidates: &[usize], budget: usize) -> Result<(Yield, Vec<usize>)> {
        if self.outstanding.is_some() {
            return Err(CfError::new(
                code::LAND_BATCH_OUTSTANDING,
                format!(
                    "generation {} already has a batch outstanding; absorb it before \
                     stepping again, or the results will be matched to the wrong batch",
                    self.generation
                ),
            ));
        }

        // Deferred first, then new candidates, with no repeats. An index that
        // was deferred and is also proposed again must not be scored twice.
        let mut queue: Vec<usize> = Vec::with_capacity(self.deferred.len() + candidates.len());
        for index in self.deferred.iter().chain(candidates.iter()) {
            if !queue.contains(index) {
                queue.push(*index);
            }
        }

        let chosen: Vec<usize> = queue.iter().copied().take(budget).collect();
        self.deferred = queue.into_iter().skip(budget).collect();

        if chosen.is_empty() {
            self.outstanding = None;
            return Ok((Yield::Ready, chosen));
        }
        for index in &self.deferred {
            self.stranded.insert(*index);
        }
        self.outstanding = Some(chosen.len());
        self.in_flight = chosen.clone();
        self.absorbed.clear();
        Ok((Yield::Score(self.generation, chosen.len()), chosen))
    }

    /// Take the results of the outstanding batch.
    ///
    /// Refuses a mismatched generation or length rather than aligning by
    /// position and hoping. A batch protocol that silently accepted the wrong
    /// results would attribute one sequence's score to another, and nothing
    /// downstream could detect it — the scores are plausible, just wrong.
    pub fn absorb(&mut self, results: &BatchResults) -> Result<()> {
        let expected = self.outstanding.ok_or_else(|| {
            CfError::new(
                code::LAND_BATCH_UNEXPECTED,
                "no batch is outstanding, so these results belong to a batch that was \
                 already absorbed or never issued",
            )
        })?;
        if results.generation != self.generation {
            return Err(CfError::new(
                code::LAND_BATCH_UNEXPECTED,
                format!(
                    "results are for generation {} but generation {} is waiting",
                    results.generation, self.generation
                ),
            ));
        }
        if results.values.len() != expected {
            return Err(CfError::new(
                code::LAND_BATCH_UNEXPECTED,
                format!(
                    "batch had {expected} items and {} results arrived; aligning them by \
                     position would attribute one sequence's score to another",
                    results.values.len()
                ),
            ));
        }

        self.evaluations += results.values.iter().filter(|v| v.is_some()).count() as u64;
        // An item is stranded only until it is actually scored. A `None` value
        // is not a score: the scorer said `UNKNOWN`, so the item is still
        // waiting for a number it never got.
        for (index, value) in self.in_flight.iter().zip(&results.values) {
            if value.is_some() {
                self.stranded.remove(index);
            }
        }
        self.absorbed = results.values.clone();
        self.outstanding = None;
        Ok(())
    }

    /// The scores absorbed for this generation, in batch order.
    pub fn absorbed(&self) -> &[Option<f64>] {
        &self.absorbed
    }

    /// Advance to the next generation.
    pub fn advance(&mut self) -> Result<()> {
        if self.outstanding.is_some() {
            return Err(CfError::new(
                code::LAND_BATCH_OUTSTANDING,
                "cannot advance with a batch outstanding; its results would be absorbed \
                 into the wrong generation",
            ));
        }
        self.generation += 1;
        self.absorbed.clear();
        Ok(())
    }

    /// A serializable snapshot.
    ///
    /// Includes the outstanding batch size, so a protocol snapshotted mid-batch
    /// resumes still waiting rather than resuming as if the batch had been
    /// answered. That is the case a checkpoint is most likely to be taken in —
    /// the search is idle precisely while the scorer runs.
    pub fn snapshot(&self) -> serde_json::Value {
        serde_json::json!({
            "generation": self.generation,
            "outstanding": self.outstanding,
            "absorbed": self.absorbed,
            "evaluations": self.evaluations,
            "deferred": self.deferred,
            // Both travel, because both decide what the run reports. Dropping
            // `stranded` would make CF-SRCH-016 vanish across a resume, which
            // is the one place a warning about silence must not itself go
            // silent. Dropping `in_flight` would leave a mid-batch resume
            // unable to clear the right items on absorb.
            "stranded": self.stranded.iter().copied().collect::<Vec<usize>>(),
            "in_flight": self.in_flight,
        })
    }

    pub fn resume(document: &serde_json::Value) -> Result<Self> {
        let bad = |what: &str| {
            CfError::new(
                code::LAND_BATCH_UNEXPECTED,
                format!("batch protocol snapshot is missing or malformed: {what}"),
            )
        };
        let map = document.as_object().ok_or_else(|| bad("not an object"))?;
        Ok(Self {
            generation: map
                .get("generation")
                .and_then(serde_json::Value::as_u64)
                .ok_or_else(|| bad("generation"))? as u32,
            outstanding: map
                .get("outstanding")
                .and_then(|v| {
                    if v.is_null() {
                        Some(None)
                    } else {
                        v.as_u64().map(|n| Some(n as usize))
                    }
                })
                .ok_or_else(|| bad("outstanding"))?,
            absorbed: map
                .get("absorbed")
                .and_then(serde_json::Value::as_array)
                .ok_or_else(|| bad("absorbed"))?
                .iter()
                .map(serde_json::Value::as_f64)
                .collect(),
            evaluations: map
                .get("evaluations")
                .and_then(serde_json::Value::as_u64)
                .ok_or_else(|| bad("evaluations"))?,
            deferred: map
                .get("deferred")
                .and_then(serde_json::Value::as_array)
                .ok_or_else(|| bad("deferred"))?
                .iter()
                .filter_map(|v| v.as_u64().map(|n| n as usize))
                .collect(),
            stranded: map
                .get("stranded")
                .and_then(serde_json::Value::as_array)
                .map(|items| items.iter().filter_map(serde_json::Value::as_u64).map(|n| n as usize).collect())
                .unwrap_or_default(),
            in_flight: map
                .get("in_flight")
                .and_then(serde_json::Value::as_array)
                .map(|items| items.iter().filter_map(serde_json::Value::as_u64).map(|n| n as usize).collect())
                .unwrap_or_default(),
        })
    }
}
