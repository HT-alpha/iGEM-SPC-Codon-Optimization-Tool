//! One feasible seed after tier 1 refuses a constraint the DP cannot carry.
//!
//! The in-loop repair operator (`repair`) is bounded so it cannot starve a
//! generation, and giving up there is correct: an individual that cannot be
//! cheaply repaired stays infeasible. That bound is the wrong behaviour at
//! the moment tier 1 has produced no seed at all. GC windows and scattered
//! k-mer repeats are post-checked, not carried in the DP, so every candidate
//! can fail `CF-SRCH-006` on a constraint set that still has a synonymous
//! solution. This pass runs once, before tier 2, and either returns a
//! sequence the checker accepts or names the rules that remain.
//!
//! It is deterministic. Ties break by higher relative adaptiveness, then by
//! the smaller codon index. The in-loop operator's random tie-break exists to
//! keep a population diverse; a single seed has nothing to diversify.
//!
//! `DEP-SRCH-015`.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering as AtomicOrdering};
use std::time::Instant;

use cf_types::codon::Codon;

/// Inclusive timings for `CF_SEED_PROFILE=1`. Off the hot path when unset:
/// each site loads one relaxed flag and returns.
struct SeedCounter {
    ns: AtomicU64,
    calls: AtomicU64,
}

impl SeedCounter {
    const fn new() -> Self {
        Self {
            ns: AtomicU64::new(0),
            calls: AtomicU64::new(0),
        }
    }

    fn record(&self, started: Instant) {
        self.ns
            .fetch_add(started.elapsed().as_nanos() as u64, AtomicOrdering::Relaxed);
        self.calls.fetch_add(1, AtomicOrdering::Relaxed);
    }
}

struct SeedProfile {
    enabled: AtomicBool,
    gc_sweeps: SeedCounter,
    fix_patterns: SeedCounter,
    fix_kmers: SeedCounter,
    fix_kmer_pairs: SeedCounter,
    polish: SeedCounter,
    ascend: SeedCounter,
    climb: SeedCounter,
    homo: SeedCounter,
    gc_delta: SeedCounter,
    prefix: SeedCounter,
    pattern_count: SeedCounter,
    pattern_violations: SeedCounter,
    is_feasible: SeedCounter,
    hot_codons: SeedCounter,
    kmer_replace: SeedCounter,
    checker_check: SeedCounter,
}

static SEED_PROFILE: SeedProfile = SeedProfile {
    enabled: AtomicBool::new(false),
    gc_sweeps: SeedCounter::new(),
    fix_patterns: SeedCounter::new(),
    fix_kmers: SeedCounter::new(),
    fix_kmer_pairs: SeedCounter::new(),
    polish: SeedCounter::new(),
    ascend: SeedCounter::new(),
    climb: SeedCounter::new(),
    homo: SeedCounter::new(),
    gc_delta: SeedCounter::new(),
    prefix: SeedCounter::new(),
    pattern_count: SeedCounter::new(),
    pattern_violations: SeedCounter::new(),
    is_feasible: SeedCounter::new(),
    hot_codons: SeedCounter::new(),
    kmer_replace: SeedCounter::new(),
    checker_check: SeedCounter::new(),
};

fn seed_profile_on() -> bool {
    SEED_PROFILE.enabled.load(AtomicOrdering::Relaxed)
}

fn arm_seed_profile() {
    if std::env::var_os("CF_SEED_PROFILE").is_some() {
        SEED_PROFILE.enabled.store(true, AtomicOrdering::Relaxed);
    }
}

fn dump_seed_profile() {
    if !seed_profile_on() {
        return;
    }
    let row = |name: &str, counter: &SeedCounter| {
        let calls = counter.calls.load(AtomicOrdering::Relaxed);
        let ns = counter.ns.load(AtomicOrdering::Relaxed);
        format!(
            "{{\"name\":\"{name}\",\"calls\":{calls},\"seconds\":{:.6}}}",
            ns as f64 / 1e9
        )
    };
    eprintln!(
        "CF_SEED_PROFILE {}",
        [
            row("gc_sweeps", &SEED_PROFILE.gc_sweeps),
            row("fix_patterns", &SEED_PROFILE.fix_patterns),
            row("fix_kmers", &SEED_PROFILE.fix_kmers),
            row("fix_kmer_pairs", &SEED_PROFILE.fix_kmer_pairs),
            row("polish", &SEED_PROFILE.polish),
            row("ascend", &SEED_PROFILE.ascend),
            row("climb", &SEED_PROFILE.climb),
            row("homo_touching", &SEED_PROFILE.homo),
            row("gc_delta", &SEED_PROFILE.gc_delta),
            row("prefix_update", &SEED_PROFILE.prefix),
            row("pattern_count", &SEED_PROFILE.pattern_count),
            row("pattern_violations", &SEED_PROFILE.pattern_violations),
            row("is_feasible", &SEED_PROFILE.is_feasible),
            row("hot_codons", &SEED_PROFILE.hot_codons),
            row("kmer_replace", &SEED_PROFILE.kmer_replace),
            row("checker_check", &SEED_PROFILE.checker_check),
        ]
        .join(",")
    );
}
use cf_types::error::Result;
use cf_types::seqrec::CodonSequence;

use crate::detect::{ConstraintChecker, ConstructView};
use crate::spec::{ConstraintClass, Rule};
use crate::spend::SpendLimit;

const MAX_CYCLES: u32 = 8;
/// Joint GC and k-mer pair attempts after the single-window escape still
/// leaves the seed infeasible. The same bound as that escape: four tries,
/// not a count fitted to one gene. `DEP-SRCH-024`.
const JOINT_ATTEMPTS: u32 = 4;
/// Codons taken from the stuck GC window, spread from the first codon to the
/// last so the edges are included.
const JOINT_GC_CODONS: usize = 8;
/// Codons taken from over-limit k-mers, spread across those sites.
const JOINT_KMER_CODONS: usize = 8;
/// The pair search runs on at most this many sites.
const JOINT_SITES: usize = 12;
const GC_SWEEPS: u32 = 8;
const PATTERN_ROUNDS: u32 = 4;
/// Pair lifts after a feasible seed. One codon at a time cannot trade a
/// higher-weight synonym for the neighbour that clears the repeat it creates.
/// `DEP-SRCH-018`.
const ASCENT_PASSES: u32 = 2;
const ASCENT_POSITIONS: usize = 24;
const ASCENT_PARTNERS: usize = 6;
const ASCENT_NEIGHBOUR: usize = 4;
/// Conflict climbs after the local pairs. A repeated 12-mer's other copy, and
/// a codon that can give a GC window back a base, are often outside the four
/// codons `DEP-SRCH-018` searches. `DEP-SRCH-021`.
const CLIMB_PASSES: u32 = 6;
const CLIMB_SITES: usize = 48;
const CLIMB_PAIR_SITES: usize = 6;
const CLIMB_PARTNER_TRIALS: usize = 48;
const KMER_ROUNDS: u32 = 12;

/// One constraint the checker still reports after seed resolution.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedResidual {
    pub constraint_id: String,
    pub kind: String,
    pub count: usize,
}

/// A rule this version does not evaluate. Named so a caller cannot present it
/// as satisfied.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnsearchedRule {
    pub id: String,
    pub reason: String,
}

/// A position whose codon is not the most frequent synonym of the loaded table.
#[derive(Clone, Debug, PartialEq)]
pub struct BaselineDeviation {
    pub codon_index: usize,
    pub codon: String,
    pub weight: f64,
    pub baseline: String,
    pub baseline_weight: f64,
    /// Constraint ids whose actionable blocking count rises if this position
    /// is put back to `baseline` on its own. Empty when that single edit does
    /// not add a violation the checker reports.
    pub constraints_if_reverted: Vec<String>,
}

/// The sequence seed resolution ended on, and what it could not clear.
#[derive(Clone)]
pub struct SeedReport {
    pub sequence: CodonSequence,
    pub feasible: bool,
    pub cycles: u32,
    pub edits: u32,
    pub residual: Vec<SeedResidual>,
    pub unsearched: Vec<UnsearchedRule>,
}

impl std::fmt::Debug for SeedReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SeedReport")
            .field("feasible", &self.feasible)
            .field("cycles", &self.cycles)
            .field("edits", &self.edits)
            .field("residual", &self.residual)
            .field("unsearched", &self.unsearched)
            .field("cds_len", &self.sequence.cds_string().len())
            .finish()
    }
}

/// Resolve one synonymous seed.
///
/// Scores four deterministic attempts: two tie orders (higher relative
/// adaptiveness, then lower codon index) and, within each, two starts (the
/// most-frequent-codon assignment, then the input). A feasible result
/// outranks an infeasible one. Two feasible results keep the higher CAI on
/// `weights` (stops and single-codon families excluded; a counted weight of
/// 0 scores 0), then the fewer committed edits. The first feasible attempt
/// does not end the search: on the assembly set that attempt is not always
/// the highest CAI. `DEP-SRCH-022`.
pub fn resolve_feasible_seed(
    sequence: &CodonSequence,
    checker: &ConstraintChecker,
    weights: &[f64; 64],
) -> Result<SeedReport> {
    resolve_feasible_seed_limited(sequence, checker, weights, SpendLimit::default())
}

/// [`resolve_feasible_seed`] with a wall-clock and resident-memory cap.
///
/// `DEP-BUDG-002`. A default limit is unlimited. A deadline that has already
/// passed returns `CF-BUDGET-003` before a cycle edits the sequence.
pub fn resolve_feasible_seed_limited(
    sequence: &CodonSequence,
    checker: &ConstraintChecker,
    weights: &[f64; 64],
    limits: SpendLimit,
) -> Result<SeedReport> {
    arm_seed_profile();
    let unsearched = unsearched_rules(checker);
    let mut best: Option<Attempt> = None;
    for tie in [TieBreak::Weight, TieBreak::LowIndex] {
        for start in [most_frequent(sequence, weights)?, sequence.clone()] {
            let mut solver = Solver::new(start, checker, weights, tie)?.with_limits(limits);
            let cycles = solver.run()?;
            let feasible = solver.is_feasible();
            let cai = cai_of(&solver.seq, weights);
            let attempt = Attempt {
                seq: solver.seq,
                feasible,
                cycles,
                edits: solver.edits,
                cai,
            };
            let improved = match &best {
                None => true,
                Some(current) => attempt.beats(current),
            };
            if improved {
                best = Some(attempt);
            }
        }
    }
    let attempt = best.expect("four starts were attempted");
    let solver = Solver::new(attempt.seq, checker, weights, TieBreak::Weight)?.with_limits(limits);
    let residual = if attempt.feasible {
        Vec::new()
    } else {
        solver.residuals()
    };
    let report = SeedReport {
        sequence: solver.seq,
        feasible: attempt.feasible,
        cycles: attempt.cycles,
        edits: attempt.edits,
        residual,
        unsearched,
    };
    dump_seed_profile();
    Ok(report)
}

/// Positions that differ from the most-frequent synonym, and which checked
/// rules that synonym would break on its own.
///
/// One [`Solver`] carries the counts. Putting the preferred codon in and
/// taking it back uses `set_codon`, which is the same edit the seed applies,
/// and the tally below is the checker's blocking actionable counts:
///
/// * pattern ids, one per ledger hit (or per `check_patterns` hit when a
///   blocking regular expression is present, because that match is not
///   confined to the motif window);
/// * one per k-mer key that exceeds its limit and touches an editable span,
///   which is one checker violation, not `count - max`;
/// * one per GC window outside its band that overlaps an editable span.
///
/// A missing id is a count of zero on both sides. The CDS is not what this
/// function returns; only the sentence that names those ids.
pub fn explain_deviations(
    sequence: &CodonSequence,
    checker: &ConstraintChecker,
    weights: &[f64; 64],
) -> Result<Vec<BaselineDeviation>> {
    let baseline = most_frequent(sequence, weights)?;
    let mut solver = Solver::new(sequence.clone(), checker, weights, TieBreak::Weight)?;
    let before = solver.actionable_tally();
    let mut out = Vec::new();
    for index in 0..solver.seq.codon_count() {
        if solver.seq.is_locked(index) {
            continue;
        }
        let codon = solver.seq.codons()[index];
        let preferred = baseline.codons()[index];
        if codon == preferred {
            continue;
        }
        solver.set_codon(index, preferred)?;
        let after = solver.actionable_tally();
        solver.set_codon(index, codon)?;
        let mut constraints_if_reverted: Vec<String> = after
            .iter()
            .filter(|(id, count)| **count > *before.get(*id).unwrap_or(&0))
            .map(|(id, _)| (*id).clone())
            .collect();
        constraints_if_reverted.sort();
        out.push(BaselineDeviation {
            codon_index: index,
            codon: codon.to_string_upper(),
            weight: weights[codon.index_usize()],
            baseline: preferred.to_string_upper(),
            baseline_weight: weights[preferred.index_usize()],
            constraints_if_reverted,
        });
    }
    Ok(out)
}

#[derive(Clone, Copy)]
enum TieBreak {
    /// Higher relative adaptiveness, then the smaller codon index.
    Weight,
    /// Smaller codon index. A second pass when the weight order gets stuck.
    LowIndex,
}

struct Attempt {
    seq: CodonSequence,
    feasible: bool,
    cycles: u32,
    edits: u32,
    cai: f64,
}

struct ClimbChoice {
    cai: f64,
    primary: Codon,
    partners: Vec<(usize, Codon)>,
}

impl ClimbChoice {
    fn beats(&self, other: &Self) -> bool {
        match self.cai.total_cmp(&other.cai) {
            Ordering::Greater => true,
            Ordering::Less => false,
            Ordering::Equal => match self.partners.len().cmp(&other.partners.len()) {
                Ordering::Less => true,
                Ordering::Greater => false,
                Ordering::Equal => {
                    let left = self
                        .partners
                        .iter()
                        .map(|(site, codon)| (*site, codon.index()))
                        .chain(std::iter::once((usize::MAX, self.primary.index())));
                    let right = other
                        .partners
                        .iter()
                        .map(|(site, codon)| (*site, codon.index()))
                        .chain(std::iter::once((usize::MAX, other.primary.index())));
                    left.cmp(right) == Ordering::Less
                }
            },
        }
    }
}

impl Attempt {
    fn beats(&self, other: &Self) -> bool {
        match (self.feasible, other.feasible) {
            (true, false) => true,
            (false, true) => false,
            _ => match self.cai.total_cmp(&other.cai) {
                Ordering::Greater => true,
                Ordering::Less => false,
                Ordering::Equal => self.edits < other.edits,
            },
        }
    }
}

fn most_frequent(sequence: &CodonSequence, weights: &[f64; 64]) -> Result<CodonSequence> {
    let mut next = sequence.clone();
    for index in 0..next.codon_count() {
        if next.is_locked(index) {
            continue;
        }
        let mut best = next.codons()[index];
        for &synonym in next.synonyms_at(index) {
            if prefer_weight(weights, synonym, best) {
                best = synonym;
            }
        }
        if best != next.codons()[index] {
            next.substitute(index, best)?;
        }
    }
    Ok(next)
}

fn prefer_weight(weights: &[f64; 64], candidate: Codon, incumbent: Codon) -> bool {
    match weights[candidate.index_usize()].total_cmp(&weights[incumbent.index_usize()]) {
        Ordering::Greater => true,
        Ordering::Less => false,
        Ordering::Equal => candidate.index() < incumbent.index(),
    }
}

fn cai_of(sequence: &CodonSequence, weights: &[f64; 64]) -> f64 {
    let mut sum = 0.0;
    let mut counted = 0usize;
    for (index, codon) in sequence.codons().iter().copied().enumerate() {
        if sequence.table().is_stop(codon) || sequence.synonyms_at(index).len() <= 1 {
            continue;
        }
        let weight = weights[codon.index_usize()];
        if weight == 0.0 {
            return 0.0;
        }
        sum += weight.ln();
        counted += 1;
    }
    if counted == 0 {
        1.0
    } else {
        (sum / counted as f64).exp()
    }
}

struct GcBand {
    id: String,
    window: usize,
    min: f64,
    max: f64,
    /// Windows of this band that sit outside `[min, max]` and overlap an
    /// editable span. That is the checker's actionable GC violation count
    /// for this rule, one window per violation.
    offenders: u32,
}

struct KmerCfg {
    id: String,
    k: usize,
    max_occurrences: u32,
    canonical: bool,
}

struct KmerBook {
    rules: Vec<KmerCfg>,
    counts: Vec<HashMap<u64, u32>>,
    /// Occurrences of each key whose span overlaps an editable nucleotide.
    ///
    /// The checker's k-mer violation is actionable when any occurrence touches
    /// an editable span. A repeat that lives entirely in a flank or a locked
    /// codon still adds to [`Self::excess`] — repair uses that total — but it
    /// does not add to [`Self::actionable_excess`].
    touching: Vec<HashMap<u64, u32>>,
    excess: Vec<u64>,
    /// Sum of `count - max` over keys that exceed the limit and have at least
    /// one editable occurrence. Zero exactly when `ConstraintChecker::check`
    /// would report no actionable blocking k-mer violation for these rules.
    actionable_excess: Vec<u64>,
    /// Number of keys that are an actionable violation. The checker emits one
    /// violation per such key, so this is the count `explain_deviations`
    /// compares. It is not [`Self::actionable_excess`]: a key that is over the
    /// limit by five still contributes one.
    offenders: Vec<u64>,
    editable: Vec<(usize, usize)>,
}

impl KmerBook {
    fn build(nt: &[u8], rules: Vec<KmerCfg>, editable: &[(usize, usize)]) -> Self {
        let mut book = Self {
            counts: vec![HashMap::new(); rules.len()],
            touching: vec![HashMap::new(); rules.len()],
            excess: vec![0; rules.len()],
            actionable_excess: vec![0; rules.len()],
            offenders: vec![0; rules.len()],
            editable: editable.to_vec(),
            rules,
        };
        let rule_count = book.rules.len();
        for rule_index in 0..rule_count {
            let k = book.rules[rule_index].k;
            let canonical = book.rules[rule_index].canonical;
            if nt.len() < k {
                continue;
            }
            for pos in 0..=nt.len() - k {
                if let Some(key) = pack(nt, pos, k, canonical) {
                    let touches = spans_overlap(&book.editable, pos, pos + k);
                    book.bump(rule_index, key, 1, touches);
                }
            }
        }
        book
    }

    fn excess(&self) -> u64 {
        self.excess.iter().sum()
    }

    fn actionable_excess(&self) -> u64 {
        self.actionable_excess.iter().sum()
    }

    fn bump(&mut self, rule: usize, key: u64, delta: i32, touches: bool) {
        let max = self.rules[rule].max_occurrences;
        let before = self.counts[rule].get(&key).copied().unwrap_or(0);
        let before_touch = self.touching[rule].get(&key).copied().unwrap_or(0);
        let before_over = u64::from(before.saturating_sub(max));
        let before_act = actionable_over(before, before_touch, max);
        let before_off = u64::from(before_act > 0);
        let next = before as i32 + delta;
        let touch_delta = if touches { delta } else { 0 };
        let next_touch = before_touch as i32 + touch_delta;
        if next <= 0 {
            debug_assert!(
                next_touch <= 0,
                "editable occurrences outlived the k-mer count"
            );
            self.counts[rule].remove(&key);
            self.touching[rule].remove(&key);
            self.excess[rule] -= before_over;
            self.actionable_excess[rule] -= before_act;
            self.offenders[rule] -= before_off;
            return;
        }
        let after = next as u32;
        let after_touch = next_touch.max(0) as u32;
        self.counts[rule].insert(key, after);
        if after_touch == 0 {
            self.touching[rule].remove(&key);
        } else {
            self.touching[rule].insert(key, after_touch);
        }
        let after_over = u64::from(after.saturating_sub(max));
        let after_act = actionable_over(after, after_touch, max);
        let after_off = u64::from(after_act > 0);
        self.excess[rule] = self.excess[rule] - before_over + after_over;
        self.actionable_excess[rule] = self.actionable_excess[rule] - before_act + after_act;
        self.offenders[rule] = self.offenders[rule] - before_off + after_off;
    }

    fn replace(&mut self, nt: &mut [u8], at: usize, bases: [u8; 3]) {
        self.adjust(nt, at, -1);
        nt[at..at + 3].copy_from_slice(&bases);
        self.adjust(nt, at, 1);
    }

    fn adjust(&mut self, nt: &[u8], at: usize, delta: i32) {
        let rule_count = self.rules.len();
        for rule_index in 0..rule_count {
            let k = self.rules[rule_index].k;
            let canonical = self.rules[rule_index].canonical;
            if nt.len() < k {
                continue;
            }
            let lo = at.saturating_sub(k - 1);
            let hi = (at + 3).min(nt.len() - k + 1);
            for pos in lo..hi {
                if let Some(key) = pack(nt, pos, k, canonical) {
                    let touches = spans_overlap(&self.editable, pos, pos + k);
                    self.bump(rule_index, key, delta, touches);
                }
            }
        }
    }
}

/// `count - max` when this key is an actionable k-mer violation, else 0.
///
/// One violation in the checker is one over-limit key, and it is actionable
/// when `touching > 0`. The sum of these contributions is zero under the same
/// condition as "no actionable blocking k-mer violation".
/// Evenly spaced indexes, always including the first and the last.
fn spread_indexes(indexes: &[usize], cap: usize) -> Vec<usize> {
    if cap == 0 || indexes.is_empty() {
        return Vec::new();
    }
    if indexes.len() <= cap {
        return indexes.to_vec();
    }
    if cap == 1 {
        return vec![indexes[0]];
    }
    let mut picked = Vec::with_capacity(cap);
    let last = indexes.len() - 1;
    for slot in 0..cap {
        picked.push(indexes[slot * last / (cap - 1)]);
    }
    picked.sort_unstable();
    picked.dedup();
    picked
}

fn actionable_over(count: u32, touching: u32, max: u32) -> u64 {
    if count > max && touching > 0 {
        u64::from(count - max)
    } else {
        0
    }
}

fn pack(seq: &[u8], pos: usize, k: usize, canonical: bool) -> Option<u64> {
    let mut key = 0u64;
    for offset in 0..k {
        let code = match seq[pos + offset] {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        key = (key << 2) | code;
    }
    if canonical {
        key = key.min(reverse_pack(key, k));
    }
    Some(key)
}

fn reverse_pack(mut packed: u64, k: usize) -> u64 {
    let mut out = 0u64;
    for _ in 0..k {
        out = (out << 2) | (3 - (packed & 3));
        packed >>= 2;
    }
    out
}

struct HomoSpec {
    base: Option<u8>,
    max_run: usize,
    id: String,
}

struct MotifHitRec {
    start: usize,
    end: usize,
    constraint: usize,
}

struct HomoRunRec {
    start: usize,
    end: usize,
    rule: usize,
}

/// Blocking actionable motif and homopolymer violations.
///
/// A codon edit can create or destroy a motif only inside
/// `MotifAutomaton::rescan_radius` of the three bases, and a homopolymer run
/// only inside the run that touches the codon (`homo_window`). The ledger
/// keeps every hit and replaces the ones that intersect that window, then
/// applies the same allowance and actionable test as `check_patterns`.
///
/// A blocking `custom_regex` is not bounded by its source length, so the
/// ledger is not used and `pattern_count` scans the construct.
struct PatternView<'a> {
    checker: &'a ConstraintChecker,
    nt: &'a [u8],
    cds_off: usize,
    codon_count: usize,
    editable: &'a [(usize, usize)],
}

struct PatternLedger {
    incremental: bool,
    radius: usize,
    violations: Vec<(usize, usize, String)>,
    motif_hits: Vec<MotifHitRec>,
    homo_runs: Vec<HomoRunRec>,
    homo_specs: Vec<HomoSpec>,
}

impl PatternLedger {
    fn build(view: PatternView<'_>) -> Self {
        let incremental = !view.checker.set().constraints.iter().any(|constraint| {
            constraint.enabled
                && constraint.class == ConstraintClass::Blocking
                && matches!(constraint.rule, Rule::CustomRegex { .. })
        });
        let automata = view.checker.motif_automata();
        let mut ledger = Self {
            incremental,
            radius: automata[0].rescan_radius().max(automata[1].rescan_radius()),
            violations: Vec::new(),
            motif_hits: Vec::new(),
            homo_runs: Vec::new(),
            homo_specs: homo_specs(view.checker),
        };
        if incremental {
            ledger.motif_hits = motif_hits_in(
                view.checker,
                view.nt,
                view.cds_off,
                view.codon_count,
                0,
                view.nt.len(),
            );
            ledger.homo_runs = homo_runs_in(view.nt, &ledger.homo_specs);
            ledger.materialize(view.checker, view.editable);
        }
        ledger
    }

    fn retile(&mut self, view: PatternView<'_>, at: usize, homo_span: Option<(usize, usize)>) {
        if !self.incremental {
            return;
        }
        let edit_lo = at;
        let edit_hi = at + 3;
        self.motif_hits
            .retain(|hit| hit.end <= edit_lo || hit.start >= edit_hi);
        let scan_lo = at.saturating_sub(self.radius);
        let scan_hi = (at + 3 + self.radius).min(view.nt.len());
        for hit in motif_hits_in(
            view.checker,
            view.nt,
            view.cds_off,
            view.codon_count,
            scan_lo,
            scan_hi,
        ) {
            if hit.start < edit_hi && hit.end > edit_lo {
                self.motif_hits.push(hit);
            }
        }
        match homo_span {
            Some((lo, hi)) if !homo_slice_truncates(view.nt, lo, hi) => {
                self.homo_runs
                    .retain(|run| run.end <= lo || run.start >= hi);
                self.homo_runs
                    .extend(homo_runs_between(view.nt, &self.homo_specs, lo, hi));
            }
            _ => self.homo_runs = homo_runs_in(view.nt, &self.homo_specs),
        }
        self.materialize(view.checker, view.editable);
    }

    fn materialize(&mut self, checker: &ConstraintChecker, editable: &[(usize, usize)]) {
        let mut grouped: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
        for hit in &self.motif_hits {
            grouped
                .entry(hit.constraint)
                .or_default()
                .push((hit.start, hit.end));
        }
        let mut violations = Vec::new();
        for (index, mut occurrences) in grouped {
            let Some(constraint) = checker.set().constraints.get(index) else {
                continue;
            };
            if constraint.class != ConstraintClass::Blocking {
                continue;
            }
            let allowed = match &constraint.rule {
                Rule::ForbiddenMotif {
                    max_occurrences, ..
                }
                | Rule::RestrictionSite {
                    max_occurrences, ..
                } => *max_occurrences as usize,
                Rule::CodonBan { .. } => 0,
                _ => 0,
            };
            if occurrences.len() <= allowed {
                continue;
            }
            occurrences.sort_unstable();
            for (start, end) in occurrences.into_iter().skip(allowed) {
                if spans_overlap(editable, start, end) {
                    violations.push((start, end, constraint.id.clone()));
                }
            }
        }
        for run in &self.homo_runs {
            let Some(spec) = self.homo_specs.get(run.rule) else {
                continue;
            };
            if spans_overlap(editable, run.start, run.end) {
                violations.push((run.start, run.end, spec.id.clone()));
            }
        }
        violations.sort_unstable();
        self.violations = violations;
    }
}

fn motif_hits_in(
    checker: &ConstraintChecker,
    nt: &[u8],
    cds_off: usize,
    codon_count: usize,
    lo: usize,
    hi: usize,
) -> Vec<MotifHitRec> {
    let mut found = Vec::new();
    if lo >= hi || lo >= nt.len() {
        return found;
    }
    let hi = hi.min(nt.len());
    let automata = checker.motif_automata();
    let mut raw = Vec::new();
    automata[0].find_all(&nt[lo..hi], lo, &mut raw);
    for hit in raw.drain(..) {
        let constraint = automata[0].patterns()[hit.pattern].constraint;
        found.push(MotifHitRec {
            start: hit.start,
            end: hit.end,
            constraint,
        });
    }
    let cds_lo = cds_off;
    let cds_hi = cds_off + codon_count * 3;
    let start = lo.max(cds_lo);
    let end = hi.min(cds_hi);
    if start < end {
        automata[1].find_all(&nt[start..end], start, &mut raw);
        for hit in raw {
            let constraint = automata[1].patterns()[hit.pattern].constraint;
            found.push(MotifHitRec {
                start: hit.start,
                end: hit.end,
                constraint,
            });
        }
    }
    found
}

fn homo_specs(checker: &ConstraintChecker) -> Vec<HomoSpec> {
    use crate::spec::HomopolymerBase;
    checker
        .set()
        .constraints
        .iter()
        .filter(|constraint| constraint.enabled && constraint.class == ConstraintClass::Blocking)
        .filter_map(|constraint| match &constraint.rule {
            Rule::Homopolymer { base, max_run } => Some(HomoSpec {
                base: match base {
                    HomopolymerBase::Any => None,
                    HomopolymerBase::A => Some(b'A'),
                    HomopolymerBase::C => Some(b'C'),
                    HomopolymerBase::G => Some(b'G'),
                    HomopolymerBase::T => Some(b'T'),
                },
                max_run: *max_run as usize,
                id: constraint.id.clone(),
            }),
            _ => None,
        })
        .collect()
}

fn homo_runs_in(nt: &[u8], specs: &[HomoSpec]) -> Vec<HomoRunRec> {
    homo_runs_between(nt, specs, 0, nt.len())
}

fn homo_runs_between(nt: &[u8], specs: &[HomoSpec], lo: usize, hi: usize) -> Vec<HomoRunRec> {
    if lo >= hi || lo >= nt.len() {
        return Vec::new();
    }
    let hi = hi.min(nt.len());
    let mut runs = Vec::new();
    for (rule, spec) in specs.iter().enumerate() {
        for run in crate::window::find_runs(&nt[lo..hi], spec.max_run + 1, spec.base) {
            runs.push(HomoRunRec {
                start: run.start + lo,
                end: run.end + lo,
                rule,
            });
        }
    }
    runs
}

/// `true` when a maximal run continues past `[lo, hi)`, so a slice scan would
/// report that run truncated.
fn homo_slice_truncates(nt: &[u8], lo: usize, hi: usize) -> bool {
    if lo > 0 && nt[lo - 1] == nt[lo] && matches!(nt[lo], b'A' | b'C' | b'G' | b'T') {
        return true;
    }
    hi < nt.len()
        && hi > 0
        && nt[hi - 1] == nt[hi]
        && matches!(nt[hi - 1], b'A' | b'C' | b'G' | b'T')
}

struct Solver<'a> {
    seq: CodonSequence,
    checker: &'a ConstraintChecker,
    weights: &'a [f64; 64],
    nt: Vec<u8>,
    cds_off: usize,
    prefix: Vec<u32>,
    /// GC-count edits not yet written into `prefix`. Each entry adds `delta`
    /// to every prefix index at or after `from`. Trials push and cancel these
    /// instead of rewriting the tail, which is the whole construct.
    prefix_adds: Vec<(usize, i32)>,
    gc_penalty: f64,
    bands: Vec<GcBand>,
    kmers: KmerBook,
    editable: Vec<(usize, usize)>,
    homo_rules: Vec<(Option<u8>, usize)>,
    homo_penalty: u64,
    patterns: PatternLedger,
    /// Every byte of `nt` is `A/C/G/T`. Codons are, and flanks are normalized
    /// to that alphabet before a solver exists. The GC proof below uses it.
    acgt: bool,
    edits: u32,
    tie: TieBreak,
    limits: SpendLimit,
}

impl<'a> Solver<'a> {
    fn new(
        seq: CodonSequence,
        checker: &'a ConstraintChecker,
        weights: &'a [f64; 64],
        tie: TieBreak,
    ) -> Result<Self> {
        let cds_off = seq.five_prime_flank().len();
        let mut nt = seq.five_prime_flank().to_vec();
        for codon in seq.codons() {
            nt.extend_from_slice(&codon.to_bytes());
        }
        nt.extend_from_slice(seq.three_prime_flank());
        let prefix = gc_prefix(&nt);
        let editable = editable_spans_of(&seq);
        let mut bands = gc_bands(checker);
        for band in &mut bands {
            band.offenders = count_gc_offenders(&nt, &prefix, band, &editable);
        }
        let gc_penalty = gc_penalty_of(&nt, &prefix, &bands, &editable);
        let kmers = KmerBook::build(&nt, kmer_rules(checker), &editable);
        let homo_rules = homo_rules(checker);
        let homo_penalty = homo_total(&nt, 0, None, &homo_rules, &editable);
        let acgt = nt
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'));
        let patterns = PatternLedger::build(PatternView {
            checker,
            nt: &nt,
            cds_off,
            codon_count: seq.codon_count(),
            editable: &editable,
        });
        Ok(Self {
            seq,
            checker,
            weights,
            nt,
            cds_off,
            prefix,
            prefix_adds: Vec::new(),
            gc_penalty,
            bands,
            kmers,
            editable,
            homo_rules,
            homo_penalty,
            patterns,
            acgt,
            edits: 0,
            tie,
            limits: SpendLimit::default(),
        })
    }

    fn with_limits(mut self, limits: SpendLimit) -> Self {
        self.limits = limits;
        self
    }

    fn enforce(&self) -> Result<()> {
        self.limits.enforce()
    }

    fn prefer(&self, candidate: Codon, incumbent: Codon) -> bool {
        match self.tie {
            TieBreak::Weight => prefer_weight(self.weights, candidate, incumbent),
            TieBreak::LowIndex => candidate.index() < incumbent.index(),
        }
    }

    fn run(&mut self) -> Result<u32> {
        let mut cycles = 0u32;
        let mut escapes = 0u32;
        self.enforce()?;
        if self.is_feasible() {
            self.polish()?;
            self.ascend_cai()?;
            self.climb_conflicts()?;
            self.finish_prefix();
            return Ok(0);
        }
        for _ in 0..MAX_CYCLES {
            self.enforce()?;
            cycles += 1;
            let before = self.score();
            self.gc_sweeps()?;
            for _ in 0..PATTERN_ROUNDS {
                if !self.fix_patterns()? {
                    break;
                }
            }
            for _ in 0..KMER_ROUNDS {
                if self.kmers.excess() == 0 {
                    break;
                }
                if self.fix_kmers()? {
                    continue;
                }
                // A single codon can sit in a repeat that no one-codon edit
                // breaks without creating another. A nearby pair can.
                if !self.fix_kmer_pairs()? {
                    break;
                }
            }
            if self.is_feasible() {
                self.polish()?;
                self.ascend_cai()?;
                self.climb_conflicts()?;
                break;
            }
            if self.score() >= before {
                // Single-codon GC sweeps stop while a window is still outside
                // the band, and k-mer repair will not edit until that penalty
                // is gone. A pair inside the window can move the count when no
                // one codon can. This runs only after a cycle made no progress,
                // so a seed that already became feasible never reaches it.
                if escapes < 4 && self.ease_gc_pairs()? {
                    escapes += 1;
                    continue;
                }
                break;
            }
        }
        if !self.is_feasible() {
            cycles += self.joint_escape()?;
        }
        self.finish_prefix();
        Ok(cycles)
    }

    /// Pair edits that can move a GC window and a repeated word together.
    ///
    /// `DEP-SRCH-024`. Reached only when the cycles above, including the
    /// middle-of-window GC pairs, left the seed infeasible. Sites are spread
    /// across the whole stuck window and across the hot k-mer codons, so the
    /// edges are eligible. A pair is kept only when the packed progress score
    /// falls and pattern hits do not rise. A feasible seed never enters.
    fn joint_escape(&mut self) -> Result<u32> {
        let mut extra = 0u32;
        for _ in 0..JOINT_ATTEMPTS {
            self.enforce()?;
            if self.is_feasible() || !self.ease_joint()? {
                break;
            }
            extra += 1;
            let before = self.score();
            self.gc_sweeps()?;
            for _ in 0..PATTERN_ROUNDS {
                if !self.fix_patterns()? {
                    break;
                }
            }
            for _ in 0..KMER_ROUNDS {
                if self.kmers.excess() == 0 {
                    break;
                }
                if self.fix_kmers()? {
                    continue;
                }
                if !self.fix_kmer_pairs()? {
                    break;
                }
            }
            if self.is_feasible() {
                self.polish()?;
                self.ascend_cai()?;
                self.climb_conflicts()?;
                break;
            }
            if self.score() >= before {
                break;
            }
        }
        Ok(extra)
    }

    fn ease_joint(&mut self) -> Result<bool> {
        if self.is_feasible() {
            return Ok(false);
        }
        let before = self.score();
        let patterns_before = self.pattern_count();
        let codons = self.joint_sites();
        if codons.len() < 2 {
            return Ok(false);
        }
        let mut best: Option<(u128, usize, usize, Codon, Codon)> = None;
        for left in 0..codons.len() {
            for right in (left + 1)..codons.len() {
                let i = codons[left];
                let j = codons[right];
                let current_i = self.seq.codons()[i];
                let current_j = self.seq.codons()[j];
                let syn_i = self.seq.synonyms_at(i).to_vec();
                let syn_j = self.seq.synonyms_at(j).to_vec();
                for &synonym_i in &syn_i {
                    for &synonym_j in &syn_j {
                        if synonym_i == current_i && synonym_j == current_j {
                            continue;
                        }
                        self.set_codon(i, synonym_i)?;
                        self.set_codon(j, synonym_j)?;
                        let trial = self.score();
                        let patterns = self.pattern_count();
                        self.set_codon(j, current_j)?;
                        self.set_codon(i, current_i)?;
                        if trial >= before || patterns > patterns_before {
                            continue;
                        }
                        let replace = best
                            .map(|(so_far, bi, bj, incumbent_i, incumbent_j)| {
                                trial < so_far
                                    || (trial == so_far
                                        && (i < bi
                                            || (i == bi
                                                && (j < bj
                                                    || (j == bj
                                                        && (self
                                                            .prefer(synonym_i, incumbent_i)
                                                            || (synonym_i == incumbent_i
                                                                && self.prefer(
                                                                    synonym_j,
                                                                    incumbent_j,
                                                                ))))))))
                            })
                            .unwrap_or(true);
                        if replace {
                            best = Some((trial, i, j, synonym_i, synonym_j));
                        }
                    }
                }
            }
        }
        let Some((_, i, j, synonym_i, synonym_j)) = best else {
            return Ok(false);
        };
        self.commit(i, synonym_i)?;
        self.commit(j, synonym_j)?;
        self.recompute_gc();
        Ok(self.score() < before)
    }

    fn joint_sites(&self) -> Vec<usize> {
        let mut sites = Vec::new();
        if let Some((start, end)) = self.first_outside_gc_window() {
            let codons = self.codons_covering(start, end);
            sites.extend(spread_indexes(&codons, JOINT_GC_CODONS));
        }
        let hot: Vec<usize> = self
            .hot_codons()
            .into_iter()
            .filter(|index| !self.seq.is_locked(*index))
            .collect();
        sites.extend(spread_indexes(&hot, JOINT_KMER_CODONS));
        sites.sort_unstable();
        sites.dedup();
        spread_indexes(&sites, JOINT_SITES)
    }

    fn score(&self) -> u128 {
        // Pattern violations dominate, then k-mer excess, then GC penalty
        // scaled so a smaller penalty sorts lower. Used only as a progress
        // stop, not as the feasibility verdict.
        let patterns = self.pattern_count() as u128;
        let kmers = self.kmers.excess() as u128;
        let gc = (self.gc_penalty * 1_000_000.0).max(0.0) as u128;
        (patterns << 96) | (kmers << 64) | gc
    }

    fn gc_offender_total(&self) -> u32 {
        self.bands.iter().map(|band| band.offenders).sum()
    }

    /// Codons that overlap `[start, end)`, skipping a locked start or stop.
    fn codons_covering(&self, start: usize, end: usize) -> Vec<usize> {
        if self.seq.codon_count() == 0 || end <= self.cds_off {
            return Vec::new();
        }
        let last = self.seq.codon_count() - 1;
        let first = start.saturating_sub(self.cds_off) / 3;
        let last_nt = end.saturating_sub(1).saturating_sub(self.cds_off);
        let end_codon = (last_nt / 3).min(last);
        (first..=end_codon)
            .filter(|index| !self.seq.is_locked(*index))
            .collect()
    }

    fn first_outside_gc_window(&self) -> Option<(usize, usize)> {
        let n = self.nt.len();
        for band in &self.bands {
            if band.offenders == 0 || band.window == 0 || n < band.window {
                continue;
            }
            for start in 0..=n - band.window {
                let end = start + band.window;
                if self.overlaps_editable(start, end)
                    && window_outside(self.prefix_at(end) - self.prefix_at(start), band)
                {
                    return Some((start, end));
                }
            }
        }
        None
    }

    /// One pair of synonymous edits inside a GC window that is still outside
    /// its band after the single-codon sweeps.
    ///
    /// The pair may raise k-mer excess. The next cycle repairs repeats once
    /// the window count has moved. Pattern hits are not allowed to rise.
    fn ease_gc_pairs(&mut self) -> Result<bool> {
        let before = self.gc_offender_total();
        if before == 0 {
            return Ok(false);
        }
        let Some((start, end)) = self.first_outside_gc_window() else {
            return Ok(false);
        };
        let mut codons = self.codons_covering(start, end);
        if codons.len() > 10 {
            let lo = codons.len() / 2 - 5;
            codons = codons[lo..lo + 10].to_vec();
        }
        if codons.len() < 2 {
            return Ok(false);
        }
        let patterns_before = self.pattern_count();
        let mut best: Option<(u32, usize, usize, Codon, Codon)> = None;
        for left in 0..codons.len() {
            for right in (left + 1)..codons.len() {
                let i = codons[left];
                let j = codons[right];
                let current_i = self.seq.codons()[i];
                let current_j = self.seq.codons()[j];
                let syn_i = self.seq.synonyms_at(i).to_vec();
                let syn_j = self.seq.synonyms_at(j).to_vec();
                for &synonym_i in &syn_i {
                    for &synonym_j in &syn_j {
                        if synonym_i == current_i && synonym_j == current_j {
                            continue;
                        }
                        self.set_codon(i, synonym_i)?;
                        self.set_codon(j, synonym_j)?;
                        let trial = self.gc_offender_total();
                        let patterns = self.pattern_count();
                        self.set_codon(j, current_j)?;
                        self.set_codon(i, current_i)?;
                        if trial >= before || patterns > patterns_before {
                            continue;
                        }
                        let replace = best
                            .map(|(so_far, bi, bj, incumbent_i, incumbent_j)| {
                                trial < so_far
                                    || (trial == so_far
                                        && (i < bi
                                            || (i == bi
                                                && (j < bj
                                                    || (j == bj
                                                        && (self
                                                            .prefer(synonym_i, incumbent_i)
                                                            || (synonym_i == incumbent_i
                                                                && self.prefer(
                                                                    synonym_j,
                                                                    incumbent_j,
                                                                ))))))))
                            })
                            .unwrap_or(true);
                        if replace {
                            best = Some((trial, i, j, synonym_i, synonym_j));
                        }
                    }
                }
            }
        }
        let Some((_, i, j, synonym_i, synonym_j)) = best else {
            return Ok(false);
        };
        self.commit(i, synonym_i)?;
        self.commit(j, synonym_j)?;
        self.recompute_gc();
        Ok(self.gc_offender_total() < before)
    }

    fn gc_sweeps(&mut self) -> Result<()> {
        let started = seed_profile_on().then(Instant::now);
        let result = self.gc_sweeps_inner();
        if let Some(started) = started {
            SEED_PROFILE.gc_sweeps.record(started);
        }
        result
    }

    fn gc_sweeps_inner(&mut self) -> Result<()> {
        for _ in 0..GC_SWEEPS {
            if self.gc_penalty <= 1e-12 {
                break;
            }
            let mut committed = false;
            for index in 0..self.seq.codon_count() {
                if self.seq.is_locked(index) {
                    continue;
                }
                let current = self.seq.codons()[index];
                let at = self.cds_off + index * 3;
                let mut best = current;
                let mut best_penalty = self.gc_penalty;
                let mut best_homo = self.homo_touching(at, None);
                for &synonym in self.seq.synonyms_at(index) {
                    if synonym == current {
                        continue;
                    }
                    let bases = synonym.to_bytes();
                    let penalty = self.gc_penalty + self.gc_delta(at, bases);
                    let homo = self.homo_touching(at, Some(bases));
                    if penalty < best_penalty - 1e-12
                        || ((penalty - best_penalty).abs() <= 1e-12
                            && (homo < best_homo
                                || (homo == best_homo && self.prefer(synonym, best))))
                    {
                        best = synonym;
                        best_penalty = penalty;
                        best_homo = homo;
                    }
                }
                if best != current && best_penalty < self.gc_penalty - 1e-12 {
                    self.commit(index, best)?;
                    committed = true;
                }
            }
            self.recompute_gc();
            if !committed {
                break;
            }
        }
        Ok(())
    }

    fn fix_patterns(&mut self) -> Result<bool> {
        let started = seed_profile_on().then(Instant::now);
        let result = self.fix_patterns_inner();
        if let Some(started) = started {
            SEED_PROFILE.fix_patterns.record(started);
        }
        result
    }

    fn fix_patterns_inner(&mut self) -> Result<bool> {
        let violations = self.pattern_violations();
        if violations.is_empty() {
            return Ok(false);
        }
        // The span list is a snapshot. The ceiling is live, so a later span
        // is accepted only when it reduces the count the earlier edit left.
        let mut ceiling = violations[0].2;
        let mut progress = false;
        for violation in violations {
            let mut best: Option<(usize, Codon, usize)> = None;
            let lo = violation.0.saturating_sub(3);
            let hi = violation.1 + 3;
            for index in 0..self.seq.codon_count() {
                if self.seq.is_locked(index) {
                    continue;
                }
                let at = self.cds_off + index * 3;
                if at + 3 <= lo || at >= hi {
                    continue;
                }
                let current = self.seq.codons()[index];
                for synonym in self.seq.synonyms_at(index).to_vec() {
                    if synonym == current {
                        continue;
                    }
                    self.set_codon(index, synonym)?;
                    let count = self.pattern_count();
                    let improved = count < ceiling;
                    self.set_codon(index, current)?;
                    let better = best.map(|(c, _, _)| count < c).unwrap_or(true);
                    let tied = best.is_some_and(|(c, incumbent, _)| {
                        count == c && self.prefer(synonym, incumbent)
                    });
                    if improved && (better || tied) {
                        best = Some((count, synonym, index));
                    }
                }
            }
            if let Some((count, synonym, index)) = best {
                self.commit(index, synonym)?;
                ceiling = count;
                progress = true;
            }
        }
        Ok(progress)
    }

    fn fix_kmers(&mut self) -> Result<bool> {
        let started = seed_profile_on().then(Instant::now);
        let result = self.fix_kmers_inner();
        if let Some(started) = started {
            SEED_PROFILE.fix_kmers.record(started);
        }
        result
    }

    fn fix_kmers_inner(&mut self) -> Result<bool> {
        let before = self.kmers.excess();
        if before == 0 {
            return Ok(false);
        }
        let hot = self.hot_codons();
        let mut progress = false;
        let mut excess = before;
        for (step, index) in hot.into_iter().enumerate() {
            if step % 32 == 0 {
                self.enforce()?;
            }
            if self.seq.is_locked(index) {
                continue;
            }
            let current = self.seq.codons()[index];
            let patterns_before = self.pattern_count();
            let mut best: Option<Codon> = None;
            let mut best_excess = excess;
            for synonym in self.seq.synonyms_at(index).to_vec() {
                if synonym == current {
                    continue;
                }
                self.set_codon(index, synonym)?;
                let trial = self.kmers.excess();
                let patterns = self.pattern_count();
                let gc_ok = self.gc_penalty <= 1e-9;
                self.set_codon(index, current)?;
                if trial < best_excess && patterns <= patterns_before && gc_ok {
                    best = Some(synonym);
                    best_excess = trial;
                } else if trial == best_excess
                    && trial < excess
                    && patterns <= patterns_before
                    && gc_ok
                    && best.is_some_and(|incumbent| self.prefer(synonym, incumbent))
                {
                    best = Some(synonym);
                }
            }
            if let Some(synonym) = best {
                self.commit(index, synonym)?;
                excess = best_excess;
                progress = true;
                if excess == 0 {
                    break;
                }
            }
        }
        Ok(progress)
    }

    /// Try two nearby codon edits when no single edit reduces k-mer excess.
    fn fix_kmer_pairs(&mut self) -> Result<bool> {
        let started = seed_profile_on().then(Instant::now);
        let result = self.fix_kmer_pairs_inner();
        if let Some(started) = started {
            SEED_PROFILE.fix_kmer_pairs.record(started);
        }
        result
    }

    fn fix_kmer_pairs_inner(&mut self) -> Result<bool> {
        let excess = self.kmers.excess();
        if excess == 0 || self.gc_penalty > 1e-9 {
            return Ok(false);
        }
        let mut seen = std::collections::BTreeSet::new();
        let hot: Vec<usize> = self
            .hot_codons()
            .into_iter()
            .filter(|index| seen.insert(*index) && !self.seq.is_locked(*index))
            .take(16)
            .collect();
        let patterns_before = self.pattern_count();
        let mut pairs = Vec::new();
        for left in 0..hot.len() {
            for right in (left + 1)..hot.len() {
                if hot[right] <= hot[left] + 4 {
                    pairs.push((hot[left], hot[right]));
                }
            }
            let index = hot[left];
            // The compensating edit often sits on the codon just outside the
            // repeated word, which `hot_codons` does not list.
            for neighbor in [index.saturating_sub(1), index + 1] {
                if neighbor == index || neighbor >= self.seq.codon_count() {
                    continue;
                }
                if self.seq.is_locked(neighbor) {
                    continue;
                }
                let pair = if neighbor < index {
                    (neighbor, index)
                } else {
                    (index, neighbor)
                };
                pairs.push(pair);
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        for (i, j) in pairs {
            let current_i = self.seq.codons()[i];
            let current_j = self.seq.codons()[j];
            let syn_i = self.seq.synonyms_at(i).to_vec();
            let syn_j = self.seq.synonyms_at(j).to_vec();
            let mut best: Option<(u64, Codon, Codon)> = None;
            for synonym_i in syn_i {
                for synonym_j in syn_j.clone() {
                    if synonym_i == current_i && synonym_j == current_j {
                        continue;
                    }
                    self.set_codon(i, synonym_i)?;
                    self.set_codon(j, synonym_j)?;
                    let trial = self.kmers.excess();
                    let patterns = self.pattern_count();
                    let gc_ok = self.gc_penalty <= 1e-9;
                    self.set_codon(j, current_j)?;
                    self.set_codon(i, current_i)?;
                    if trial < excess && patterns <= patterns_before && gc_ok {
                        let replace = best
                            .map(|(so_far, incumbent_i, incumbent_j)| {
                                trial < so_far
                                    || (trial == so_far
                                        && (self.prefer(synonym_i, incumbent_i)
                                            || (synonym_i == incumbent_i
                                                && self.prefer(synonym_j, incumbent_j))))
                            })
                            .unwrap_or(true);
                        if replace {
                            best = Some((trial, synonym_i, synonym_j));
                        }
                    }
                }
            }
            if let Some((trial, synonym_i, synonym_j)) = best {
                self.commit(i, synonym_i)?;
                self.commit(j, synonym_j)?;
                if trial == 0 {
                    return Ok(true);
                }
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn polish(&mut self) -> Result<()> {
        let started = seed_profile_on().then(Instant::now);
        let result = self.polish_inner();
        if let Some(started) = started {
            SEED_PROFILE.polish.record(started);
        }
        result
    }

    fn polish_inner(&mut self) -> Result<()> {
        if !self.is_feasible() {
            return Ok(());
        }
        let saved = self.seq.clone();
        let edits_before = self.edits;
        for index in 0..saved.codon_count() {
            if self.seq.is_locked(index) {
                continue;
            }
            let current = self.seq.codons()[index];
            let mut best = current;
            for synonym in self.seq.synonyms_at(index).to_vec() {
                if synonym == current || !prefer_weight(self.weights, synonym, best) {
                    continue;
                }
                if self.weights[synonym.index_usize()] == 0.0
                    && self.weights[best.index_usize()] > 0.0
                {
                    continue;
                }
                self.set_codon(index, synonym)?;
                let ok = self.gc_penalty <= 1e-9
                    && self.kmers.excess() == 0
                    && self.pattern_count() == 0
                    && self.is_feasible();
                self.set_codon(index, current)?;
                if ok {
                    best = synonym;
                }
            }
            if best != current {
                self.commit(index, best)?;
            }
        }
        if !self.is_feasible() {
            let limits = self.limits;
            *self = Self::new(saved, self.checker, self.weights, self.tie)?.with_limits(limits);
            self.edits = edits_before;
        }
        Ok(())
    }

    /// Raise CAI inside the feasible set by coordinated pair edits.
    ///
    /// Single-codon polish already ran. A higher-weight codon that recreates a
    /// repeated word is rejected there even when a second codon, a few positions
    /// away, would clear that word. Two passes, each trying the positions with
    /// the largest remaining weight gap against a few nearby partners.
    fn ascend_cai(&mut self) -> Result<()> {
        let started = seed_profile_on().then(Instant::now);
        let result = self.ascend_cai_inner();
        if let Some(started) = started {
            SEED_PROFILE.ascend.record(started);
        }
        result
    }

    fn ascend_cai_inner(&mut self) -> Result<()> {
        if !self.is_feasible() {
            return Ok(());
        }
        for _ in 0..ASCENT_PASSES {
            if !self.lift_pairs()? {
                break;
            }
            self.polish()?;
        }
        Ok(())
    }

    fn ascent_positions(&self) -> Vec<usize> {
        let mut ranked = Vec::new();
        for index in 0..self.seq.codon_count() {
            if self.seq.is_locked(index) {
                continue;
            }
            let current = self.weights[self.seq.codons()[index].index_usize()];
            let mut best = current;
            for synonym in self.seq.synonyms_at(index) {
                let weight = self.weights[synonym.index_usize()];
                if weight > best {
                    best = weight;
                }
            }
            let gap = best - current;
            if gap > 0.0 {
                ranked.push((gap, index));
            }
        }
        ranked.sort_by(|left, right| match right.0.total_cmp(&left.0) {
            Ordering::Equal => left.1.cmp(&right.1),
            other => other,
        });
        ranked
            .into_iter()
            .take(ASCENT_POSITIONS)
            .map(|(_, index)| index)
            .collect()
    }

    fn partners_for(&self, hot: &[usize], index: usize) -> Vec<usize> {
        let mut pool = Vec::new();
        for &partner in hot {
            if partner != index {
                pool.push(partner);
            }
        }
        let last = self.seq.codon_count().saturating_sub(1);
        let lo = index.saturating_sub(ASCENT_NEIGHBOUR);
        let hi = (index + ASCENT_NEIGHBOUR).min(last);
        for partner in lo..=hi {
            if partner != index && !self.seq.is_locked(partner) {
                pool.push(partner);
            }
        }
        pool.sort_unstable();
        pool.dedup();
        let mut ranked: Vec<(usize, usize)> = pool
            .into_iter()
            .map(|partner| (index.abs_diff(partner), partner))
            .collect();
        ranked.sort_unstable();
        ranked
            .into_iter()
            .take(ASCENT_PARTNERS)
            .map(|(_, partner)| partner)
            .collect()
    }

    fn lift_pairs(&mut self) -> Result<bool> {
        if !self.is_feasible() {
            return Ok(false);
        }
        let hot = self.ascent_positions();
        if hot.is_empty() {
            return Ok(false);
        }
        let mut cai_before = cai_of(&self.seq, self.weights);
        let mut improved = false;
        for &index in &hot {
            if self.seq.is_locked(index) {
                continue;
            }
            let partners = self.partners_for(&hot, index);
            if partners.is_empty() {
                continue;
            }
            let current_i = self.seq.codons()[index];
            let syn_i = self.seq.synonyms_at(index).to_vec();
            let mut best: Option<(f64, usize, Codon, Codon)> = None;
            for &partner in &partners {
                let current_j = self.seq.codons()[partner];
                let syn_j = self.seq.synonyms_at(partner).to_vec();
                for &codon_i in &syn_i {
                    for &codon_j in &syn_j {
                        if codon_i == current_i && codon_j == current_j {
                            continue;
                        }
                        self.set_codon(index, codon_i)?;
                        self.set_codon(partner, codon_j)?;
                        let cai = cai_of(&self.seq, self.weights);
                        let ok = cai > cai_before
                            && self.gc_penalty <= 1e-9
                            && self.kmers.excess() == 0
                            && self.is_feasible();
                        self.set_codon(partner, current_j)?;
                        self.set_codon(index, current_i)?;
                        if !ok {
                            continue;
                        }
                        let replace = match &best {
                            None => true,
                            Some((so_far, so_partner, so_i, so_j)) => match cai.total_cmp(so_far) {
                                Ordering::Greater => true,
                                Ordering::Less => false,
                                Ordering::Equal => {
                                    (partner, codon_i.index(), codon_j.index())
                                        < (*so_partner, so_i.index(), so_j.index())
                                }
                            },
                        };
                        if replace {
                            best = Some((cai, partner, codon_i, codon_j));
                        }
                    }
                }
            }
            if let Some((_, partner, codon_i, codon_j)) = best {
                self.commit(index, codon_i)?;
                self.commit(partner, codon_j)?;
                cai_before = cai_of(&self.seq, self.weights);
                improved = true;
            }
        }
        Ok(improved)
    }

    /// Raise CAI by repairing the constraint a higher-weight codon breaks.
    ///
    /// `DEP-SRCH-018` only pairs a hot codon with a neighbour a few positions
    /// away. A unique 12-mer is broken by editing the other copy, which can sit
    /// anywhere in the CDS, and a GC window of 50 nt can be repaired by a codon
    /// that is inside the window and outside that neighbourhood. Each pass tries
    /// every counted position that is not already on its heaviest synonym. The
    /// edit is kept only when the geometric mean rises and the sequence stays
    /// feasible. `DEP-SRCH-021`.
    fn climb_conflicts(&mut self) -> Result<()> {
        let started = seed_profile_on().then(Instant::now);
        let result = self.climb_conflicts_inner();
        if let Some(started) = started {
            SEED_PROFILE.climb.record(started);
        }
        result
    }

    fn climb_conflicts_inner(&mut self) -> Result<()> {
        if !self.is_feasible() {
            return Ok(());
        }
        for _ in 0..CLIMB_PASSES {
            self.enforce()?;
            let order = self.climb_positions();
            if order.is_empty() {
                break;
            }
            let mut improved = false;
            for index in order {
                if self.lift_at(index)? {
                    improved = true;
                }
            }
            if !improved {
                break;
            }
        }
        Ok(())
    }

    fn climb_positions(&self) -> Vec<usize> {
        let mut ranked = Vec::new();
        for index in 0..self.seq.codon_count() {
            if self.seq.is_locked(index) || !self.position_counted(index) {
                continue;
            }
            let current = self.seq.codons()[index];
            let current_weight = self.weights[current.index_usize()];
            let mut best = current_weight;
            for synonym in self.seq.synonyms_at(index) {
                let weight = self.weights[synonym.index_usize()];
                if weight > best {
                    best = weight;
                }
            }
            let gap = best - current_weight;
            if gap > 0.0 {
                ranked.push((gap, index));
            }
        }
        ranked.sort_by(|left, right| match right.0.total_cmp(&left.0) {
            Ordering::Equal => left.1.cmp(&right.1),
            other => other,
        });
        ranked.into_iter().map(|(_, index)| index).collect()
    }

    fn position_counted(&self, index: usize) -> bool {
        let codon = self.seq.codons()[index];
        !self.seq.table().is_stop(codon) && self.seq.synonyms_at(index).len() > 1
    }

    /// Log-weight change of replacing `from` with `to` at `index`.
    ///
    /// `None` when the replacement introduces a counted weight of 0, which
    /// makes the geometric mean 0. An uncounted position contributes 0.
    /// Leaving a zero weight is an infinite rise.
    fn log_shift(&self, index: usize, from: Codon, to: Codon) -> Option<f64> {
        if !self.position_counted(index) {
            return Some(0.0);
        }
        let old_weight = self.weights[from.index_usize()];
        let new_weight = self.weights[to.index_usize()];
        if new_weight == 0.0 && old_weight > 0.0 {
            return None;
        }
        if old_weight == 0.0 && new_weight > 0.0 {
            return Some(f64::INFINITY);
        }
        if old_weight == 0.0 || new_weight == 0.0 {
            return Some(0.0);
        }
        Some(new_weight.ln() - old_weight.ln())
    }

    fn lift_at(&mut self, index: usize) -> Result<bool> {
        if self.seq.is_locked(index) || !self.position_counted(index) {
            return Ok(false);
        }
        let current = self.seq.codons()[index];
        let cai_before = cai_of(&self.seq, self.weights);
        let mut synonyms: Vec<Codon> = self.seq.synonyms_at(index).to_vec();
        synonyms.retain(|codon| {
            *codon != current
                && self
                    .log_shift(index, current, *codon)
                    .is_some_and(|shift| shift > 0.0)
        });
        synonyms.sort_by(|left, right| {
            let left_shift = self.log_shift(index, current, *left).unwrap_or(0.0);
            let right_shift = self.log_shift(index, current, *right).unwrap_or(0.0);
            match right_shift.total_cmp(&left_shift) {
                Ordering::Equal => left.index().cmp(&right.index()),
                other => other,
            }
        });
        let mut best: Option<ClimbChoice> = None;
        for synonym in synonyms {
            let primary_shift = self.log_shift(index, current, synonym).unwrap_or(0.0);
            self.set_codon(index, synonym)?;
            if self.is_feasible() {
                self.note_climb(synonym, Vec::new(), cai_before, &mut best);
            } else {
                let sites = self.violation_sites(index);
                let singles = self.single_repairs(primary_shift, &sites);
                let mut repaired = false;
                for (site, codon) in singles {
                    let partner = self.seq.codons()[site];
                    self.set_codon(site, codon)?;
                    if self.is_feasible() {
                        repaired = true;
                        self.note_climb(synonym, vec![(site, codon)], cai_before, &mut best);
                    }
                    self.set_codon(site, partner)?;
                }
                if !repaired {
                    for partners in self.pair_repairs(primary_shift, &sites) {
                        let saved: Vec<(usize, Codon)> = partners
                            .iter()
                            .map(|(site, _)| (*site, self.seq.codons()[*site]))
                            .collect();
                        for &(site, codon) in &partners {
                            self.set_codon(site, codon)?;
                        }
                        if self.is_feasible() {
                            self.note_climb(synonym, partners.clone(), cai_before, &mut best);
                        }
                        for (site, codon) in saved.into_iter().rev() {
                            self.set_codon(site, codon)?;
                        }
                    }
                }
            }
            self.set_codon(index, current)?;
        }
        let Some(choice) = best else {
            return Ok(false);
        };
        self.commit(index, choice.primary)?;
        for (site, codon) in choice.partners {
            self.commit(site, codon)?;
        }
        Ok(true)
    }

    fn note_climb(
        &mut self,
        primary: Codon,
        partners: Vec<(usize, Codon)>,
        cai_before: f64,
        best: &mut Option<ClimbChoice>,
    ) {
        let cai = cai_of(&self.seq, self.weights);
        if cai <= cai_before {
            return;
        }
        let choice = ClimbChoice {
            cai,
            primary,
            partners,
        };
        let replace = match best.as_ref() {
            None => true,
            Some(incumbent) => choice.beats(incumbent),
        };
        if replace {
            *best = Some(choice);
        }
    }

    fn single_repairs(&self, primary_shift: f64, sites: &[usize]) -> Vec<(usize, Codon)> {
        let mut candidates = Vec::new();
        for &site in sites {
            let current = self.seq.codons()[site];
            for &codon in self.seq.synonyms_at(site) {
                if codon == current {
                    continue;
                }
                let Some(partner_shift) = self.log_shift(site, current, codon) else {
                    continue;
                };
                let total = primary_shift + partner_shift;
                if total <= 0.0 {
                    continue;
                }
                candidates.push((total, site, codon));
            }
        }
        candidates.sort_by(|left, right| match right.0.total_cmp(&left.0) {
            Ordering::Equal => (left.1, left.2.index()).cmp(&(right.1, right.2.index())),
            other => other,
        });
        candidates.truncate(CLIMB_PARTNER_TRIALS);
        candidates
            .into_iter()
            .map(|(_, site, codon)| (site, codon))
            .collect()
    }

    fn pair_repairs(&self, primary_shift: f64, sites: &[usize]) -> Vec<Vec<(usize, Codon)>> {
        let sites = &sites[..sites.len().min(CLIMB_PAIR_SITES)];
        let mut candidates = Vec::new();
        for left in 0..sites.len() {
            for right in (left + 1)..sites.len() {
                let site_a = sites[left];
                let site_b = sites[right];
                let current_a = self.seq.codons()[site_a];
                let current_b = self.seq.codons()[site_b];
                for &codon_a in self.seq.synonyms_at(site_a) {
                    if codon_a == current_a {
                        continue;
                    }
                    let Some(shift_a) = self.log_shift(site_a, current_a, codon_a) else {
                        continue;
                    };
                    for &codon_b in self.seq.synonyms_at(site_b) {
                        if codon_b == current_b {
                            continue;
                        }
                        let Some(shift_b) = self.log_shift(site_b, current_b, codon_b) else {
                            continue;
                        };
                        let total = primary_shift + shift_a + shift_b;
                        if total <= 0.0 {
                            continue;
                        }
                        candidates.push((total, vec![(site_a, codon_a), (site_b, codon_b)]));
                    }
                }
            }
        }
        candidates.sort_by(|left, right| match right.0.total_cmp(&left.0) {
            Ordering::Equal => {
                let key = |partners: &Vec<(usize, Codon)>| -> Vec<(usize, u8)> {
                    partners
                        .iter()
                        .map(|(site, codon)| (*site, codon.index()))
                        .collect()
                };
                key(&left.1).cmp(&key(&right.1))
            }
            other => other,
        });
        candidates.truncate(CLIMB_PARTNER_TRIALS);
        candidates
            .into_iter()
            .map(|(_, partners)| partners)
            .collect()
    }

    /// Codons that overlap a violation created by the edit at `primary`.
    ///
    /// The sequence was feasible before that edit, so every current violation
    /// was introduced by it. Sites are nearest-first and capped. A 50 nt GC
    /// window fits under the cap.
    fn violation_sites(&self, primary: usize) -> Vec<usize> {
        let mut sites = BTreeSet::new();
        let at = self.cds_off + primary * 3;
        for (start, end, _) in self.pattern_spans() {
            self.collect_codons(start, end, primary, &mut sites);
        }
        self.collect_gc_sites(at, primary, &mut sites);
        self.collect_kmer_sites(at, primary, &mut sites);
        let mut ranked: Vec<(usize, usize)> = sites
            .into_iter()
            .map(|index| (index.abs_diff(primary), index))
            .collect();
        ranked.sort_unstable();
        ranked
            .into_iter()
            .take(CLIMB_SITES)
            .map(|(_, index)| index)
            .collect()
    }

    fn collect_gc_sites(&self, at: usize, primary: usize, sites: &mut BTreeSet<usize>) {
        let n = self.nt.len();
        for band in &self.bands {
            if band.window == 0 || n < band.window || band.offenders == 0 {
                continue;
            }
            let lo = at.saturating_sub(band.window - 1);
            let hi = (at + 3).min(n - band.window + 1);
            for start in lo..hi {
                let end = start + band.window;
                if !self.overlaps_editable(start, end) {
                    continue;
                }
                let count = self.prefix_at(end) - self.prefix_at(start);
                if window_outside(count, band) {
                    self.collect_codons(start, end, primary, sites);
                }
            }
        }
    }

    fn collect_kmer_sites(&self, at: usize, primary: usize, sites: &mut BTreeSet<usize>) {
        let n = self.nt.len();
        for (rule_index, rule) in self.kmers.rules.iter().enumerate() {
            let k = rule.k;
            if k == 0 || k > 32 || n < k {
                continue;
            }
            let lo = at.saturating_sub(k - 1);
            let hi = (at + 3).min(n - k + 1);
            let mut keys = Vec::new();
            for pos in lo..hi {
                if let Some(key) = pack(&self.nt, pos, k, rule.canonical) {
                    let count = self.kmers.counts[rule_index]
                        .get(&key)
                        .copied()
                        .unwrap_or(0);
                    let touch = self.kmers.touching[rule_index]
                        .get(&key)
                        .copied()
                        .unwrap_or(0);
                    if count > rule.max_occurrences && touch > 0 {
                        keys.push(key);
                    }
                }
            }
            if keys.is_empty() {
                continue;
            }
            keys.sort_unstable();
            keys.dedup();
            for pos in 0..=n - k {
                if let Some(key) = pack(&self.nt, pos, k, rule.canonical) {
                    if keys.binary_search(&key).is_ok() {
                        self.collect_codons(pos, pos + k, primary, sites);
                    }
                }
            }
        }
    }

    fn collect_codons(
        &self,
        start: usize,
        end: usize,
        primary: usize,
        sites: &mut BTreeSet<usize>,
    ) {
        let cds_end = self.cds_off + self.seq.codon_count() * 3;
        if end <= self.cds_off || start >= cds_end || end <= start {
            return;
        }
        let from_nt = start.max(self.cds_off);
        let to_nt = end.min(cds_end);
        if to_nt <= from_nt {
            return;
        }
        let first = (from_nt - self.cds_off) / 3;
        let last = (to_nt - 1 - self.cds_off) / 3;
        for index in first..=last {
            if index != primary && !self.seq.is_locked(index) {
                sites.insert(index);
            }
        }
    }

    fn set_codon(&mut self, index: usize, codon: Codon) -> Result<Codon> {
        let old = self.seq.substitute(index, codon)?;
        if old == codon {
            return Ok(old);
        }
        let at = self.cds_off + index * 3;
        let bases = codon.to_bytes();
        let homo_span = homo_window(&self.nt, at, bases);
        self.homo_penalty = self.homo_after(at, bases);
        self.gc_penalty += self.gc_delta(at, bases);
        self.revise_gc_offenders(at, bases);
        let kmer_started = seed_profile_on().then(Instant::now);
        self.kmers.replace(&mut self.nt, at, bases);
        if let Some(kmer_started) = kmer_started {
            SEED_PROFILE.kmer_replace.record(kmer_started);
        }
        let prefix_started = seed_profile_on().then(Instant::now);
        self.shift_prefix(at, old, bases);
        if let Some(prefix_started) = prefix_started {
            SEED_PROFILE.prefix.record(prefix_started);
        }
        self.patterns.retile(
            PatternView {
                checker: self.checker,
                nt: &self.nt,
                cds_off: self.cds_off,
                codon_count: self.seq.codon_count(),
                editable: &self.editable,
            },
            at,
            homo_span,
        );
        Ok(old)
    }

    fn commit(&mut self, index: usize, codon: Codon) -> Result<()> {
        let old = self.set_codon(index, codon)?;
        if old != codon {
            self.edits += 1;
        }
        // A committed edit is not a trial. Materializing keeps the add list
        // at the depth of the edits currently being probed.
        self.flush_prefix();
        Ok(())
    }

    fn gc_delta(&self, at: usize, bases: [u8; 3]) -> f64 {
        let started = seed_profile_on().then(Instant::now);
        let delta = self.gc_delta_inner(at, bases);
        if let Some(started) = started {
            SEED_PROFILE.gc_delta.record(started);
        }
        delta
    }

    fn gc_delta_inner(&self, at: usize, bases: [u8; 3]) -> f64 {
        let mut delta = 0.0;
        let n = self.nt.len();
        for band in &self.bands {
            if n < band.window {
                continue;
            }
            let lo = at.saturating_sub(band.window - 1);
            let hi = (at + 3).min(n - band.window + 1);
            for start in lo..hi {
                let end = start + band.window;
                if !self.overlaps_editable(start, end) {
                    continue;
                }
                let before = self.prefix_at(end) - self.prefix_at(start);
                let after = adjusted_gc(before, &self.nt, at, bases, start, end);
                delta += band_penalty(after, band) - band_penalty(before, band);
            }
        }
        delta
    }

    /// Record the GC-count change of a codon without walking the suffix.
    ///
    /// `prefix[i]` is the GC count of `nt[..i]` at the last flush. A base
    /// whose GC bit flips adds that delta to every later prefix index. Those
    /// adds sit in `prefix_adds` until [`Self::flush_prefix`], so a trial can
    /// change a codon and put it back without rewriting the whole construct.
    /// Counts stay exact integers.
    fn shift_prefix(&mut self, at: usize, old: Codon, bases: [u8; 3]) {
        let old_bases = old.to_bytes();
        for offset in 0..3 {
            let delta = gc_bit(bases[offset]) as i32 - gc_bit(old_bases[offset]) as i32;
            if delta != 0 {
                self.push_suffix_add(at + 1 + offset, delta);
            }
        }
    }

    fn prefix_at(&self, index: usize) -> u32 {
        prefix_value(&self.prefix, &self.prefix_adds, index)
    }

    fn suffix_extra(&self, index: usize) -> i32 {
        let mut total = 0i32;
        for &(from, delta) in &self.prefix_adds {
            if from <= index {
                total += delta;
            }
        }
        total
    }

    fn push_suffix_add(&mut self, from: usize, delta: i32) {
        if delta == 0 {
            return;
        }
        if let Some(index) = self
            .prefix_adds
            .iter()
            .position(|(start, _)| *start == from)
        {
            let next = self.prefix_adds[index].1 + delta;
            if next == 0 {
                self.prefix_adds.swap_remove(index);
            } else {
                self.prefix_adds[index].1 = next;
            }
        } else {
            self.prefix_adds.push((from, delta));
        }
    }

    fn flush_prefix(&mut self) {
        if self.prefix_adds.is_empty() {
            return;
        }
        for index in 1..=self.nt.len() {
            let extra = self.suffix_extra(index);
            if extra != 0 {
                let next = self.prefix[index] as i32 + extra;
                debug_assert!(next >= 0, "gc prefix went negative");
                self.prefix[index] = next as u32;
            }
        }
        self.prefix_adds.clear();
    }

    fn finish_prefix(&mut self) {
        self.flush_prefix();
        debug_assert!(
            self.prefix_sound(),
            "gc prefix drifted from the nucleotide buffer"
        );
    }

    fn prefix_sound(&self) -> bool {
        let rebuilt = gc_prefix(&self.nt);
        rebuilt.len() == self.prefix.len()
            && rebuilt
                .iter()
                .enumerate()
                .all(|(index, count)| self.prefix_at(index) == *count)
    }

    fn recompute_gc(&mut self) {
        self.flush_prefix();
        self.gc_penalty = gc_penalty_of(&self.nt, &self.prefix, &self.bands, &self.editable);
    }

    fn homo_touching(&self, at: usize, bases: Option<[u8; 3]>) -> u64 {
        let started = seed_profile_on().then(Instant::now);
        let total = match bases {
            None => self.homo_penalty,
            Some(bases) => self.homo_after(at, bases),
        };
        if let Some(started) = started {
            SEED_PROFILE.homo.record(started);
        }
        total
    }

    /// Homopolymer excess after a virtual codon substitution.
    ///
    /// The excess is the same integer a full scan of the edited sequence
    /// returns. Only runs that intersect the codon can change, and those runs
    /// are bounded by walking out from the codon while the base repeats.
    fn homo_after(&self, at: usize, bases: [u8; 3]) -> u64 {
        let Some((left, right)) = homo_window(&self.nt, at, bases) else {
            return homo_total(&self.nt, at, Some(bases), &self.homo_rules, &self.editable);
        };
        let old = homo_between(
            &self.nt,
            at,
            None,
            left,
            right,
            &self.homo_rules,
            &self.editable,
        );
        let new = homo_between(
            &self.nt,
            at,
            Some(bases),
            left,
            right,
            &self.homo_rules,
            &self.editable,
        );
        match (old, new) {
            (Some(old), Some(new)) => {
                let next = self.homo_penalty as i64 + new as i64 - old as i64;
                if next < 0 {
                    homo_total(&self.nt, at, Some(bases), &self.homo_rules, &self.editable)
                } else {
                    next as u64
                }
            }
            _ => homo_total(&self.nt, at, Some(bases), &self.homo_rules, &self.editable),
        }
    }

    fn overlaps_editable(&self, start: usize, end: usize) -> bool {
        self.editable.iter().any(|&(lo, hi)| start < hi && lo < end)
    }

    fn view(&self) -> ConstructView<'_> {
        let start = self.cds_off;
        let end = start + self.seq.codon_count() * 3;
        ConstructView::new(&self.nt, vec![(start, end)]).with_editable_spans(self.editable.clone())
    }

    /// Whether the checker would accept this sequence.
    ///
    /// Equal to `ConstraintChecker::is_feasible` on this solver's view. The
    /// checker reports a candidate infeasible exactly when some blocking
    /// violation is actionable. Phase 1 puts every such violation in one of
    /// three buckets, and each bucket has one incremental source:
    ///
    /// * GC windows. `band.offenders` counts windows `window_outside` would
    ///   emit, maintained by an integer delta on each codon edit. That is the
    ///   checker's actionable test. The float `gc_penalty` still drives the
    ///   sweeps, but it is not the feasibility bit: after incremental edits it
    ///   can sit in `(0, 1e-9]` and a full re-scan of every window was the
    ///   cost of trusting it.
    /// * K-mer repeats. `KmerBook` packs the same overlapping words the
    ///   checker does, including flanks, and `actionable_excess` counts only
    ///   keys that exceed `max_occurrences` and have an occurrence touching an
    ///   editable span. A flank repeat that never touches the CDS keeps
    ///   `excess > 0` and `actionable_excess == 0`, so it does not flip this
    ///   predicate. Raw `excess` stays the repair objective.
    /// * Patterns (motifs, restriction sites, codon bans, regex, homopolymers).
    ///   `pattern_count` is the blocking actionable subset of `check_patterns`.
    ///   Motifs and homopolymers are kept current by a windowed rescan; a
    ///   blocking regular expression still scans the construct.
    ///
    /// Scope is not a fourth bucket. The checker applies CDS scope to motifs
    /// and ignores it for k-mers, GC, and homopolymers; the three sources
    /// above do the same. Advisory and disabled rules are in none of them,
    /// and `is_feasible` ignores those too. Undetectable rules are listed on
    /// the report rather than turned into violations.
    ///
    /// A construct that is not `ACGT` cannot use the GC proof (`fraction`
    /// drops masked bases from the denominator and the prefix does not), so
    /// that case asks the checker. Codon bytes and normalized flanks are
    /// `ACGT`, which is every sequence this solver builds.
    fn is_feasible(&self) -> bool {
        let started = seed_profile_on().then(Instant::now);
        let feasible = if self.acgt {
            // Offender counts are the integer windows `window_outside` would
            // emit. The float penalty can sit in (0, 1e-9] after incremental
            // edits and used to re-scan every window on each call.
            self.gc_offender_total() == 0
                && self.kmers.actionable_excess() == 0
                && self.pattern_count() == 0
        } else {
            let view = self.view();
            let check_started = seed_profile_on().then(Instant::now);
            let feasible = self.checker.is_feasible(&view);
            if let Some(check_started) = check_started {
                SEED_PROFILE.checker_check.record(check_started);
            }
            feasible
        };
        if let Some(started) = started {
            SEED_PROFILE.is_feasible.record(started);
        }
        feasible
    }

    /// Blocking actionable counts, keyed by constraint id.
    ///
    /// The same multiset `ConstraintChecker::check` reports for this view:
    /// pattern hits, one entry per over-limit k-mer key that touches an
    /// editable span, and one entry per GC window outside its band that
    /// overlaps an editable span. Ids with a zero count are absent.
    fn actionable_tally(&self) -> BTreeMap<String, usize> {
        let mut counts = BTreeMap::new();
        for (_start, _end, id) in self.pattern_spans() {
            *counts.entry(id).or_insert(0) += 1;
        }
        for (rule, offenders) in self.kmers.rules.iter().zip(self.kmers.offenders.iter()) {
            if *offenders > 0 {
                *counts.entry(rule.id.clone()).or_insert(0) += *offenders as usize;
            }
        }
        for band in &self.bands {
            if band.offenders > 0 {
                *counts.entry(band.id.clone()).or_insert(0) += band.offenders as usize;
            }
        }
        counts
    }

    fn revise_gc_offenders(&mut self, at: usize, bases: [u8; 3]) {
        let n = self.nt.len();
        for band in &mut self.bands {
            if band.window == 0 || n < band.window {
                continue;
            }
            let lo = at.saturating_sub(band.window - 1);
            let hi = (at + 3).min(n - band.window + 1);
            let mut delta: i32 = 0;
            for start in lo..hi {
                let end = start + band.window;
                if !spans_overlap(&self.editable, start, end) {
                    continue;
                }
                let before = prefix_value(&self.prefix, &self.prefix_adds, end)
                    - prefix_value(&self.prefix, &self.prefix_adds, start);
                let after = adjusted_gc(before, &self.nt, at, bases, start, end);
                delta += i32::from(window_outside(after, band))
                    - i32::from(window_outside(before, band));
            }
            let next = band.offenders as i32 + delta;
            debug_assert!(next >= 0, "gc offender count went negative");
            band.offenders = next.max(0) as u32;
        }
    }

    fn pattern_count(&self) -> usize {
        let started = seed_profile_on().then(Instant::now);
        let count = if self.patterns.incremental {
            self.patterns.violations.len()
        } else {
            self.pattern_spans().len()
        };
        if let Some(started) = started {
            SEED_PROFILE.pattern_count.record(started);
        }
        count
    }

    fn pattern_violations(&self) -> Vec<(usize, usize, usize)> {
        let started = seed_profile_on().then(Instant::now);
        let hits = self.pattern_spans();
        let count = hits.len();
        let found = hits
            .into_iter()
            .map(|(start, end, _id)| (start, end, count))
            .collect();
        if let Some(started) = started {
            SEED_PROFILE.pattern_violations.record(started);
        }
        found
    }

    /// Blocking actionable pattern hits, in checker order.
    ///
    /// The incremental ledger is the motif and homopolymer subset. It matches
    /// `check_patterns` filtered to blocking actionable pattern rules whenever
    /// no blocking regular expression is compiled. With one, the checker is
    /// scanned, because a regex match is not confined to `radius`.
    fn pattern_spans(&self) -> Vec<(usize, usize, String)> {
        if self.patterns.incremental {
            return self.patterns.violations.clone();
        }
        let view = self.view();
        let check_started = seed_profile_on().then(Instant::now);
        let violations = self.checker.check_patterns(&view);
        if let Some(check_started) = check_started {
            SEED_PROFILE.checker_check.record(check_started);
        }
        violations
            .into_iter()
            .filter(|violation| {
                violation.class == ConstraintClass::Blocking
                    && violation.actionable
                    && is_pattern(&self.checker.set().constraints, &violation.constraint_id)
            })
            .map(|violation| {
                (
                    violation.start_nt,
                    violation.end_nt,
                    violation.constraint_id,
                )
            })
            .collect()
    }

    fn hot_codons(&self) -> Vec<usize> {
        let started = seed_profile_on().then(Instant::now);
        let indexes = self.hot_codons_inner();
        if let Some(started) = started {
            SEED_PROFILE.hot_codons.record(started);
        }
        indexes
    }

    fn hot_codons_inner(&self) -> Vec<usize> {
        let mut indexes = Vec::new();
        for rule in &self.kmers.rules {
            if self.nt.len() < rule.k {
                continue;
            }
            // Both copies. Marking only the occurrence that pushed the count
            // over the limit left the earlier copy out of the pair search, and
            // the edit that breaks the repeat can sit on that earlier copy.
            let mut seen: HashMap<u64, Vec<usize>> = HashMap::new();
            for pos in 0..=self.nt.len() - rule.k {
                let Some(key) = pack(&self.nt, pos, rule.k, rule.canonical) else {
                    continue;
                };
                seen.entry(key).or_default().push(pos);
            }
            let last = self.seq.codon_count().saturating_sub(1);
            for positions in seen.values() {
                if positions.len() <= rule.max_occurrences as usize {
                    continue;
                }
                for pos in positions {
                    let cds_pos = pos.saturating_sub(self.cds_off);
                    let start_codon = cds_pos / 3;
                    let end_codon = (cds_pos + rule.k - 1) / 3;
                    for codon in start_codon..=end_codon.min(last) {
                        indexes.push(codon);
                    }
                }
            }
        }
        indexes.sort_unstable();
        indexes.dedup();
        indexes
    }

    fn residuals(&self) -> Vec<SeedResidual> {
        let mut counts: BTreeMap<String, usize> = BTreeMap::new();
        for violation in self.checker.check(&self.view()) {
            if violation.class == ConstraintClass::Blocking && violation.actionable {
                *counts.entry(violation.constraint_id).or_insert(0) += 1;
            }
        }
        counts
            .into_iter()
            .map(|(constraint_id, count)| SeedResidual {
                kind: kind_of(self.checker, &constraint_id),
                constraint_id,
                count,
            })
            .collect()
    }
}

pub fn unsearched_rules(checker: &ConstraintChecker) -> Vec<UnsearchedRule> {
    checker
        .undetectable()
        .iter()
        .map(|id| UnsearchedRule {
            reason: checker
                .set()
                .constraints
                .iter()
                .find(|constraint| &constraint.id == id)
                .map(|constraint| unsearched_reason(&constraint.rule))
                .unwrap_or_else(|| "this rule is not searched in this version".to_string()),
            id: id.clone(),
        })
        .collect()
}

fn unsearched_reason(rule: &Rule) -> String {
    match rule {
        Rule::LocalStructure { .. } => {
            "local_structure is not searched in this version; it needs the Phase 2 structure engine"
                .to_string()
        }
        Rule::MinCodonFrequency { .. } => {
            "min_codon_frequency is not searched in this version".to_string()
        }
        _ => "this rule is not searched in this version".to_string(),
    }
}

fn kind_of(checker: &ConstraintChecker, id: &str) -> String {
    checker
        .set()
        .constraints
        .iter()
        .find(|constraint| constraint.id == id)
        .map(|constraint| rule_kind(&constraint.rule).to_string())
        .unwrap_or_else(|| "unknown".to_string())
}

fn rule_kind(rule: &Rule) -> &'static str {
    match rule {
        Rule::ForbiddenMotif { .. } => "forbidden_motif",
        Rule::RestrictionSite { .. } => "restriction_site",
        Rule::Homopolymer { .. } => "homopolymer",
        Rule::KmerRepeat { .. } => "kmer_repeat",
        Rule::GcWindow { .. } => "gc_window",
        Rule::LocalStructure { .. } => "local_structure",
        Rule::CodonBan { .. } => "codon_ban",
        Rule::MinCodonFrequency { .. } => "min_codon_frequency",
        Rule::CustomRegex { .. } => "custom_regex",
    }
}

fn is_pattern(constraints: &[crate::spec::Constraint], id: &str) -> bool {
    constraints.iter().any(|constraint| {
        constraint.id == id
            && matches!(
                constraint.rule,
                Rule::ForbiddenMotif { .. }
                    | Rule::RestrictionSite { .. }
                    | Rule::Homopolymer { .. }
                    | Rule::CustomRegex { .. }
                    | Rule::CodonBan { .. }
            )
    })
}

fn gc_bands(checker: &ConstraintChecker) -> Vec<GcBand> {
    checker
        .set()
        .constraints
        .iter()
        .filter(|constraint| constraint.enabled && constraint.class == ConstraintClass::Blocking)
        .filter_map(|constraint| match &constraint.rule {
            Rule::GcWindow {
                window_nt,
                min,
                max,
            } => Some(GcBand {
                id: constraint.id.clone(),
                window: *window_nt as usize,
                min: *min,
                max: *max,
                offenders: 0,
            }),
            _ => None,
        })
        .collect()
}

fn kmer_rules(checker: &ConstraintChecker) -> Vec<KmerCfg> {
    checker
        .set()
        .constraints
        .iter()
        .filter(|constraint| constraint.enabled && constraint.class == ConstraintClass::Blocking)
        .filter_map(|constraint| match &constraint.rule {
            Rule::KmerRepeat {
                k,
                max_occurrences,
                canonical_strand,
                ..
            } => Some(KmerCfg {
                id: constraint.id.clone(),
                k: *k as usize,
                max_occurrences: *max_occurrences,
                canonical: *canonical_strand,
            }),
            _ => None,
        })
        .collect()
}

fn homo_rules(checker: &ConstraintChecker) -> Vec<(Option<u8>, usize)> {
    use crate::spec::HomopolymerBase;
    checker
        .set()
        .constraints
        .iter()
        .filter(|constraint| constraint.enabled && constraint.class == ConstraintClass::Blocking)
        .filter_map(|constraint| match &constraint.rule {
            Rule::Homopolymer { base, max_run } => Some((
                match base {
                    HomopolymerBase::Any => None,
                    HomopolymerBase::A => Some(b'A'),
                    HomopolymerBase::C => Some(b'C'),
                    HomopolymerBase::G => Some(b'G'),
                    HomopolymerBase::T => Some(b'T'),
                },
                *max_run as usize,
            )),
            _ => None,
        })
        .collect()
}

fn editable_spans_of(sequence: &CodonSequence) -> Vec<(usize, usize)> {
    let offset = sequence.five_prime_flank().len();
    let mut spans: Vec<(usize, usize)> = Vec::new();
    for index in 0..sequence.codon_count() {
        if sequence.is_locked(index) {
            continue;
        }
        let start = offset + index * 3;
        match spans.last_mut() {
            Some(last) if last.1 == start => last.1 = start + 3,
            _ => spans.push((start, start + 3)),
        }
    }
    spans
}

fn gc_prefix(nt: &[u8]) -> Vec<u32> {
    let mut prefix = Vec::with_capacity(nt.len() + 1);
    prefix.push(0);
    for (index, base) in nt.iter().enumerate() {
        prefix.push(prefix[index] + gc_bit(*base));
    }
    prefix
}

fn gc_bit(base: u8) -> u32 {
    u32::from(matches!(base, b'G' | b'C'))
}

fn prefix_value(prefix: &[u32], adds: &[(usize, i32)], index: usize) -> u32 {
    let mut extra = 0i32;
    for &(from, delta) in adds {
        if from <= index {
            extra += delta;
        }
    }
    let value = prefix[index] as i32 + extra;
    debug_assert!(value >= 0, "gc prefix went negative");
    value as u32
}

fn gc_penalty_of(nt: &[u8], prefix: &[u32], bands: &[GcBand], editable: &[(usize, usize)]) -> f64 {
    let mut penalty = 0.0;
    for band in bands {
        if nt.len() < band.window {
            continue;
        }
        for start in 0..=nt.len() - band.window {
            let end = start + band.window;
            if !editable.iter().any(|&(lo, hi)| start < hi && lo < end) {
                continue;
            }
            penalty += band_penalty(prefix[end] - prefix[start], band);
        }
    }
    penalty
}

fn band_penalty(count: u32, band: &GcBand) -> f64 {
    let fraction = count as f64 / band.window as f64;
    if fraction > band.max {
        fraction - band.max
    } else if fraction < band.min {
        band.min - fraction
    } else {
        0.0
    }
}

/// A window the checker would emit: strict `fraction < min` or `fraction > max`.
///
/// On an `ACGT` construct the prefix count over `window` bases is the same
/// fraction `GcIndex::fraction` returns, so this predicate matches
/// `windows_outside`.
fn window_outside(count: u32, band: &GcBand) -> bool {
    let fraction = count as f64 / band.window as f64;
    fraction < band.min || fraction > band.max
}

fn count_gc_offenders(
    nt: &[u8],
    prefix: &[u32],
    band: &GcBand,
    editable: &[(usize, usize)],
) -> u32 {
    if band.window == 0 || nt.len() < band.window {
        return 0;
    }
    let mut count = 0u32;
    for start in 0..=nt.len() - band.window {
        let end = start + band.window;
        if spans_overlap(editable, start, end) && window_outside(prefix[end] - prefix[start], band)
        {
            count += 1;
        }
    }
    count
}

fn adjusted_gc(before: u32, nt: &[u8], at: usize, bases: [u8; 3], start: usize, end: usize) -> u32 {
    let mut count = before;
    for (offset, base) in bases.iter().enumerate() {
        let pos = at + offset;
        if pos >= start && pos < end {
            count = count - gc_bit(nt[pos]) + gc_bit(*base);
        }
    }
    count
}

fn spans_overlap(editable: &[(usize, usize)], start: usize, end: usize) -> bool {
    editable.iter().any(|&(lo, hi)| start < hi && lo < end)
}

fn homo_total(
    nt: &[u8],
    at: usize,
    bases: Option<[u8; 3]>,
    rules: &[(Option<u8>, usize)],
    editable: &[(usize, usize)],
) -> u64 {
    let mut total = 0u64;
    for &(base, max_run) in rules {
        let mut pos = 0;
        while pos < nt.len() {
            let here = base_at(nt, at, bases, pos);
            if base.is_some_and(|wanted| wanted != here)
                || !matches!(here, b'A' | b'C' | b'G' | b'T')
            {
                pos += 1;
                continue;
            }
            let mut end = pos + 1;
            while end < nt.len() && base_at(nt, at, bases, end) == here {
                end += 1;
            }
            if end - pos > max_run && spans_overlap(editable, pos, end) {
                total += (end - pos - max_run) as u64;
            }
            pos = end;
        }
    }
    total
}

/// Inclusive-exclusive window that contains every homopolymer run intersecting
/// the codon, in both the current sequence and the sequence with `bases`
/// written at `at`. `None` when the codon does not fit.
fn homo_window(nt: &[u8], at: usize, bases: [u8; 3]) -> Option<(usize, usize)> {
    if at.checked_add(3)? > nt.len() {
        return None;
    }
    let mut left = at;
    for base in [nt[at], bases[0]] {
        let mut pos = at;
        while pos > 0 && nt[pos - 1] == base {
            pos -= 1;
        }
        left = left.min(pos);
    }
    let mut right = at + 3;
    for base in [nt[at + 2], bases[2]] {
        let mut pos = at + 2;
        while pos + 1 < nt.len() && nt[pos + 1] == base {
            pos += 1;
        }
        right = right.max(pos + 1);
    }
    Some((left, right))
}

/// Penalty of runs that lie inside `[start, end)`.
///
/// `None` when a run crosses the window. That means the window is too small
/// and the caller must scan the whole sequence instead of trusting a delta.
fn homo_between(
    nt: &[u8],
    at: usize,
    bases: Option<[u8; 3]>,
    start: usize,
    end: usize,
    rules: &[(Option<u8>, usize)],
    editable: &[(usize, usize)],
) -> Option<u64> {
    let mut total = 0u64;
    for &(base, max_run) in rules {
        let mut pos = start;
        while pos < end {
            let here = base_at(nt, at, bases, pos);
            if base.is_some_and(|wanted| wanted != here)
                || !matches!(here, b'A' | b'C' | b'G' | b'T')
            {
                pos += 1;
                continue;
            }
            let mut run_end = pos + 1;
            while run_end < nt.len() && base_at(nt, at, bases, run_end) == here {
                run_end += 1;
            }
            if run_end > end {
                return None;
            }
            if run_end - pos > max_run && spans_overlap(editable, pos, run_end) {
                total += (run_end - pos - max_run) as u64;
            }
            pos = run_end;
        }
    }
    Some(total)
}

fn base_at(nt: &[u8], at: usize, bases: Option<[u8; 3]>, pos: usize) -> u8 {
    if let Some(bases) = bases {
        if pos >= at && pos < at + 3 {
            return bases[pos - at];
        }
    }
    nt[pos]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::{ConstraintSet, RepairPolicy, Scope};

    fn checker(rules: Vec<(&str, Rule)>) -> ConstraintChecker {
        let set = ConstraintSet {
            schema_version: Some("1.0.0".into()),
            constraint_set_id: "cst_seed".into(),
            name: None,
            derived_from: vec![],
            constraints: rules
                .into_iter()
                .map(|(id, rule)| crate::spec::Constraint {
                    id: id.into(),
                    class: ConstraintClass::Blocking,
                    scope: Scope::Cds,
                    cistron_id: None,
                    weight: 1.0,
                    enabled: true,
                    rule,
                })
                .collect(),
            repair_policy: RepairPolicy::default(),
            digest: None,
        };
        ConstraintChecker::compile(&set).unwrap()
    }

    fn weights_one() -> [f64; 64] {
        [1.0; 64]
    }

    fn tiny_sequence() -> CodonSequence {
        CodonSequence::from_cds("ATGAAATAA", 11).unwrap()
    }

    #[test]
    fn a_past_deadline_is_budget_exceeded_and_a_future_one_still_resolves() {
        let sequence = tiny_sequence();
        let compiled = checker(vec![]);
        let weights = weights_one();
        let past = SpendLimit {
            deadline: Some(std::time::Instant::now() - std::time::Duration::from_secs(1)),
            max_rss_bytes: None,
        };
        let err = resolve_feasible_seed_limited(&sequence, &compiled, &weights, past)
            .expect_err("a deadline that has already passed must refuse");
        assert_eq!(err.code, "CF-BUDGET-003");
        assert!(err.detail.contains("wall-clock"));

        let future = SpendLimit {
            deadline: Some(std::time::Instant::now() + std::time::Duration::from_secs(60)),
            max_rss_bytes: None,
        };
        let report = resolve_feasible_seed_limited(&sequence, &compiled, &weights, future).unwrap();
        assert!(report.feasible);

        let tight = SpendLimit {
            deadline: None,
            max_rss_bytes: Some(1),
        };
        let err = resolve_feasible_seed_limited(&sequence, &compiled, &weights, tight)
            .expect_err("a one-byte resident cap must refuse on this process");
        assert_eq!(err.code, "CF-BUDGET-003");
        assert!(err.detail.contains("memory cap"));
    }

    fn weight_table(placed: &CodonSequence, overrides: &[(&str, f64)]) -> [f64; 64] {
        let mut weights = [0.05; 64];
        for codon in placed.codons() {
            weights[codon.index_usize()] = 1.0;
        }
        for (label, weight) in overrides {
            weights[Codon::parse(label).unwrap().index_usize()] = *weight;
        }
        weights
    }

    #[test]
    fn unit_seed_climb_repairs_a_gc_window_outside_the_local_neighbourhood() {
        // Codon 1 wants CTG (GC 2, weight 1) over TTA (GC 0, weight 0.3).
        // That step pushes the 24 nt window of codons 1..=8 over 0.5. The only
        // codon that can give two GC bases back is Arg at index 8, seven
        // positions away. Alanines in between are GC-invariant. DEP-SRCH-018
        // stops at four codons, so this stays on TTA unless the conflict climb
        // runs.
        let cds = format!(
            "ATG{}{}{}{}TAA",
            "TTA",
            "GCT".repeat(4),
            "TTT".repeat(2),
            "CGG"
        );
        let sequence = CodonSequence::from_cds(&cds, 11).unwrap();
        let weights = weight_table(
            &sequence,
            &[
                ("CTG", 1.0),
                ("TTA", 0.3),
                ("TTG", 0.1),
                ("CTT", 0.1),
                ("CTC", 0.1),
                ("CTA", 0.1),
                ("AGA", 0.97),
                ("AGG", 0.1),
                ("CGT", 0.1),
                ("CGC", 0.1),
                ("CGA", 0.1),
                ("TTC", 0.5),
                ("GCC", 0.2),
                ("GCA", 0.2),
                ("GCG", 0.2),
            ],
        );
        let compiled = checker(vec![(
            "c_gc",
            Rule::GcWindow {
                window_nt: 24,
                min: 0.0,
                max: 0.5,
            },
        )]);
        let report = resolve_feasible_seed(&sequence, &compiled, &weights).unwrap();
        assert!(report.feasible, "{:?}", report.residual);
        assert_eq!(report.sequence.protein(), sequence.protein());
        let codons: Vec<String> = report
            .sequence
            .codons()
            .iter()
            .map(|codon| codon.to_string_upper())
            .collect();
        assert_eq!(
            codons[1], "CTG",
            "high-weight leucine was not restored: {codons:?}"
        );
        assert_eq!(
            codons[8], "AGA",
            "distant arginine did not pay the GC cost: {codons:?}"
        );
        let sacrificed = cai_of(&sequence, &weights);
        assert!(
            cai_of(&report.sequence, &weights) > sacrificed,
            "climb did not raise CAI above the low-weight leucine"
        );
    }

    #[test]
    fn unit_seed_climb_breaks_the_distant_copy_of_a_12mer() {
        // AAACCCGGGTTT twice. The second copy is shifted by one nucleotide, so
        // the three overlapping 12-mers share an arginine (CGG) whose CGT
        // synonym breaks all three for a small weight drop. Every synonym
        // inside the first copy is expensive. Greedy repair edits the first
        // lysine. A partner search limited to four codons cannot see the
        // arginine.
        let fillers = [
            "GAT", "GAC", "GAA", "GAG", "CAT", "CAC", "TAT", "TAC", "TGG", "TGT", "TGC", "ATT",
        ];
        let mut codons = vec!["ATG".to_string(), "CTA".to_string()];
        codons.extend(["AAA", "CCC", "GGG", "TTT"].map(str::to_string));
        while codons.len() < 30 {
            codons.push(fillers[(codons.len() - 6) % fillers.len()].to_string());
        }
        assert_eq!(codons.len(), 30);
        // Codon 30 starts one base before the second copy, so the word is not
        // on the same codon boundaries as the first copy.
        let mut nt = codons.concat();
        assert_eq!(nt.len(), 90);
        nt.push('A');
        nt.push_str("AAACCCGGGTTT");
        nt.push_str("TT");
        nt.push_str("TAA");
        assert_eq!(nt.len() % 3, 0, "{}", nt.len());
        let sequence = CodonSequence::from_cds(&nt, 11).unwrap();
        let weights = weight_table(
            &sequence,
            &[
                ("AAG", 0.05),
                ("CCT", 0.05),
                ("CCA", 0.05),
                ("CCG", 0.05),
                ("GGA", 0.05),
                ("GGT", 0.05),
                ("GGC", 0.05),
                ("TTC", 0.05),
                ("CGT", 0.99),
                ("CGC", 0.05),
                ("CGA", 0.05),
                ("AGA", 0.05),
                ("AGG", 0.05),
            ],
        );
        let compiled = checker(vec![(
            "c_kmer",
            Rule::KmerRepeat {
                k: 12,
                max_occurrences: 1,
                count_overlapping: true,
                canonical_strand: false,
            },
        )]);
        let report = resolve_feasible_seed(&sequence, &compiled, &weights).unwrap();
        assert!(report.feasible, "{:?}", report.residual);
        assert_eq!(report.sequence.protein(), sequence.protein());
        let designed: Vec<String> = report
            .sequence
            .codons()
            .iter()
            .map(|codon| codon.to_string_upper())
            .collect();
        assert_eq!(
            designed[2], "AAA",
            "first lysine was sacrificed instead of the distant arginine: {designed:?}"
        );
        assert_eq!(
            &designed[3..6],
            &["CCC".to_string(), "GGG".to_string(), "TTT".to_string()]
        );
        assert_eq!(
            designed[32], "CGT",
            "distant arginine did not break the repeat: {designed:?}"
        );
    }

    #[test]
    fn unit_seed_repair_clears_a_gc_window_and_keeps_the_protein() {
        // Leucine CTG is 2/3 GC. TTA is 0/3, so a 12-nt window can enter 0.3–0.6.
        // Alanine cannot: every synonym is at least 2/3 GC, and four of them
        // already exceed a 0.6 ceiling.
        let cds = format!("ATG{}TAA", "CTG".repeat(12));
        let sequence = CodonSequence::from_cds(&cds, 11).unwrap();
        let compiled = checker(vec![(
            "c_gc",
            Rule::GcWindow {
                window_nt: 12,
                min: 0.3,
                max: 0.6,
            },
        )]);
        let report = resolve_feasible_seed(&sequence, &compiled, &weights_one()).unwrap();
        assert!(
            report.feasible,
            "seed repair returned the GC problem unresolved: {:?}",
            report.residual
        );
        assert_ne!(report.sequence.cds_string(), cds);
        assert_eq!(report.sequence.protein(), sequence.protein());
    }

    #[test]
    fn unit_seed_repair_breaks_a_custom_regex_cryptic_start() {
        // ATG|AGG|AGG|AAA|AAT|GTA|TAA contains AGGAGGAAAAATG.
        let cds = "ATGAGGAGGAAAAATGTATAA";
        let sequence = CodonSequence::from_cds(cds, 11).unwrap();
        let compiled = checker(vec![(
            "c_cryptic_start",
            Rule::CustomRegex {
                pattern: "AGGAGG[ACGT]{4,12}ATG".into(),
                max_occurrences: 0,
                flags: vec!["both_strands".into()],
            },
        )]);
        let report = resolve_feasible_seed(&sequence, &compiled, &weights_one()).unwrap();
        assert!(
            report.feasible,
            "cryptic-start regex was not repaired: {:?}",
            report.residual
        );
        assert_eq!(report.sequence.protein(), sequence.protein());
        assert!(!report.sequence.cds_string().contains("AGGAGGAAAAATG"));
    }

    #[test]
    fn pattern_check_matches_the_filtered_full_check_and_skips_kmer_and_gc() {
        use crate::spec::{HomopolymerBase, Strand};

        let compiled = checker(vec![
            (
                "c_gc",
                Rule::GcWindow {
                    window_nt: 6,
                    min: 0.0,
                    max: 0.2,
                },
            ),
            (
                "c_kmer",
                Rule::KmerRepeat {
                    k: 4,
                    max_occurrences: 1,
                    count_overlapping: true,
                    canonical_strand: false,
                },
            ),
            (
                "c_ecori",
                Rule::RestrictionSite {
                    enzyme: "EcoRI".into(),
                    site: "GAATTC".into(),
                    strand: Strand::Both,
                    max_occurrences: 0,
                },
            ),
            (
                "c_homo",
                Rule::Homopolymer {
                    base: HomopolymerBase::Any,
                    max_run: 3,
                },
            ),
        ]);
        let sequence = CodonSequence::from_cds("ATGGAATTCGGGGGGAAAAAATAA", 11).unwrap();
        let weights = weights_one();
        let solver = Solver::new(sequence, &compiled, &weights, TieBreak::Weight).unwrap();
        let view = solver.view();
        let project = |violation: &crate::detect::Violation| {
            (
                violation.constraint_id.clone(),
                violation.start_nt,
                violation.end_nt,
                violation.actionable,
                violation.class,
            )
        };
        let filtered: Vec<_> = compiled
            .check(&view)
            .into_iter()
            .filter(|violation| is_pattern(&compiled.set().constraints, &violation.constraint_id))
            .map(|violation| project(&violation))
            .collect();
        let patterns: Vec<_> = compiled
            .check_patterns(&view)
            .into_iter()
            .map(|violation| project(&violation))
            .collect();
        assert_eq!(patterns, filtered);
        assert!(patterns.iter().any(|hit| hit.0 == "c_ecori"));
        assert!(patterns.iter().any(|hit| hit.0 == "c_homo"));
        assert!(patterns
            .iter()
            .all(|hit| hit.0 != "c_gc" && hit.0 != "c_kmer"));
        let full = compiled.check(&view);
        assert!(full
            .iter()
            .any(|hit| hit.constraint_id == "c_gc" || hit.constraint_id == "c_kmer"));
        assert_eq!(solver.pattern_count(), filtered.len());
    }

    #[test]
    fn homopolymer_delta_matches_a_full_scan_including_merges() {
        use crate::spec::HomopolymerBase;

        let compiled = checker(vec![
            (
                "c_any",
                Rule::Homopolymer {
                    base: HomopolymerBase::Any,
                    max_run: 4,
                },
            ),
            (
                "c_a",
                Rule::Homopolymer {
                    base: HomopolymerBase::A,
                    max_run: 3,
                },
            ),
        ]);
        let block = [
            "AAA", "AAA", "CTG", "CCC", "CCC", "GCT", "TTT", "GGG", "AAC",
        ];
        let cds = format!("ATG{}TAA", block.repeat(6).join(""));
        assert_eq!(cds.len() % 3, 0, "{}", cds.len());
        let sequence = CodonSequence::from_cds(&cds, 11).unwrap();
        let weights = weights_one();
        let mut solver = Solver::new(sequence, &compiled, &weights, TieBreak::Weight).unwrap();
        for index in 0..solver.seq.codon_count() {
            let current = solver.seq.codons()[index];
            for &synonym in solver.seq.synonyms_at(index) {
                let at = solver.cds_off + index * 3;
                let bases = synonym.to_bytes();
                let expected = homo_total(
                    &solver.nt,
                    at,
                    Some(bases),
                    &solver.homo_rules,
                    &solver.editable,
                );
                assert_eq!(
                    solver.homo_after(at, bases),
                    expected,
                    "virtual edit {index} {current} -> {synonym}"
                );
                let (left, right) = homo_window(&solver.nt, at, bases).unwrap();
                assert!(
                    homo_between(
                        &solver.nt,
                        at,
                        None,
                        left,
                        right,
                        &solver.homo_rules,
                        &solver.editable,
                    )
                    .is_some(),
                    "old window crossed at {index}"
                );
                assert!(
                    homo_between(
                        &solver.nt,
                        at,
                        Some(bases),
                        left,
                        right,
                        &solver.homo_rules,
                        &solver.editable,
                    )
                    .is_some(),
                    "new window crossed at {index}"
                );
            }
        }
        let index = 4;
        let original = solver.seq.codons()[index];
        let replacement = solver
            .seq
            .synonyms_at(index)
            .iter()
            .copied()
            .find(|codon| *codon != original)
            .unwrap();
        solver.set_codon(index, replacement).unwrap();
        let rescanned = homo_total(&solver.nt, 0, None, &solver.homo_rules, &solver.editable);
        assert_eq!(solver.homo_penalty, rescanned);
        solver.set_codon(index, original).unwrap();
        let restored = homo_total(&solver.nt, 0, None, &solver.homo_rules, &solver.editable);
        assert_eq!(solver.homo_penalty, restored);
    }

    fn assert_feasible_matches(solver: &Solver<'_>, label: &str) {
        let view = solver.view();
        let violations = solver.checker.check(&view);
        let checker_feasible = !violations
            .iter()
            .any(|violation| violation.class == ConstraintClass::Blocking && violation.actionable);
        assert_eq!(
            solver.is_feasible(),
            checker_feasible,
            "{label}: incremental feasible={} gc={} actionable_kmer_excess={} raw_excess={} patterns={} violations={violations:?}",
            solver.is_feasible(),
            solver.gc_penalty,
            solver.kmers.actionable_excess(),
            solver.kmers.excess(),
            solver.pattern_count(),
        );
    }

    /// A repeat that never touches an editable base is not the candidate's
    /// problem. `excess` still counts it. `is_feasible` must not.
    #[test]
    fn flank_repeat_is_not_actionable_and_stays_feasible() {
        use crate::spec::{HomopolymerBase, Strand};

        let compiled = checker(vec![
            (
                "c_gc",
                Rule::GcWindow {
                    window_nt: 12,
                    min: 0.0,
                    max: 1.0,
                },
            ),
            (
                "c_kmer",
                Rule::KmerRepeat {
                    k: 8,
                    max_occurrences: 1,
                    count_overlapping: true,
                    canonical_strand: false,
                },
            ),
            (
                "c_homo",
                Rule::Homopolymer {
                    base: HomopolymerBase::Any,
                    max_run: 20,
                },
            ),
            (
                "c_ecori",
                Rule::RestrictionSite {
                    enzyme: "EcoRI".into(),
                    site: "GAATTC".into(),
                    strand: Strand::Both,
                    max_occurrences: 0,
                },
            ),
        ]);
        // `ACGTACGT` three times, all inside the 16 nt flank. The CDS is short
        // enough that none of its own 8-mers repeat, and the repeated flank
        // words end before the first editable base.
        let sequence = CodonSequence::from_cds("ATGGCTTGGCAGTAA", 11)
            .unwrap()
            .with_five_prime_flank("ACGTACGTACGTACGT")
            .unwrap();
        let weights = weights_one();
        let solver = Solver::new(sequence, &compiled, &weights, TieBreak::Weight).unwrap();
        assert!(
            solver.kmers.excess() > 0,
            "the flank repeat must be counted"
        );
        assert_eq!(solver.kmers.actionable_excess(), 0);
        assert!(
            solver.is_feasible(),
            "a flank-only repeat is not infeasible"
        );
        assert_feasible_matches(&solver, "flank repeat");
    }

    #[test]
    fn incremental_feasibility_matches_the_checker_over_random_edits() {
        use crate::spec::{HomopolymerBase, Strand};
        use cf_types::LockReason;

        let compiled = checker(vec![
            (
                "c_gc",
                Rule::GcWindow {
                    window_nt: 12,
                    min: 0.1,
                    max: 0.7,
                },
            ),
            (
                "c_kmer",
                Rule::KmerRepeat {
                    k: 8,
                    max_occurrences: 1,
                    count_overlapping: true,
                    canonical_strand: true,
                },
            ),
            (
                "c_homo",
                Rule::Homopolymer {
                    base: HomopolymerBase::Any,
                    max_run: 6,
                },
            ),
            (
                "c_ecori",
                Rule::RestrictionSite {
                    enzyme: "EcoRI".into(),
                    site: "GAATTC".into(),
                    strand: Strand::Both,
                    max_occurrences: 0,
                },
            ),
            (
                "c_motif",
                Rule::ForbiddenMotif {
                    motif: "CCC".into(),
                    strand: Strand::Both,
                    max_occurrences: 0,
                },
            ),
            (
                "c_regex",
                Rule::CustomRegex {
                    pattern: "CCA".into(),
                    max_occurrences: 0,
                    flags: vec![],
                },
            ),
        ]);
        // One of each editable family, so the 8-mers do not already repeat.
        // Alanine can enter the GC band, proline can enter CCC or CCA.
        let cds = "ATGGCTGCATGGCAGTTACGTGGACCTAAC TAA".replace(' ', "");
        let mut sequence = CodonSequence::from_cds(&cds, 11)
            .unwrap()
            .with_five_prime_flank("ACGTACGTACGTACGT")
            .unwrap()
            .with_three_prime_flank("TTAA")
            .unwrap();
        sequence.add_lock(5, 6, LockReason::UserLocked).unwrap();
        let weights = weights_one();
        let mut solver = Solver::new(sequence, &compiled, &weights, TieBreak::Weight).unwrap();
        assert_feasible_matches(&solver, "initial");
        assert!(
            solver.kmers.excess() > solver.kmers.actionable_excess(),
            "the flank repeat has to keep raw excess above the actionable excess"
        );

        let mut state = 0x1234_5678_9ABC_DEF0u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut saw_feasible = solver.is_feasible();
        let mut saw_infeasible = !saw_feasible;
        let mut saw_flank_only = solver.is_feasible() && solver.kmers.excess() > 0;
        let positions = solver.seq.codon_count();
        for step in 0..400 {
            let index = (next() as usize) % positions;
            if solver.seq.is_locked(index) {
                continue;
            }
            let synonyms = solver.seq.synonyms_at(index);
            let codon = synonyms[(next() as usize) % synonyms.len()];
            solver.set_codon(index, codon).unwrap();
            assert_feasible_matches(&solver, &format!("step {step} codon {index}"));
            if solver.is_feasible() {
                saw_feasible = true;
                if solver.kmers.excess() > 0 {
                    saw_flank_only = true;
                }
            } else {
                saw_infeasible = true;
            }
        }
        assert!(
            saw_feasible && saw_infeasible && saw_flank_only,
            "walk missed a verdict: feasible={saw_feasible} infeasible={saw_infeasible} flank_only={saw_flank_only}"
        );
    }

    fn checker_pattern_spans(solver: &Solver<'_>) -> Vec<(usize, usize, String)> {
        let view = solver.view();
        solver
            .checker
            .check_patterns(&view)
            .into_iter()
            .filter(|violation| {
                violation.class == ConstraintClass::Blocking
                    && violation.actionable
                    && is_pattern(&solver.checker.set().constraints, &violation.constraint_id)
            })
            .map(|violation| {
                (
                    violation.start_nt,
                    violation.end_nt,
                    violation.constraint_id,
                )
            })
            .collect()
    }

    /// The EcoRI site begins exactly `rescan_radius` bases before the edited
    /// codon, and the homopolymer merge crosses that codon. A window that is
    /// short by one base reports a different span list.
    #[test]
    fn windowed_pattern_spans_match_a_full_scan_over_random_edits() {
        use crate::spec::{HomopolymerBase, Strand};
        use cf_types::LockReason;

        let compiled = checker(vec![
            (
                "c_ecori",
                Rule::RestrictionSite {
                    enzyme: "EcoRI".into(),
                    site: "GAATTC".into(),
                    strand: Strand::Both,
                    max_occurrences: 0,
                },
            ),
            (
                "c_motif",
                Rule::ForbiddenMotif {
                    motif: "GCG".into(),
                    strand: Strand::Sense,
                    max_occurrences: 2,
                },
            ),
            (
                "c_homo",
                Rule::Homopolymer {
                    base: HomopolymerBase::Any,
                    max_run: 4,
                },
            ),
        ]);
        let weights = weights_one();
        // ATG CGA ATT TTA TAA. Codon 3 TTA -> CTG completes GAATTC at nt 4,
        // five bases before the codon, which is the EcoRI rescan radius.
        let boundary = CodonSequence::from_cds("ATGCGAATTTTATAA", 11).unwrap();
        let mut solver = Solver::new(boundary, &compiled, &weights, TieBreak::Weight).unwrap();
        assert_eq!(solver.pattern_spans(), checker_pattern_spans(&solver));
        let leucine = solver.seq.codons()[3];
        let ctg = solver
            .seq
            .synonyms_at(3)
            .iter()
            .copied()
            .find(|codon| codon.to_bytes() == *b"CTG")
            .unwrap();
        solver.set_codon(3, ctg).unwrap();
        let spans = solver.pattern_spans();
        assert_eq!(spans, checker_pattern_spans(&solver));
        assert!(
            spans
                .iter()
                .any(|(start, end, id)| { *start == 4 && *end == 10 && id == "c_ecori" }),
            "the boundary EcoRI site was not recorded: {spans:?}"
        );
        solver.set_codon(3, leucine).unwrap();
        assert_eq!(solver.pattern_spans(), checker_pattern_spans(&solver));

        // AAA AAG AAA. AAG -> AAA joins the two lysine runs across the codon.
        let merge = CodonSequence::from_cds("ATGAAAAAGAAATAA", 11).unwrap();
        let mut solver = Solver::new(merge, &compiled, &weights, TieBreak::Weight).unwrap();
        assert_eq!(solver.pattern_spans(), checker_pattern_spans(&solver));
        let lysine = solver
            .seq
            .synonyms_at(2)
            .iter()
            .copied()
            .find(|codon| codon.to_bytes() == *b"AAA")
            .unwrap();
        solver.set_codon(2, lysine).unwrap();
        let spans = solver.pattern_spans();
        assert_eq!(spans, checker_pattern_spans(&solver), "merged run");
        assert!(
            spans
                .iter()
                .any(|(start, end, id)| { *start == 3 && *end == 12 && id == "c_homo" }),
            "merged homopolymer span drifted: {spans:?}"
        );

        let cds = format!(
            "ATG{}{}{}TAA",
            "GCT".repeat(5),
            "CGT".repeat(4),
            "GGA".repeat(3)
        );
        let mut sequence = CodonSequence::from_cds(&cds, 11)
            .unwrap()
            .with_five_prime_flank("GAATT")
            .unwrap()
            .with_three_prime_flank("AAAAAA")
            .unwrap();
        sequence.add_lock(4, 5, LockReason::UserLocked).unwrap();
        let mut solver = Solver::new(sequence, &compiled, &weights, TieBreak::Weight).unwrap();
        assert_eq!(
            solver.pattern_spans(),
            checker_pattern_spans(&solver),
            "initial"
        );
        let mut state = 0x0F1E_2D3C_4B5A_6978u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let positions = solver.seq.codon_count();
        for step in 0..360 {
            let index = (next() as usize) % positions;
            if solver.seq.is_locked(index) {
                continue;
            }
            let synonyms = solver.seq.synonyms_at(index);
            let codon = synonyms[(next() as usize) % synonyms.len()];
            solver.set_codon(index, codon).unwrap();
            assert_eq!(
                solver.pattern_spans(),
                checker_pattern_spans(&solver),
                "step {step} codon {index}"
            );
        }
    }

    /// The per-codon explanation used to call `ConstraintChecker::check` once
    /// per deviant codon. The shared tally has to name the same ids, including
    /// a revert that changes a k-mer key from legal to over-limit. A counter
    /// stuck at zero drops `c_kmer` and this comparison fails.
    #[test]
    fn explain_deviations_matches_a_full_check_per_codon() {
        use crate::spec::{HomopolymerBase, Strand};
        use cf_types::LockReason;

        let compiled = checker(vec![
            (
                "c_gc",
                Rule::GcWindow {
                    window_nt: 6,
                    min: 0.0,
                    max: 0.3,
                },
            ),
            (
                "c_kmer",
                Rule::KmerRepeat {
                    k: 6,
                    max_occurrences: 1,
                    count_overlapping: true,
                    canonical_strand: false,
                },
            ),
            (
                "c_ecori",
                Rule::RestrictionSite {
                    enzyme: "EcoRI".into(),
                    site: "GAATTC".into(),
                    strand: Strand::Both,
                    max_occurrences: 0,
                },
            ),
            (
                "c_homo",
                Rule::Homopolymer {
                    base: HomopolymerBase::Any,
                    max_run: 4,
                },
            ),
            (
                "c_regex",
                Rule::CustomRegex {
                    pattern: "CCA".into(),
                    max_occurrences: 0,
                    flags: vec![],
                },
            ),
        ]);
        // Codon 1 TTA -> CTG duplicates CTGAAA (codon 3 + codon 4) and pushes
        // a 6-nt window over 0.3 GC. CTG is weighted above the other leucines
        // so that revert is the one the k-mer counter has to see; equal
        // weights would pick CTA and the repeat would never form. Codon 5
        // GAG -> GAA joins a run of five A's. Codon 8 TTT -> TTC completes
        // GAATTC. Codon 9 CCC -> CCA matches the regex. Codon 10 is a real
        // deviation and is locked, so neither side may report it.
        let cds = [
            "ATG", "TTA", "AAA", "CTG", "AAA", "GAG", "AAA", "GAA", "TTT", "CCC", "GCT", "TAA",
        ]
        .join("");
        let mut sequence = CodonSequence::from_cds(&cds, 11)
            .unwrap()
            .with_five_prime_flank("ACGTACGT")
            .unwrap()
            .with_three_prime_flank("TTAA")
            .unwrap();
        sequence.add_lock(10, 11, LockReason::UserLocked).unwrap();
        let mut weights = weights_one();
        weights[Codon::parse("CTG").unwrap().index_usize()] = 2.0;

        let mut solver =
            Solver::new(sequence.clone(), &compiled, &weights, TieBreak::Weight).unwrap();
        assert_eq!(
            solver.actionable_tally(),
            actionable_counts(&solver.seq, &compiled),
            "initial tally"
        );
        for index in 0..solver.seq.codon_count() {
            if solver.seq.is_locked(index) {
                continue;
            }
            let original = solver.seq.codons()[index];
            let synonyms = solver.seq.synonyms_at(index).to_vec();
            for synonym in synonyms {
                solver.set_codon(index, synonym).unwrap();
                assert_eq!(
                    solver.actionable_tally(),
                    actionable_counts(&solver.seq, &compiled),
                    "tally after codon {index} -> {synonym}"
                );
            }
            solver.set_codon(index, original).unwrap();
            assert_eq!(
                solver.actionable_tally(),
                actionable_counts(&solver.seq, &compiled),
                "tally restored at {index}"
            );
        }

        let found = explain_deviations(&sequence, &compiled, &weights).unwrap();
        let expected = deviations_from_full_checks(&sequence, &compiled, &weights);
        assert_eq!(found, expected);
        let indexes: Vec<usize> = found.iter().map(|item| item.codon_index).collect();
        assert!(
            !indexes
                .iter()
                .any(|index| *index == 0 || *index == 10 || *index == 11),
            "a locked codon was explained: {indexes:?}"
        );
        let mut covered = std::collections::BTreeSet::new();
        for item in &found {
            covered.extend(item.constraints_if_reverted.iter().cloned());
        }
        for id in ["c_kmer", "c_gc", "c_homo", "c_ecori", "c_regex"] {
            assert!(
                covered.contains(id),
                "{id} never rose on a single revert: {found:?}"
            );
        }
    }

    /// Blocking actionable counts from a full `ConstraintChecker::check`.
    ///
    /// The oracle for the shared tally. It stays in the test module so a
    /// release build cannot call it from the explanation path.
    fn actionable_counts(
        sequence: &CodonSequence,
        checker: &ConstraintChecker,
    ) -> BTreeMap<String, usize> {
        let mut nt = sequence.five_prime_flank().to_vec();
        for codon in sequence.codons() {
            nt.extend_from_slice(&codon.to_bytes());
        }
        nt.extend_from_slice(sequence.three_prime_flank());
        let start = sequence.five_prime_flank().len();
        let editable = editable_spans_of(sequence);
        let view = ConstructView::new(&nt, vec![(start, start + sequence.codon_count() * 3)])
            .with_editable_spans(editable);
        let mut counts = BTreeMap::new();
        for violation in checker.check(&view) {
            if violation.class == ConstraintClass::Blocking && violation.actionable {
                *counts.entry(violation.constraint_id).or_insert(0) += 1;
            }
        }
        counts
    }

    /// A GC-changing edit shifts every later prefix entry by the same integer,
    /// and a GC-neutral edit that moves GC inside the codon does not. Either
    /// mistake leaves the value [`Solver::prefix_at`] reads disagreeing with
    /// a rebuild. The stored array lags until `flush_prefix`; that is the
    /// trial path, and the flushed array has to agree too.
    fn assert_logical_prefix(solver: &Solver<'_>, label: &str) {
        let rebuilt = gc_prefix(&solver.nt);
        assert_eq!(rebuilt.len(), solver.prefix.len(), "{label}");
        for (index, count) in rebuilt.iter().enumerate() {
            assert_eq!(solver.prefix_at(index), *count, "{label} at {index}");
        }
    }

    #[test]
    fn gc_prefix_tail_matches_a_full_rebuild() {
        let compiled = checker(vec![(
            "c_gc",
            Rule::GcWindow {
                window_nt: 12,
                min: 0.0,
                max: 1.0,
            },
        )]);
        let cds = format!("ATG{}TAA", "CTG".repeat(40));
        let sequence = CodonSequence::from_cds(&cds, 11).unwrap();
        let weights = weights_one();
        let mut solver = Solver::new(sequence, &compiled, &weights, TieBreak::Weight).unwrap();
        assert_logical_prefix(&solver, "initial");
        let index = 20;
        let find = |solver: &Solver<'_>, bases: &[u8; 3]| {
            solver
                .seq
                .synonyms_at(index)
                .iter()
                .copied()
                .find(|codon| codon.to_bytes() == *bases)
                .unwrap()
        };
        // CTG (GC 2) -> TTA (GC 0). The tail has to move.
        let tta = find(&solver, b"TTA");
        solver.set_codon(index, tta).unwrap();
        assert_logical_prefix(&solver, "tail after CTG -> TTA");
        solver.flush_prefix();
        assert_eq!(
            solver.prefix,
            gc_prefix(&solver.nt),
            "flushed tail after CTG -> TTA"
        );
        // TTA (0,0,0) -> CTT (1,0,0) is not neutral. CTT (1,0,0) -> TTG (0,0,1)
        // keeps the codon GC count at 1 and moves it to the last base.
        let ctt = find(&solver, b"CTT");
        solver.set_codon(index, ctt).unwrap();
        let ttg = find(&solver, b"TTG");
        solver.set_codon(index, ttg).unwrap();
        assert_logical_prefix(&solver, "prefix after GC moved inside the codon");
        solver.flush_prefix();
        assert_eq!(
            solver.prefix,
            gc_prefix(&solver.nt),
            "flushed prefix after GC moved inside the codon"
        );
        let positions = solver.seq.codon_count();
        let mut state = 0xC0FF_EE00_1234_ABCD_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for step in 0..80 {
            let at = (next() as usize) % positions;
            if solver.seq.is_locked(at) {
                continue;
            }
            let synonyms = solver.seq.synonyms_at(at).to_vec();
            let codon = synonyms[(next() as usize) % synonyms.len()];
            solver.set_codon(at, codon).unwrap();
            assert_logical_prefix(&solver, &format!("step {step} codon {at}"));
        }
        solver.flush_prefix();
        assert_eq!(
            solver.prefix,
            gc_prefix(&solver.nt),
            "flushed prefix after the random edits"
        );
    }

    /// The explanation `explain_deviations` used to build: one full check, then
    /// another full check for each codon put back on its own.
    fn deviations_from_full_checks(
        sequence: &CodonSequence,
        checker: &ConstraintChecker,
        weights: &[f64; 64],
    ) -> Vec<BaselineDeviation> {
        let baseline = most_frequent(sequence, weights).unwrap();
        let mut probe = sequence.clone();
        let before = actionable_counts(&probe, checker);
        let mut out = Vec::new();
        for index in 0..probe.codon_count() {
            if probe.is_locked(index) {
                continue;
            }
            let codon = probe.codons()[index];
            let preferred = baseline.codons()[index];
            if codon == preferred {
                continue;
            }
            probe.substitute(index, preferred).unwrap();
            let after = actionable_counts(&probe, checker);
            probe.substitute(index, codon).unwrap();
            let mut constraints_if_reverted: Vec<String> = after
                .iter()
                .filter(|(id, count)| **count > *before.get(*id).unwrap_or(&0))
                .map(|(id, _)| (*id).clone())
                .collect();
            constraints_if_reverted.sort();
            out.push(BaselineDeviation {
                codon_index: index,
                codon: codon.to_string_upper(),
                weight: weights[codon.index_usize()],
                baseline: preferred.to_string_upper(),
                baseline_weight: weights[preferred.index_usize()],
                constraints_if_reverted,
            });
        }
        out
    }
}
