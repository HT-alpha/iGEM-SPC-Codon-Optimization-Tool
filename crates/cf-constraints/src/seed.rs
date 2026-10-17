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
use std::collections::{BTreeMap, HashMap};
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

const MAX_CYCLES: u32 = 8;
const GC_SWEEPS: u32 = 8;
const PATTERN_ROUNDS: u32 = 4;
/// Pair lifts after a feasible seed. One codon at a time cannot trade a
/// higher-weight synonym for the neighbour that clears the repeat it creates.
/// `DEP-SRCH-018`.
const ASCENT_PASSES: u32 = 2;
const ASCENT_POSITIONS: usize = 24;
const ASCENT_PARTNERS: usize = 6;
const ASCENT_NEIGHBOUR: usize = 4;
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
/// Tries the most-frequent-codon assignment first, then the input. A feasible
/// result outranks an infeasible one. Two feasible results keep the higher
/// CAI on `weights` (stops and single-codon families excluded; a counted
/// weight of 0 scores 0).
pub fn resolve_feasible_seed(
    sequence: &CodonSequence,
    checker: &ConstraintChecker,
    weights: &[f64; 64],
) -> Result<SeedReport> {
    arm_seed_profile();
    let unsearched = unsearched_rules(checker);
    let mut best: Option<Attempt> = None;
    let mut solved = false;
    for tie in [TieBreak::Weight, TieBreak::LowIndex] {
        if solved {
            break;
        }
        for start in [most_frequent(sequence, weights)?, sequence.clone()] {
            let mut solver = Solver::new(start, checker, weights, tie)?;
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
            if best.as_ref().is_some_and(|item| item.feasible) {
                solved = true;
                break;
            }
        }
    }
    let attempt = best.expect("two starts were attempted");
    let solver = Solver::new(attempt.seq, checker, weights, TieBreak::Weight)?;
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
pub fn explain_deviations(
    sequence: &CodonSequence,
    checker: &ConstraintChecker,
    weights: &[f64; 64],
) -> Result<Vec<BaselineDeviation>> {
    let baseline = most_frequent(sequence, weights)?;
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
        probe.substitute(index, preferred)?;
        let after = actionable_counts(&probe, checker);
        probe.substitute(index, codon)?;
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
    window: usize,
    min: f64,
    max: f64,
}

struct KmerCfg {
    k: usize,
    max_occurrences: u32,
    canonical: bool,
}

struct KmerBook {
    rules: Vec<KmerCfg>,
    counts: Vec<HashMap<u64, u32>>,
    excess: Vec<u64>,
}

impl KmerBook {
    fn build(nt: &[u8], rules: Vec<KmerCfg>) -> Self {
        let mut book = Self {
            counts: vec![HashMap::new(); rules.len()],
            excess: vec![0; rules.len()],
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
                    book.bump(rule_index, key, 1);
                }
            }
        }
        book
    }

    fn excess(&self) -> u64 {
        self.excess.iter().sum()
    }

    fn bump(&mut self, rule: usize, key: u64, delta: i32) {
        let max = self.rules[rule].max_occurrences;
        let before = self.counts[rule].get(&key).copied().unwrap_or(0);
        let before_over = u64::from(before.saturating_sub(max));
        let next = before as i32 + delta;
        if next <= 0 {
            self.counts[rule].remove(&key);
            self.excess[rule] -= before_over;
            return;
        }
        let after = next as u32;
        self.counts[rule].insert(key, after);
        let after_over = u64::from(after.saturating_sub(max));
        self.excess[rule] = self.excess[rule] - before_over + after_over;
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
                    self.bump(rule_index, key, delta);
                }
            }
        }
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

struct Solver<'a> {
    seq: CodonSequence,
    checker: &'a ConstraintChecker,
    weights: &'a [f64; 64],
    nt: Vec<u8>,
    cds_off: usize,
    prefix: Vec<u32>,
    gc_penalty: f64,
    bands: Vec<GcBand>,
    kmers: KmerBook,
    editable: Vec<(usize, usize)>,
    homo_rules: Vec<(Option<u8>, usize)>,
    homo_penalty: u64,
    edits: u32,
    tie: TieBreak,
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
        let bands = gc_bands(checker);
        let gc_penalty = gc_penalty_of(&nt, &prefix, &bands, &editable_spans_of(&seq));
        let kmers = KmerBook::build(&nt, kmer_rules(checker));
        let editable = editable_spans_of(&seq);
        let homo_rules = homo_rules(checker);
        let homo_penalty = homo_total(&nt, 0, None, &homo_rules, &editable);
        Ok(Self {
            seq,
            checker,
            weights,
            nt,
            cds_off,
            prefix,
            gc_penalty,
            bands,
            kmers,
            editable,
            homo_rules,
            homo_penalty,
            edits: 0,
            tie,
        })
    }

    fn prefer(&self, candidate: Codon, incumbent: Codon) -> bool {
        match self.tie {
            TieBreak::Weight => prefer_weight(self.weights, candidate, incumbent),
            TieBreak::LowIndex => candidate.index() < incumbent.index(),
        }
    }

    fn run(&mut self) -> Result<u32> {
        let mut cycles = 0u32;
        if self.is_feasible() {
            self.polish()?;
            self.ascend_cai()?;
            return Ok(0);
        }
        for _ in 0..MAX_CYCLES {
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
                break;
            }
            if self.score() >= before {
                break;
            }
        }
        Ok(cycles)
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
        for index in hot {
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
            *self = Self::new(saved, self.checker, self.weights, self.tie)?;
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

    fn set_codon(&mut self, index: usize, codon: Codon) -> Result<Codon> {
        let old = self.seq.substitute(index, codon)?;
        if old == codon {
            return Ok(old);
        }
        let at = self.cds_off + index * 3;
        let bases = codon.to_bytes();
        self.homo_penalty = self.homo_after(at, bases);
        self.gc_penalty += self.gc_delta(at, bases);
        let kmer_started = seed_profile_on().then(Instant::now);
        self.kmers.replace(&mut self.nt, at, bases);
        if let Some(kmer_started) = kmer_started {
            SEED_PROFILE.kmer_replace.record(kmer_started);
        }
        let prefix_started = seed_profile_on().then(Instant::now);
        for pos in at..self.nt.len() {
            self.prefix[pos + 1] = self.prefix[pos] + gc_bit(self.nt[pos]);
        }
        if let Some(prefix_started) = prefix_started {
            SEED_PROFILE.prefix.record(prefix_started);
        }
        Ok(old)
    }

    fn commit(&mut self, index: usize, codon: Codon) -> Result<()> {
        let old = self.set_codon(index, codon)?;
        if old != codon {
            self.edits += 1;
        }
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
                let before = self.prefix[end] - self.prefix[start];
                let after = adjusted_gc(before, &self.nt, at, bases, start, end);
                delta += band_penalty(after, band) - band_penalty(before, band);
            }
        }
        delta
    }

    fn recompute_gc(&mut self) {
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

    fn is_feasible(&self) -> bool {
        let started = seed_profile_on().then(Instant::now);
        let view = self.view();
        let check_started = seed_profile_on().then(Instant::now);
        let feasible = self.checker.is_feasible(&view);
        if let Some(check_started) = check_started {
            SEED_PROFILE.checker_check.record(check_started);
        }
        if let Some(started) = started {
            SEED_PROFILE.is_feasible.record(started);
        }
        feasible
    }

    fn pattern_count(&self) -> usize {
        let started = seed_profile_on().then(Instant::now);
        let count = self.pattern_hits().len();
        if let Some(started) = started {
            SEED_PROFILE.pattern_count.record(started);
        }
        count
    }

    fn pattern_violations(&self) -> Vec<(usize, usize, usize)> {
        let started = seed_profile_on().then(Instant::now);
        let hits = self.pattern_hits();
        let count = hits.len();
        let found = hits
            .into_iter()
            .map(|violation| (violation.start_nt, violation.end_nt, count))
            .collect();
        if let Some(started) = started {
            SEED_PROFILE.pattern_violations.record(started);
        }
        found
    }

    /// Blocking actionable pattern hits.
    ///
    /// `check_patterns` omits k-mer and GC scans. Those rules are not patterns,
    /// so the list matches `checker.check` filtered by `is_pattern`.
    fn pattern_hits(&self) -> Vec<crate::detect::Violation> {
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
                window: *window_nt as usize,
                min: *min,
                max: *max,
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
}
