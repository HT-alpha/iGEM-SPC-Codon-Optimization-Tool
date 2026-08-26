//! Constraint feasibility analysis: five checks, run once at submission.
//!
//! Ticket P1-19. Contract: `planning/03-algorithms/06-constraint-repair.md` §5.
//! Amendment implemented: `AMD-010` (jointly unsatisfiable constraint sets).
//!
//! WHY THIS EXISTS. The prototype in `codon_6.0.py` reports
//! `Exceeded max restarts (1000) without finding a solution` and stops. The
//! user learns nothing about why, or which constraint to relax. Every finding
//! here names a constraint and, where the check is exact, how far it is from
//! being satisfiable.
//!
//! WHY IT RUNS BEFORE REPAIR. See §5.1 of the contract. Repair's bounded budget
//! makes each call terminate; it does not stop repair from spending that budget
//! on every candidate forever when a constraint is unsatisfiable in principle.
//! The two halves together are what make repair's termination principled rather
//! than a restart cap with better manners.
//!
//! NON-ACTIONABLE VIOLATIONS LAND HERE. Check 1 is the destination that
//! `DEP-CONS-002` routes them to. That reasoning is not re-derived here —
//! `DEP-CONS-002` is the source of truth, and [`Violation::actionable`] is the
//! mechanism.

use std::collections::{BTreeMap, BTreeSet};

use cf_types::codon::Codon;
use cf_types::error::{code, Result};
use cf_types::seqrec::CodonSequence;

use crate::detect::{ConstraintChecker, ConstructView, KmerMasking};
use crate::motif::MotifAutomaton;
use crate::spec::{ConstraintClass, HomopolymerBase, Rule};

/// Which of the five checks produced a finding.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum CheckId {
    LockedSpan,
    CountingBound,
    PerPosition,
    AutomatonReachability,
    GcAttainability,
}

impl CheckId {
    pub const ALL: [CheckId; 5] = [
        CheckId::LockedSpan,
        CheckId::CountingBound,
        CheckId::PerPosition,
        CheckId::AutomatonReachability,
        CheckId::GcAttainability,
    ];

    pub fn as_str(self) -> &'static str {
        match self {
            CheckId::LockedSpan => "locked_span",
            CheckId::CountingBound => "counting_bound",
            CheckId::PerPosition => "per_position",
            CheckId::AutomatonReachability => "automaton_reachability",
            CheckId::GcAttainability => "gc_attainability",
        }
    }

    pub fn error_code(self) -> &'static str {
        match self {
            CheckId::LockedSpan => code::FEAS_LOCKED_SPAN,
            CheckId::CountingBound => code::FEAS_COUNTING_BOUND,
            CheckId::PerPosition => code::FEAS_PER_POSITION,
            CheckId::AutomatonReachability => code::FEAS_UNREACHABLE,
            CheckId::GcAttainability => code::FEAS_GC_UNATTAINABLE,
        }
    }

    /// Whether a finding from this check is a proof of unsatisfiability, or
    /// only evidence of it.
    ///
    /// Exact checks may be acted on automatically; necessary-condition checks
    /// are reported for a human. Conflating the two is how a feasibility gate
    /// starts refusing satisfiable jobs.
    pub fn is_exact(self) -> bool {
        match self {
            CheckId::CountingBound | CheckId::AutomatonReachability | CheckId::GcAttainability => {
                true
            }
            CheckId::LockedSpan | CheckId::PerPosition => false,
        }
    }
}

/// One problem, with the constraint to relax.
#[derive(Clone, Debug, PartialEq)]
pub struct Finding {
    pub check: CheckId,
    pub code: &'static str,
    pub constraint_id: Option<String>,
    pub span_nt: Option<(usize, usize)>,
    pub detail: String,
    /// What would have to change for this finding to go away. `None` only when
    /// the check genuinely cannot say.
    pub suggested_relaxation: Option<String>,
    /// True when this finding proves unsatisfiability rather than suggesting it.
    pub exact: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FeasibilityReport {
    /// No finding was produced. Not a guarantee of satisfiability: two of the
    /// five checks are necessary conditions only.
    pub feasible: bool,
    /// Always true, and stated in the report rather than assumed by the reader.
    pub necessary_conditions_only: bool,
    pub findings: Vec<Finding>,
    /// Every check that ran, whether or not it found anything, so a silent
    /// check is distinguishable from a passing one.
    pub checks_run: Vec<CheckId>,
}

impl FeasibilityReport {
    pub fn findings_for(&self, check: CheckId) -> impl Iterator<Item = &Finding> {
        self.findings.iter().filter(move |f| f.check == check)
    }

    /// Findings that prove unsatisfiability.
    pub fn proofs(&self) -> impl Iterator<Item = &Finding> {
        self.findings.iter().filter(|f| f.exact)
    }
}

/// What can vary: the construct, plus the coding regions that may be edited.
#[derive(Clone, Debug)]
pub struct DesignSpace<'a> {
    /// The assembled construct as it currently stands.
    pub construct: &'a [u8],
    /// `(nucleotide offset, sequence)` per editable coding region.
    pub regions: Vec<(usize, &'a CodonSequence)>,
}

impl<'a> DesignSpace<'a> {
    pub fn new(construct: &'a [u8], regions: Vec<(usize, &'a CodonSequence)>) -> Self {
        Self { construct, regions }
    }

    pub fn cds_spans(&self) -> Vec<(usize, usize)> {
        self.regions
            .iter()
            .map(|(offset, seq)| (*offset, offset + seq.codon_count() * 3))
            .collect()
    }

    fn view(&self, masking: KmerMasking) -> ConstructView<'a> {
        ConstructView::new(self.construct, self.cds_spans())
            .with_masking(masking)
            .with_editable_spans(crate::repair::editable_spans(
                &self
                    .regions
                    .iter()
                    .map(|(offset, sequence)| (*offset, (*sequence).clone()))
                    .collect::<Vec<_>>(),
            ))
    }

    /// Alternatives at a nucleotide offset, or `None` if it is not the start of
    /// an editable codon.
    fn alternatives_at(&self, offset: usize) -> Option<(&'a CodonSequence, usize, &'a [Codon])> {
        for (start, seq) in &self.regions {
            if offset < *start {
                continue;
            }
            let delta = offset - start;
            if delta % 3 != 0 {
                continue;
            }
            let index = delta / 3;
            if index < seq.codon_count() && !seq.is_locked(index) {
                return Some((seq, index, seq.synonyms_at(index)));
            }
        }
        None
    }
}

/// Run all five checks.
pub fn analyse(
    space: &DesignSpace,
    checker: &ConstraintChecker,
    masking: KmerMasking,
) -> Result<FeasibilityReport> {
    let mut findings = Vec::new();

    check_locked_span(space, checker, masking, &mut findings);
    check_counting_bound(space, checker, &mut findings);
    check_per_position(space, checker, &mut findings);
    check_automaton_reachability(space, checker, &mut findings)?;
    check_gc_attainability(space, checker, &mut findings);

    findings.sort_by(|a, b| {
        (a.check, a.span_nt, &a.constraint_id).cmp(&(b.check, b.span_nt, &b.constraint_id))
    });
    Ok(FeasibilityReport {
        feasible: findings.is_empty(),
        necessary_conditions_only: true,
        findings,
        checks_run: CheckId::ALL.to_vec(),
    })
}

// --------------------------------------------------------------------------- //
// 1. locked span — where DEP-CONS-002 routes non-actionable violations
// --------------------------------------------------------------------------- //

fn check_locked_span(
    space: &DesignSpace,
    checker: &ConstraintChecker,
    masking: KmerMasking,
    out: &mut Vec<Finding>,
) {
    let view = space.view(masking);
    for violation in checker.unfixable(&view) {
        out.push(Finding {
            check: CheckId::LockedSpan,
            code: CheckId::LockedSpan.error_code(),
            constraint_id: Some(violation.constraint_id.clone()),
            span_nt: Some((violation.start_nt, violation.end_nt)),
            detail: format!(
                "{} is violated entirely within immutable sequence: {}",
                violation.constraint_id, violation.detail
            ),
            suggested_relaxation: Some(format!(
                "no synonymous assignment can fix this; relax or disable {}, \
                 or change the immutable sequence itself",
                violation.constraint_id
            )),
            exact: CheckId::LockedSpan.is_exact(),
        });
    }
}

// --------------------------------------------------------------------------- //
// 2. counting bound — pigeonhole over k-mers
// --------------------------------------------------------------------------- //

fn check_counting_bound(space: &DesignSpace, checker: &ConstraintChecker, out: &mut Vec<Finding>) {
    let n = space.construct.len() as u128;
    for constraint in checker
        .set()
        .constraints
        .iter()
        .filter(|c| c.enabled && c.class == ConstraintClass::Blocking)
    {
        let Rule::KmerRepeat {
            k,
            max_occurrences,
            canonical_strand,
            ..
        } = &constraint.rule
        else {
            continue;
        };
        let k = *k as u128;
        if n < k {
            continue;
        }
        let windows = n - k + 1;
        // Distinct keys available. Canonical-strand folding roughly halves the
        // key space: a k-mer and its reverse complement share a bucket, except
        // for the palindromes, which are their own. Rounding up keeps this a
        // sound upper bound on capacity, so the check never fires spuriously.
        let alphabet: u128 = 1u128 << (2 * k.min(40));
        let distinct = if *canonical_strand {
            alphabet.div_ceil(2)
        } else {
            alphabet
        };
        let capacity = distinct.saturating_mul(*max_occurrences as u128);
        if capacity < windows {
            let min_k = smallest_feasible_k(n, *max_occurrences as u128, *canonical_strand);
            out.push(Finding {
                check: CheckId::CountingBound,
                code: CheckId::CountingBound.error_code(),
                constraint_id: Some(constraint.id.clone()),
                span_nt: None,
                detail: format!(
                    "{}: a {n}-nt construct has {windows} overlapping {k}-mer windows, but at \
                     most {max_occurrences} occurrence(s) of each of {distinct} distinct keys \
                     gives capacity {capacity}. No sequence of this length can satisfy it.",
                    constraint.id
                ),
                suggested_relaxation: Some(match min_k {
                    Some(min_k) => format!(
                        "raise k to at least {min_k}, or raise max_occurrences to at least {}",
                        windows.div_ceil(distinct)
                    ),
                    None => format!(
                        "raise max_occurrences to at least {}",
                        windows.div_ceil(distinct)
                    ),
                }),
                exact: CheckId::CountingBound.is_exact(),
            });
        }
    }
}

fn smallest_feasible_k(n: u128, max_occurrences: u128, canonical: bool) -> Option<u128> {
    for k in 4u128..=32 {
        if n < k {
            return None;
        }
        let alphabet = 1u128 << (2 * k);
        let distinct = if canonical {
            alphabet.div_ceil(2)
        } else {
            alphabet
        };
        // `windows` is `n - k + 1`; written this way because clippy objects to
        // the `>= x + 1` form and the intent is "capacity covers the windows".
        let windows = n - k + 1;
        if distinct.saturating_mul(max_occurrences) >= windows {
            return Some(k);
        }
    }
    None
}

// --------------------------------------------------------------------------- //
// 3. per-position exhaustion
// --------------------------------------------------------------------------- //

/// A position is exhausted when every synonym produces a blocking violation
/// whose span touches only immutable positions plus this one.
///
/// The restriction to such violations is what makes this sound: if a violation
/// also covers another editable position, that position might change too, so
/// the finding would be a guess. Confined violations cannot be escaped by any
/// other choice, which is why this reports a real problem rather than a
/// suspicion — though it remains a necessary condition, since it examines one
/// position at a time.
fn check_per_position(space: &DesignSpace, checker: &ConstraintChecker, out: &mut Vec<Finding>) {
    let cds_spans = space.cds_spans();
    for (offset, seq) in &space.regions {
        for index in seq.editable_positions() {
            let alternatives = seq.synonyms_at(index);
            if alternatives.len() <= 1 {
                continue;
            }
            let position_nt = offset + index * 3;
            let mut blocked_by: BTreeSet<String> = BTreeSet::new();
            let mut all_blocked = true;

            for alternative in alternatives {
                let mut construct = space.construct.to_vec();
                construct[position_nt..position_nt + 3].copy_from_slice(&alternative.to_bytes());
                let view = ConstructView::new(&construct, cds_spans.clone());
                let confined: Vec<_> = checker
                    .check(&view)
                    .into_iter()
                    .filter(|v| {
                        v.class == ConstraintClass::Blocking
                            && v.start_nt < position_nt + 3
                            && position_nt < v.end_nt
                            && confined_to(v.start_nt, v.end_nt, position_nt, &cds_spans)
                    })
                    .collect();
                if confined.is_empty() {
                    all_blocked = false;
                    break;
                }
                blocked_by.extend(confined.into_iter().map(|v| v.constraint_id));
            }

            if all_blocked {
                out.push(Finding {
                    check: CheckId::PerPosition,
                    code: CheckId::PerPosition.error_code(),
                    constraint_id: blocked_by.iter().next().cloned(),
                    span_nt: Some((position_nt, position_nt + 3)),
                    detail: format!(
                        "codon {index} of the region at {offset} has {} synonym(s), and every \
                         one violates {:?} within immutable context",
                        alternatives.len(),
                        blocked_by
                    ),
                    suggested_relaxation: Some(format!(
                        "relax {:?}, or unlock surrounding context so another position can move",
                        blocked_by
                    )),
                    exact: CheckId::PerPosition.is_exact(),
                });
            }
        }
    }
}

/// Whether `[start, end)` touches no editable position other than `position_nt`.
fn confined_to(start: usize, end: usize, position_nt: usize, cds_spans: &[(usize, usize)]) -> bool {
    for nt in start..end {
        if nt >= position_nt && nt < position_nt + 3 {
            continue;
        }
        if cds_spans.iter().any(|&(lo, hi)| nt >= lo && nt < hi) {
            return false;
        }
    }
    true
}

// --------------------------------------------------------------------------- //
// 4. automaton reachability — exact (DEP-CONS-003)
// --------------------------------------------------------------------------- //

/// Exact forward reachability over `(codon position) x (motif DFA state)`.
///
/// The originally-specified beam DP is unsound here: an empty beam proves the
/// heuristic lost the surviving paths, not that none exist, and a feasibility
/// gate that refuses satisfiable jobs is worse than one that misses some
/// infeasible ones. The exact form is `O(L · |Q| · S)`, needs no beam, and
/// needs nothing from `cf-search`.
///
/// A `CF-FEAS-004` finding is therefore a proof: no synonymous assignment of
/// the editable region avoids the regular constraints.
fn check_automaton_reachability(
    space: &DesignSpace,
    checker: &ConstraintChecker,
    out: &mut Vec<Finding>,
) -> Result<()> {
    let Some(automaton) = regular_automaton(checker)? else {
        return Ok(());
    };

    for (offset, seq) in &space.regions {
        // Feed the immutable prefix to fix the entry state. A motif already
        // present in the prefix is check 1's business, not this one, so an
        // accepting state here is not treated specially; only the resulting
        // state matters. A non-base (a masked position) resets, since no match
        // can span it.
        let mut state = automaton.start();
        for &byte in &space.construct[..*offset] {
            state = automaton
                .step(state, byte)
                .unwrap_or_else(|| automaton.start());
        }

        let mut reachable: BTreeSet<u32> = BTreeSet::from([state]);
        for index in 0..seq.codon_count() {
            let alternatives: Vec<Codon> = if seq.is_locked(index) {
                vec![seq.codons()[index]]
            } else {
                seq.synonyms_at(index).to_vec()
            };
            let mut next: BTreeSet<u32> = BTreeSet::new();
            for &from in &reachable {
                for alternative in &alternatives {
                    if let Some(to) = automaton.step_codon(from, alternative.to_bytes()) {
                        next.insert(to);
                    }
                }
            }
            if next.is_empty() {
                let closed = automaton.pattern_literals();
                out.push(Finding {
                    check: CheckId::AutomatonReachability,
                    code: CheckId::AutomatonReachability.error_code(),
                    constraint_id: None,
                    span_nt: Some((offset + index * 3, offset + index * 3 + 3)),
                    detail: format!(
                        "no assignment of synonymous codons reaches codon {index} of the region \
                         at {offset} without matching one of {closed:?}. Exhaustive over all \
                         {} synonym(s) at that position and every reachable automaton state.",
                        alternatives.len()
                    ),
                    suggested_relaxation: Some(
                        "relax one of the listed motif or restriction-site constraints, or \
                         unlock a codon in the closing window"
                            .into(),
                    ),
                    exact: CheckId::AutomatonReachability.is_exact(),
                });
                break;
            }
            reachable = next;
        }
    }
    Ok(())
}

/// Number of reachable automaton states at each codon of each region.
///
/// Exposed so that the *exactness* of check 4 is observable rather than merely
/// asserted. A truncating implementation — the beam DP `DEP-CONS-003` replaced
/// — would report 1 everywhere; exact reachability reports the true branching.
/// `TST-CONS-006`'s exactness test reads this, because a scenario-based test
/// cannot reliably distinguish the two: whether a width-1 beam survives depends
/// on which state it happens to keep, so such a test passes or fails by luck.
pub fn reachable_state_counts(
    space: &DesignSpace,
    checker: &ConstraintChecker,
) -> Result<Vec<Vec<usize>>> {
    let Some(automaton) = regular_automaton(checker)? else {
        return Ok(vec![Vec::new(); space.regions.len()]);
    };
    let mut out = Vec::with_capacity(space.regions.len());
    for (offset, seq) in &space.regions {
        let mut state = automaton.start();
        for &byte in &space.construct[..*offset] {
            state = automaton
                .step(state, byte)
                .unwrap_or_else(|| automaton.start());
        }
        let mut reachable: BTreeSet<u32> = BTreeSet::from([state]);
        let mut counts = Vec::with_capacity(seq.codon_count());
        for index in 0..seq.codon_count() {
            let alternatives: Vec<Codon> = if seq.is_locked(index) {
                vec![seq.codons()[index]]
            } else {
                seq.synonyms_at(index).to_vec()
            };
            let mut next: BTreeSet<u32> = BTreeSet::new();
            for &from in &reachable {
                for alternative in &alternatives {
                    if let Some(to) = automaton.step_codon(from, alternative.to_bytes()) {
                        next.insert(to);
                    }
                }
            }
            counts.push(next.len());
            if next.is_empty() {
                break;
            }
            reachable = next;
        }
        out.push(counts);
    }
    Ok(out)
}

/// Compile the regular (finite-automaton-expressible) blocking constraints.
///
/// Homopolymer limits are included: `max_run = n` is the regular language "no
/// run of n+1", i.e. four literal motifs. Excluding them would make the check
/// weaker for no reason.
fn regular_automaton(checker: &ConstraintChecker) -> Result<Option<ReachabilityAutomaton>> {
    let mut specs: Vec<(usize, String, bool)> = Vec::new();
    for (index, constraint) in checker.set().constraints.iter().enumerate() {
        if !constraint.enabled || constraint.class != ConstraintClass::Blocking {
            continue;
        }
        match &constraint.rule {
            // `max_occurrences > 0` is a counting constraint, not a regular
            // one, so it is out of scope for this check.
            Rule::ForbiddenMotif {
                motif,
                strand,
                max_occurrences: 0,
            } => specs.push((index, motif.clone(), *strand == crate::spec::Strand::Both)),
            Rule::RestrictionSite {
                site,
                strand,
                max_occurrences: 0,
                ..
            } => specs.push((index, site.clone(), *strand == crate::spec::Strand::Both)),
            Rule::Homopolymer { base, max_run } => {
                let bases: &[&str] = match base {
                    HomopolymerBase::Any => &["A", "C", "G", "T"],
                    HomopolymerBase::A => &["A"],
                    HomopolymerBase::C => &["C"],
                    HomopolymerBase::G => &["G"],
                    HomopolymerBase::T => &["T"],
                };
                for b in bases {
                    specs.push((index, b.repeat(*max_run as usize + 1), false));
                }
            }
            _ => {}
        }
    }
    if specs.is_empty() {
        return Ok(None);
    }
    Ok(Some(ReachabilityAutomaton::new(MotifAutomaton::build(
        &specs,
    )?)))
}

/// Thin wrapper exposing the automaton's transition function for reachability.
struct ReachabilityAutomaton {
    inner: MotifAutomaton,
}

impl ReachabilityAutomaton {
    fn new(inner: MotifAutomaton) -> Self {
        Self { inner }
    }

    fn start(&self) -> u32 {
        self.inner.start_state() as u32
    }

    fn step(&self, state: u32, byte: u8) -> Option<u32> {
        self.inner
            .transition(state as usize, byte)
            .map(|s| s as u32)
    }

    fn is_accepting(&self, state: u32) -> bool {
        self.inner.is_accepting(state as usize)
    }

    /// Feed three bases, rejecting if any intermediate state accepts.
    fn step_codon(&self, state: u32, codon: [u8; 3]) -> Option<u32> {
        let mut current = state;
        for byte in codon {
            current = self.step(current, byte)?;
            if self.is_accepting(current) {
                return None;
            }
        }
        Some(current)
    }

    fn pattern_literals(&self) -> Vec<String> {
        let mut literals: Vec<String> = self
            .inner
            .patterns()
            .iter()
            .map(|p| p.literal.clone())
            .collect();
        literals.sort();
        literals.dedup();
        literals
    }
}

// --------------------------------------------------------------------------- //
// 5. GC-window attainability
// --------------------------------------------------------------------------- //

/// For each window, the minimum and maximum GC fraction any synonymous
/// assignment can reach. Exact for a window in isolation: within one window
/// each codon appears once, so the extremes are attained by choosing each
/// codon's extreme independently.
fn check_gc_attainability(
    space: &DesignSpace,
    checker: &ConstraintChecker,
    out: &mut Vec<Finding>,
) {
    let n = space.construct.len();
    for constraint in checker
        .set()
        .constraints
        .iter()
        .filter(|c| c.enabled && c.class == ConstraintClass::Blocking)
    {
        let Rule::GcWindow {
            window_nt,
            min,
            max,
        } = &constraint.rule
        else {
            continue;
        };
        let window = *window_nt as usize;
        if window == 0 || n < window {
            continue;
        }
        for start in 0..=n - window {
            let end = start + window;
            let (lo, hi, denominator) = attainable_gc(space, start, end);
            if denominator == 0 {
                continue;
            }
            let lo_fraction = lo as f64 / denominator as f64;
            let hi_fraction = hi as f64 / denominator as f64;
            if hi_fraction < *min || lo_fraction > *max {
                out.push(Finding {
                    check: CheckId::GcAttainability,
                    code: CheckId::GcAttainability.error_code(),
                    constraint_id: Some(constraint.id.clone()),
                    span_nt: Some((start, end)),
                    detail: format!(
                        "{}: window {start}..{end} can reach GC in [{lo_fraction:.3}, \
                         {hi_fraction:.3}], which does not intersect the required \
                         [{min:.3}, {max:.3}]",
                        constraint.id
                    ),
                    suggested_relaxation: Some(format!(
                        "widen the band to include [{lo_fraction:.3}, {hi_fraction:.3}], \
                         or unlock codons overlapping {start}..{end}"
                    )),
                    exact: CheckId::GcAttainability.is_exact(),
                });
                // One finding per constraint is enough to act on; reporting
                // every window of a long unsatisfiable stretch is noise.
                break;
            }
        }
    }
}

/// `(min GC count, max GC count, denominator)` over a window.
fn attainable_gc(space: &DesignSpace, start: usize, end: usize) -> (usize, usize, usize) {
    let mut min = 0usize;
    let mut max = 0usize;
    let mut denominator = 0usize;
    let mut nt = start;
    while nt < end {
        if let Some((_seq, _index, alternatives)) = space.alternatives_at(nt) {
            // A whole codon starts here; count the part inside the window.
            let overlap = (end - nt).min(3);
            let mut lo = usize::MAX;
            let mut hi = 0usize;
            for alternative in alternatives {
                let bytes = alternative.to_bytes();
                let gc = bytes[..overlap]
                    .iter()
                    .filter(|b| matches!(b, b'G' | b'C'))
                    .count();
                lo = lo.min(gc);
                hi = hi.max(gc);
            }
            min += lo;
            max += hi;
            denominator += overlap;
            nt += overlap;
        } else {
            let byte = space.construct[nt];
            if matches!(byte.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T') {
                let gc = usize::from(matches!(byte.to_ascii_uppercase(), b'G' | b'C'));
                min += gc;
                max += gc;
                denominator += 1;
            }
            nt += 1;
        }
    }
    (min, max, denominator)
}

/// Findings grouped by check, for reporting.
pub fn summarize(report: &FeasibilityReport) -> BTreeMap<&'static str, usize> {
    let mut counts = BTreeMap::new();
    for check in CheckId::ALL {
        counts.insert(check.as_str(), report.findings_for(check).count());
    }
    counts
}
