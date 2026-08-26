//! The constraint checker: compile a set, detect violations in a construct.
//!
//! Tickets P1-14, P1-15, P1-16. Contract:
//! `planning/02-schemas/constraint_set.schema.json`.
//!
//! Blocking and advisory violations are both reported, separately. A rule the
//! Phase 1 checker cannot evaluate (local structure, LM-dependent codon
//! frequency, custom regex) is listed by [`ConstraintChecker::undetectable`]
//! rather than being treated as satisfied — a constraint that is silently not
//! checked is worse than one that fails loudly.

use cf_types::error::Result;

use crate::kmer::{unpack, KmerIndex};
use crate::motif::MotifAutomaton;
use crate::spec::{
    Constraint, ConstraintClass, ConstraintSet, HomopolymerBase, Rule, Scope, Strand,
};
use crate::window::{find_runs, GcIndex};

/// How k-mer occurrences are counted across region boundaries.
///
/// This is a property of the *view*, not of the constraint: it describes what
/// the checker is looking at rather than what the user asked for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum KmerMasking {
    /// Count over the assembled construct. A repeat straddling the boundary
    /// between a coding region and the backbone is counted, which is correct:
    /// the synthesised molecule contains it.
    #[default]
    Assembled,
    /// Count coding regions and backbone as independent segments, so a k-mer
    /// spanning a boundary is counted by neither.
    ///
    /// This reproduces `codon_6.0.py`, which masks coding regions out of the
    /// backbone with `N` and then adds a separate within-CDS count. It exists
    /// so the legacy prototype's decisions can be reproduced exactly for the
    /// differential test in `P1-27`; see `docs/legacy-parity.md` for why it is
    /// not the default.
    SegmentedByRegion,
}

/// What the checker is looking at.
#[derive(Clone, Debug)]
pub struct ConstructView<'a> {
    /// The assembled construct, uppercase. Positions outside `ACGT` are treated
    /// as masked and break k-mers, motifs and runs.
    pub construct: &'a [u8],
    /// Half-open nucleotide spans of the coding regions.
    pub cds_spans: Vec<(usize, usize)>,
    pub kmer_masking: KmerMasking,
    /// Position of `construct[0]` within the enclosing molecule.
    ///
    /// Non-zero only when the view is a window onto a larger construct, which
    /// [`ConstraintChecker::check_region`] uses to bound its work. Every
    /// reported coordinate — in the span *and* in the human-readable detail —
    /// is offset by it, so a windowed check and a whole-construct check produce
    /// violations that compare equal. Offsetting the span alone would leave a
    /// detail string quoting window-local coordinates, which is the kind of
    /// discrepancy that survives review because both numbers look plausible.
    pub origin_nt: usize,
}

impl<'a> ConstructView<'a> {
    pub fn new(construct: &'a [u8], cds_spans: Vec<(usize, usize)>) -> Self {
        Self {
            construct,
            cds_spans,
            kmer_masking: KmerMasking::default(),
            origin_nt: 0,
        }
    }

    pub fn with_masking(mut self, masking: KmerMasking) -> Self {
        self.kmer_masking = masking;
        self
    }

    pub fn with_origin(mut self, origin_nt: usize) -> Self {
        self.origin_nt = origin_nt;
        self
    }

    /// Whether a nucleotide span touches any editable region.
    fn touches_editable(&self, start: usize, end: usize) -> bool {
        self.cds_spans
            .iter()
            .any(|&(lo, hi)| start < hi && lo < end)
    }

    /// Segments a k-mer rule counts over, given the masking policy.
    fn kmer_segments(&self) -> Vec<(usize, Vec<u8>)> {
        match self.kmer_masking {
            KmerMasking::Assembled => vec![(0, self.construct.to_vec())],
            KmerMasking::SegmentedByRegion => {
                let mut segments = Vec::new();
                // Backbone with every coding region masked out.
                let mut backbone = self.construct.to_vec();
                for &(start, end) in &self.cds_spans {
                    for byte in &mut backbone[start..end] {
                        *byte = b'N';
                    }
                }
                segments.push((0, backbone));
                for &(start, end) in &self.cds_spans {
                    segments.push((start, self.construct[start..end].to_vec()));
                }
                segments
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Violation {
    pub constraint_id: String,
    pub class: ConstraintClass,
    pub start_nt: usize,
    pub end_nt: usize,
    pub detail: String,
    /// Whether any editable position is involved, and therefore whether the
    /// search engine could do anything about it.
    ///
    /// A blocking violation confined to immutable sequence — a repeat between
    /// two backbone regions, a restriction site in the vector — is real, but it
    /// is a property of the *input*, not of the candidate. Counting it against
    /// every candidate would report every design as infeasible and give the
    /// search nothing to optimize toward. It belongs to the feasibility
    /// analyser (`P1-19` check 1, `CF-FEAS-001`), which reports it once at
    /// submission with the constraint to relax.
    ///
    /// This is also what the legacy prototype does, though implicitly: it only
    /// ever examines k-mers that occur in the coding sequence, so a
    /// backbone-only repeat is invisible to it. The differential test in
    /// `P1-27` is what surfaced the distinction.
    pub actionable: bool,
}

struct KmerRule {
    index: usize,
    k: usize,
    max_occurrences: u32,
    canonical: bool,
}

struct HomopolymerRule {
    index: usize,
    base: Option<u8>,
    max_run: usize,
}

struct GcRule {
    index: usize,
    window_nt: usize,
    min: f64,
    max: f64,
}

/// A compiled constraint set.
pub struct ConstraintChecker {
    set: ConstraintSet,
    cds_motifs: MotifAutomaton,
    construct_motifs: MotifAutomaton,
    kmer_rules: Vec<KmerRule>,
    homopolymer_rules: Vec<HomopolymerRule>,
    gc_rules: Vec<GcRule>,
    undetectable: Vec<String>,
}

impl std::fmt::Debug for ConstraintChecker {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ConstraintChecker({}, {} constraints, {} undetectable)",
            self.set.constraint_set_id,
            self.set.constraints.len(),
            self.undetectable.len()
        )
    }
}

impl ConstraintChecker {
    pub fn compile(set: &ConstraintSet) -> Result<Self> {
        set.validate()?;

        let mut cds_specs = Vec::new();
        let mut construct_specs = Vec::new();
        let mut kmer_rules = Vec::new();
        let mut homopolymer_rules = Vec::new();
        let mut gc_rules = Vec::new();
        let mut undetectable = Vec::new();

        for (index, constraint) in set.constraints.iter().enumerate() {
            if !constraint.enabled {
                continue;
            }
            if !constraint.rule.is_detectable_in_phase_1() {
                undetectable.push(constraint.id.clone());
                continue;
            }
            let specs = match constraint.scope {
                Scope::Cds | Scope::Cistron => &mut cds_specs,
                Scope::Construct => &mut construct_specs,
            };
            match &constraint.rule {
                Rule::ForbiddenMotif { motif, strand, .. } => {
                    specs.push((index, motif.clone(), *strand == Strand::Both));
                }
                Rule::RestrictionSite { site, strand, .. } => {
                    specs.push((index, site.clone(), *strand == Strand::Both));
                }
                Rule::CodonBan { codons, .. } => {
                    for codon in codons {
                        specs.push((index, codon.clone(), false));
                    }
                }
                Rule::KmerRepeat {
                    k,
                    max_occurrences,
                    canonical_strand,
                    ..
                } => kmer_rules.push(KmerRule {
                    index,
                    k: *k as usize,
                    max_occurrences: *max_occurrences,
                    canonical: *canonical_strand,
                }),
                Rule::Homopolymer { base, max_run } => homopolymer_rules.push(HomopolymerRule {
                    index,
                    base: match base {
                        HomopolymerBase::Any => None,
                        HomopolymerBase::A => Some(b'A'),
                        HomopolymerBase::C => Some(b'C'),
                        HomopolymerBase::G => Some(b'G'),
                        HomopolymerBase::T => Some(b'T'),
                    },
                    // `max_run` is the largest permitted run, so a violation is
                    // a run of max_run + 1.
                    max_run: *max_run as usize,
                }),
                Rule::GcWindow {
                    window_nt,
                    min,
                    max,
                } => gc_rules.push(GcRule {
                    index,
                    window_nt: *window_nt as usize,
                    min: *min,
                    max: *max,
                }),
                _ => undetectable.push(constraint.id.clone()),
            }
        }

        Ok(Self {
            set: set.clone(),
            cds_motifs: MotifAutomaton::build(&cds_specs)?,
            construct_motifs: MotifAutomaton::build(&construct_specs)?,
            kmer_rules,
            homopolymer_rules,
            gc_rules,
            undetectable,
        })
    }

    /// Constraints the Phase 1 checker cannot evaluate.
    pub fn undetectable(&self) -> &[String] {
        &self.undetectable
    }

    pub fn set(&self) -> &ConstraintSet {
        &self.set
    }

    fn constraint(&self, index: usize) -> &Constraint {
        &self.set.constraints[index]
    }

    /// Every violation, blocking and advisory.
    pub fn check(&self, view: &ConstructView) -> Vec<Violation> {
        let mut violations = Vec::new();
        self.check_motifs(view, &mut violations);
        self.check_kmers(view, &mut violations);
        self.check_homopolymers(view, &mut violations);
        self.check_gc(view, &mut violations);
        violations.sort_by(|a, b| {
            (a.start_nt, a.end_nt, &a.constraint_id).cmp(&(b.start_nt, b.end_nt, &b.constraint_id))
        });
        violations
    }

    /// Whether this candidate is buildable: no blocking violation the search
    /// could act on.
    ///
    /// Blocking violations confined to immutable sequence are excluded — see
    /// [`Violation::actionable`]. Use [`ConstraintChecker::check`] and filter on
    /// `!actionable` to find those; they are the feasibility analyser's job.
    pub fn is_feasible(&self, view: &ConstructView) -> bool {
        !self
            .check(view)
            .iter()
            .any(|v| v.class == ConstraintClass::Blocking && v.actionable)
    }

    /// How far beyond a changed span a violation can be created or destroyed.
    ///
    /// `None` means unbounded, which happens exactly when a k-mer repeat rule
    /// is compiled in: that rule counts occurrences over the whole construct,
    /// so editing three nucleotides here can push a k-mer somewhere else across
    /// its threshold. Every other Phase 1 rule has a bounded reach — a motif
    /// cannot be longer than its pattern, a GC window cannot be wider than
    /// `window_nt`, and a homopolymer run is bounded by finding its own ends.
    ///
    /// Callers use this to decide whether they can compute an exact local
    /// delta or must rescore. Returning `None` rather than a large number is
    /// the point: a caller that treated "unbounded" as "big" would produce
    /// deltas that are usually right, which is the worst of the three
    /// outcomes.
    pub fn influence_radius_nt(&self) -> Option<usize> {
        if !self.kmer_rules.is_empty() {
            return None;
        }
        let motifs = self
            .cds_motifs
            .max_pattern_len()
            .max(self.construct_motifs.max_pattern_len());
        let windows = self.gc_rules.iter().map(|r| r.window_nt).max().unwrap_or(0);
        Some(motifs.max(windows))
    }

    /// Violations intersecting `[start_nt, end_nt)`, computed over a bounded
    /// neighbourhood rather than the whole construct.
    ///
    /// Exact, not approximate: the neighbourhood is widened until every rule's
    /// reach is contained in it, so a violation this misses is one that could
    /// not have intersected the region. Returns `None` when
    /// [`ConstraintChecker::influence_radius_nt`] is `None`, because there is
    /// no such neighbourhood.
    ///
    /// The homopolymer case is why the neighbourhood is grown rather than
    /// simply offset by the radius. A run has no length bound, so a run
    /// intersecting the region can extend arbitrarily far past any fixed
    /// margin; truncating it would split one violation into two or hide it
    /// entirely. The window is therefore extended to the run's own ends, which
    /// costs the length of the run and nothing more.
    pub fn check_region(
        &self,
        view: &ConstructView,
        start_nt: usize,
        end_nt: usize,
    ) -> Option<Vec<Violation>> {
        let radius = self.influence_radius_nt()?;
        let len = view.construct.len();
        let start = start_nt.min(len);
        let end = end_nt.clamp(start, len);

        let mut lo = start.saturating_sub(radius);
        let mut hi = (end + radius).min(len);
        if !self.homopolymer_rules.is_empty() {
            while lo > 0 && view.construct[lo - 1] == view.construct[lo] {
                lo -= 1;
            }
            while hi < len && view.construct[hi - 1] == view.construct[hi] {
                hi += 1;
            }
        }

        let slice = &view.construct[lo..hi];
        let spans: Vec<(usize, usize)> = view
            .cds_spans
            .iter()
            .filter(|&&(s, e)| s < hi && lo < e)
            .map(|&(s, e)| (s.max(lo) - lo, e.min(hi) - lo))
            .collect();
        let sub = ConstructView::new(slice, spans)
            .with_masking(view.kmer_masking)
            .with_origin(view.origin_nt + lo);

        let mut local = Vec::new();
        self.check_motifs(&sub, &mut local);
        self.check_homopolymers(&sub, &mut local);
        self.check_gc(&sub, &mut local);

        let origin = view.origin_nt;
        let mut out: Vec<Violation> = local
            .into_iter()
            .filter(|v| v.start_nt < origin + end && origin + start < v.end_nt)
            .collect();
        out.sort_by(|a, b| {
            (a.start_nt, a.end_nt, &a.constraint_id).cmp(&(b.start_nt, b.end_nt, &b.constraint_id))
        });
        Some(out)
    }

    /// Blocking violations that lie entirely in immutable sequence.
    pub fn unfixable(&self, view: &ConstructView) -> Vec<Violation> {
        self.check(view)
            .into_iter()
            .filter(|v| v.class == ConstraintClass::Blocking && !v.actionable)
            .collect()
    }

    fn check_motifs(&self, view: &ConstructView, out: &mut Vec<Violation>) {
        let mut hits = Vec::new();

        if !self.construct_motifs.is_empty() {
            self.construct_motifs.find_all(view.construct, 0, &mut hits);
            self.tally_motifs(&self.construct_motifs, &hits, view, out);
            hits.clear();
        }
        if !self.cds_motifs.is_empty() {
            for &(start, end) in &view.cds_spans {
                self.cds_motifs
                    .find_all(&view.construct[start..end], start, &mut hits);
            }
            self.tally_motifs(&self.cds_motifs, &hits, view, out);
        }
    }

    fn tally_motifs(
        &self,
        automaton: &MotifAutomaton,
        hits: &[crate::motif::MotifHit],
        view: &ConstructView,
        out: &mut Vec<Violation>,
    ) {
        // Count per owning constraint, then report the occurrences beyond the
        // allowance. Reporting every hit for a `max_occurrences: 2` rule that
        // saw three would be three violations for one problem.
        let mut per_constraint: std::collections::BTreeMap<usize, Vec<&crate::motif::MotifHit>> =
            Default::default();
        for hit in hits {
            per_constraint
                .entry(automaton.patterns()[hit.pattern].constraint)
                .or_default()
                .push(hit);
        }
        for (index, mut occurrences) in per_constraint {
            let constraint = self.constraint(index);
            let allowed = match &constraint.rule {
                Rule::ForbiddenMotif {
                    max_occurrences, ..
                }
                | Rule::RestrictionSite {
                    max_occurrences, ..
                } => *max_occurrences as usize,
                _ => 0,
            };
            if occurrences.len() <= allowed {
                continue;
            }
            occurrences.sort_by_key(|hit| (hit.start, hit.end));
            for hit in occurrences.iter().skip(allowed) {
                let pattern = &automaton.patterns()[hit.pattern];
                let (start, end) = (view.origin_nt + hit.start, view.origin_nt + hit.end);
                out.push(Violation {
                    constraint_id: constraint.id.clone(),
                    class: constraint.class,
                    start_nt: start,
                    end_nt: end,
                    actionable: view.touches_editable(hit.start, hit.end),
                    detail: format!(
                        "{} occurs at {}..{} ({} strand); at most {allowed} permitted",
                        pattern.literal,
                        start,
                        end,
                        if pattern.reverse {
                            "antisense"
                        } else {
                            "sense"
                        }
                    ),
                });
            }
        }
    }

    fn check_kmers(&self, view: &ConstructView, out: &mut Vec<Violation>) {
        if self.kmer_rules.is_empty() {
            return;
        }
        let segments = view.kmer_segments();
        for rule in &self.kmer_rules {
            let constraint = self.constraint(rule.index);
            // Counts are summed across segments; positions are reported in
            // construct coordinates.
            let mut totals: std::collections::HashMap<u64, u32> = Default::default();
            let mut positions: std::collections::HashMap<u64, Vec<usize>> = Default::default();
            for (offset, bytes) in &segments {
                let Ok(index) = KmerIndex::build(bytes, rule.k, rule.canonical) else {
                    continue;
                };
                for (key, count) in index.counts() {
                    *totals.entry(*key).or_insert(0) += count;
                }
                if bytes.len() >= rule.k {
                    for pos in 0..=bytes.len() - rule.k {
                        if let Some(packed) = crate::kmer::pack_at(bytes, pos, rule.k) {
                            let key = if rule.canonical {
                                canonical_of(packed, rule.k)
                            } else {
                                packed
                            };
                            positions.entry(key).or_default().push(offset + pos);
                        }
                    }
                }
            }
            let mut offenders: Vec<(u64, u32)> = totals
                .into_iter()
                .filter(|(_, count)| *count > rule.max_occurrences)
                .collect();
            offenders.sort();
            for (key, count) in offenders {
                let mut where_ = positions.remove(&key).unwrap_or_default();
                where_.sort_unstable();
                let start = where_.first().copied().unwrap_or(0);
                out.push(Violation {
                    constraint_id: constraint.id.clone(),
                    class: constraint.class,
                    start_nt: start,
                    end_nt: start + rule.k,
                    // Actionable if ANY occurrence touches editable sequence:
                    // moving one occurrence is enough to bring the count down.
                    actionable: where_
                        .iter()
                        .any(|&pos| view.touches_editable(pos, pos + rule.k)),
                    detail: format!(
                        "{}-mer {} occurs {count} times (limit {}), at {:?}",
                        rule.k,
                        unpack(key, rule.k),
                        rule.max_occurrences,
                        where_
                    ),
                });
            }
        }
    }

    fn check_homopolymers(&self, view: &ConstructView, out: &mut Vec<Violation>) {
        for rule in &self.homopolymer_rules {
            let constraint = self.constraint(rule.index);
            for run in find_runs(view.construct, rule.max_run + 1, rule.base) {
                let (start, end) = (view.origin_nt + run.start, view.origin_nt + run.end);
                out.push(Violation {
                    constraint_id: constraint.id.clone(),
                    class: constraint.class,
                    start_nt: start,
                    end_nt: end,
                    actionable: view.touches_editable(run.start, run.end),
                    detail: format!(
                        "run of {} {}s at {}..{} exceeds the limit of {}",
                        run.len(),
                        run.base as char,
                        start,
                        end,
                        rule.max_run
                    ),
                });
            }
        }
    }

    fn check_gc(&self, view: &ConstructView, out: &mut Vec<Violation>) {
        if self.gc_rules.is_empty() {
            return;
        }
        let index = GcIndex::build(view.construct);
        for rule in &self.gc_rules {
            let constraint = self.constraint(rule.index);
            for (start, end, fraction) in index.windows_outside(rule.window_nt, rule.min, rule.max)
            {
                let (global_start, global_end) = (view.origin_nt + start, view.origin_nt + end);
                out.push(Violation {
                    constraint_id: constraint.id.clone(),
                    class: constraint.class,
                    start_nt: global_start,
                    end_nt: global_end,
                    actionable: view.touches_editable(start, end),
                    detail: format!(
                        "GC {fraction:.3} in {global_start}..{global_end} is outside [{:.2}, {:.2}]",
                        rule.min, rule.max
                    ),
                });
            }
        }
    }
}

fn canonical_of(packed: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    let mut value = packed;
    for _ in 0..k {
        rc = (rc << 2) | (3 - (value & 3));
        value >>= 2;
    }
    packed.min(rc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::RepairPolicy;

    fn set_with(rules: Vec<(&str, Scope, Rule)>) -> ConstraintSet {
        ConstraintSet {
            schema_version: Some("1.0.0".into()),
            constraint_set_id: "cst_test".into(),
            name: None,
            derived_from: vec![],
            constraints: rules
                .into_iter()
                .map(|(id, scope, rule)| Constraint {
                    id: id.into(),
                    class: ConstraintClass::Blocking,
                    scope,
                    cistron_id: None,
                    weight: 1.0,
                    enabled: true,
                    rule,
                })
                .collect(),
            repair_policy: RepairPolicy::default(),
            digest: None,
        }
    }

    #[test]
    fn unit_restriction_site_scope_matters() {
        let set = set_with(vec![(
            "c_bsai",
            Scope::Construct,
            Rule::RestrictionSite {
                enzyme: "BsaI".into(),
                site: "GGTCTC".into(),
                strand: Strand::Both,
                max_occurrences: 0,
            },
        )]);
        let checker = ConstraintChecker::compile(&set).unwrap();

        // Site straddling the CDS/backbone junction: visible at construct scope.
        let construct = b"AAAAGGTCTCAAAA";
        let view = ConstructView::new(construct, vec![(6, 14)]);
        assert_eq!(checker.check(&view).len(), 1);
        assert!(!checker.is_feasible(&view));

        // The same rule scoped to the CDS cannot see it.
        let mut cds_set = set.clone();
        cds_set.constraints[0].scope = Scope::Cds;
        let cds_checker = ConstraintChecker::compile(&cds_set).unwrap();
        assert!(cds_checker.check(&view).is_empty());
    }

    #[test]
    fn unit_kmer_masking_changes_junction_visibility() {
        let set = set_with(vec![(
            "c_kmer",
            Scope::Construct,
            Rule::KmerRepeat {
                k: 4,
                max_occurrences: 1,
                count_overlapping: true,
                canonical_strand: false,
            },
        )]);
        let checker = ConstraintChecker::compile(&set).unwrap();

        // ACGT twice: once wholly inside the backbone at 0..4, once straddling
        // the CDS boundary at 6..10. The flanks are chosen so that no other
        // 4-mer repeats, which keeps the assertion about exactly one offender.
        //
        //   index  0 1 2 3 4 5 6 7 | 8 9 10 11
        //   base   A C G T A A A C | G T  C  C
        //          \--ACGT--/   \--ACGT--/
        let construct = b"ACGTAAACGTCC";
        let cds = vec![(8usize, 12usize)];

        let assembled = ConstructView::new(construct, cds.clone());
        let violations = checker.check(&assembled);
        assert_eq!(
            violations.len(),
            1,
            "the assembled view sees both occurrences: {violations:?}"
        );
        assert!(violations[0].detail.contains("ACGT"));

        let segmented =
            ConstructView::new(construct, cds).with_masking(KmerMasking::SegmentedByRegion);
        assert!(
            checker.check(&segmented).is_empty(),
            "the legacy segmented view is blind to the straddling occurrence"
        );
    }

    /// A violation the search cannot act on is reported, but is not counted
    /// against the candidate.
    ///
    /// Surfaced by the `P1-27` differential: the legacy prototype only examines
    /// k-mers that occur in the coding sequence, so a repeat between two
    /// backbone regions is invisible to it. That turned out to be the right
    /// question to ask of a *candidate*; the violation itself belongs to the
    /// feasibility analyser.
    #[test]
    fn unit_backbone_only_violation_is_reported_but_not_actionable() {
        let set = set_with(vec![(
            "c_kmer",
            Scope::Construct,
            Rule::KmerRepeat {
                k: 4,
                max_occurrences: 1,
                count_overlapping: true,
                canonical_strand: false,
            },
        )]);
        let checker = ConstraintChecker::compile(&set).unwrap();

        // ACGT twice, at 0..4 and 10..14, with the coding region at 6..10 in
        // between. Flanks chosen so ACGT is the only repeated 4-mer.
        //
        //   index  0 1 2 3 4 5 | 6 7 8 9 | 10 11 12 13 14 15
        //   base   A C G T T T | G G C C |  A  C  G  T  A  A
        //          \--ACGT--/    \-CDS-/    \---ACGT---/
        let construct = b"ACGTTTGGCCACGTAA";
        let view = ConstructView::new(construct, vec![(6, 10)]);

        let violations = checker.check(&view);
        assert_eq!(violations.len(), 1, "the violation is still reported");
        assert!(
            !violations[0].actionable,
            "no occurrence touches an editable position, so nothing can fix it"
        );
        assert!(
            checker.is_feasible(&view),
            "a candidate is not infeasible because of a pre-existing backbone repeat"
        );
        assert_eq!(checker.unfixable(&view).len(), 1);

        // Move the coding region so that one occurrence overlaps it: now the
        // same violation is the candidate's to fix.
        let overlapping = ConstructView::new(construct, vec![(10, 14)]);
        let violations = checker.check(&overlapping);
        assert_eq!(violations.len(), 1);
        assert!(violations[0].actionable);
        assert!(!checker.is_feasible(&overlapping));
    }

    #[test]
    fn unit_homopolymer_and_gc() {
        let set = set_with(vec![
            (
                "c_hp",
                Scope::Construct,
                Rule::Homopolymer {
                    base: HomopolymerBase::Any,
                    max_run: 5,
                },
            ),
            (
                "c_gc",
                Scope::Construct,
                Rule::GcWindow {
                    window_nt: 8,
                    min: 0.2,
                    max: 0.8,
                },
            ),
        ]);
        let checker = ConstraintChecker::compile(&set).unwrap();

        let ok = ConstructView::new(b"ACGTACGTACGTACGT", vec![]);
        assert!(checker.check(&ok).is_empty());

        let run = ConstructView::new(b"ACGTAAAAAAGTACGT", vec![]);
        let violations = checker.check(&run);
        assert!(violations.iter().any(|v| v.constraint_id == "c_hp"));

        let gc = ConstructView::new(b"GGGGGGGGGGGGGGGG", vec![]);
        let violations = checker.check(&gc);
        assert!(violations.iter().any(|v| v.constraint_id == "c_gc"));
    }

    #[test]
    fn unit_undetectable_rules_are_listed_not_ignored() {
        let example: serde_json::Value = serde_json::from_str(include_str!(
            "../../../planning/02-schemas/examples/constraint_set.example.json"
        ))
        .unwrap();
        let set = ConstraintSet::from_json(&example).unwrap();
        let checker = ConstraintChecker::compile(&set).unwrap();
        assert_eq!(checker.undetectable(), ["c_initiation_structure"]);
    }

    #[test]
    fn unit_max_occurrences_allowance_is_respected() {
        let set = set_with(vec![(
            "c_motif",
            Scope::Construct,
            Rule::ForbiddenMotif {
                motif: "ACGT".into(),
                strand: Strand::Sense,
                max_occurrences: 2,
            },
        )]);
        let checker = ConstraintChecker::compile(&set).unwrap();
        // Two occurrences: allowed.
        assert!(checker
            .check(&ConstructView::new(b"ACGTTTACGT", vec![]))
            .is_empty());
        // Three: one violation, reported at the third.
        let violations = checker.check(&ConstructView::new(b"ACGTTTACGTTTACGT", vec![]));
        assert_eq!(violations.len(), 1);
        assert_eq!(violations[0].start_nt, 12);
    }
}
