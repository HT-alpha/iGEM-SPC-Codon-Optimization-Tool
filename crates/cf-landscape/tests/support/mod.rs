//! Fixtures shared by the `P1-21` and `P1-22` tests.

use cf_constraints::spec::{
    Constraint, ConstraintClass, ConstraintSet, HomopolymerBase, RepairPolicy, Rule, Scope, Strand,
};
use cf_landscape::{CodonUsageTable, ConstructLayout};
use cf_landscape_api::{Budget, BudgetTier, EvalContext, HostContext};
use cf_types::{Codon, CodonSequence, Digest, GeneticCode, HostId};
use serde_json::{json, Map, Value};

/// A coding sequence with real synonymous choice at most positions.
///
/// `M-L-R-S-G-V-A-P-T-K-*`: leucine, arginine and serine have six codons each,
/// and the rest have two or four, so a delta table has something to iterate.
pub const CDS: &str = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAATAA";

pub fn sequence() -> CodonSequence {
    CodonSequence::from_cds(CDS, 1).unwrap()
}

/// The same coding sequence inside a plasmid backbone.
pub fn flanked_sequence() -> CodonSequence {
    sequence()
        .with_five_prime_flank("TTGACAATTAATCATCGGCTCGTATAATGTGTGGA")
        .unwrap()
        .with_three_prime_flank("AAAGCCCGCTCTAGAGGATCCCCG")
        .unwrap()
}

pub fn host() -> HostContext {
    HostContext::new(
        HostId::new("escherichia_coli:511145/BL21").unwrap(),
        Digest::of_bytes(b"host-profile-fixture"),
    )
}

pub fn ctx<'a>(sequence: &'a CodonSequence, host: &'a HostContext) -> EvalContext<'a> {
    EvalContext::new(sequence, host, Budget::new(BudgetTier::Standard), 0)
}

/// A deterministic but non-uniform codon usage table over all 64 codons.
///
/// Built by assigning each codon a pseudo-random raw frequency and normalizing
/// within its synonymous family, which is how `relative_adaptiveness` is
/// defined. Non-uniform on purpose: a table where every weight was 1 would make
/// every CAI delta zero and the delta tests vacuous.
pub fn usage_table_json(subset: &str) -> Value {
    let code = GeneticCode::get(1).unwrap();
    let mut raw = [0.0f64; 64];
    for codon in Codon::all() {
        let index = codon.index_usize();
        // Spread across (0, 1] with no zeros: an unobserved codon is a separate
        // case with its own fixture.
        raw[index] = ((index as f64 * 37.0) % 61.0 + 1.0) / 61.0;
    }

    let mut table = Map::new();
    for codon in Codon::all() {
        let family = code.synonyms(codon);
        let max = family
            .iter()
            .map(|c| raw[c.index_usize()])
            .fold(f64::MIN, f64::max);
        table.insert(
            codon.to_string_upper(),
            json!({
                "aa": (code.amino_acid(codon) as char).to_string(),
                "count": (raw[codon.index_usize()] * 1000.0) as u64,
                "relative_adaptiveness": raw[codon.index_usize()] / max,
            }),
        );
    }

    json!({
        "codon_usage": {
            "subset": subset,
            "n_cds": 128,
            "table": Value::Object(table),
        }
    })
}

pub fn usage_table(subset: &str) -> CodonUsageTable {
    CodonUsageTable::from_host_profile(&usage_table_json(subset), 1).unwrap()
}

pub fn layout() -> ConstructLayout {
    ConstructLayout::from_sequence(&flanked_sequence())
}

fn set_with(rules: Vec<(&str, Rule)>) -> ConstraintSet {
    ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_landscape_test".into(),
        name: None,
        derived_from: vec![],
        constraints: rules
            .into_iter()
            .map(|(id, rule)| Constraint {
                id: id.into(),
                class: ConstraintClass::Blocking,
                scope: Scope::Construct,
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

/// Motif, homopolymer and GC rules — every reach bounded, so the violation
/// counter is L1.
pub fn local_constraints() -> ConstraintSet {
    set_with(vec![
        (
            "no_bamhi",
            Rule::RestrictionSite {
                enzyme: "BamHI".into(),
                site: "GGATCC".into(),
                strand: Strand::Both,
                max_occurrences: 0,
            },
        ),
        (
            "no_ggcgc",
            Rule::ForbiddenMotif {
                motif: "GGCGC".into(),
                strand: Strand::Both,
                max_occurrences: 0,
            },
        ),
        (
            "runs",
            Rule::Homopolymer {
                base: HomopolymerBase::Any,
                max_run: 3,
            },
        ),
    ])
}

/// The same rules plus a k-mer repeat limit, whose reach is the whole
/// construct — so the violation counter degrades to L0.
pub fn global_constraints() -> ConstraintSet {
    let mut set = local_constraints();
    set.constraints.push(Constraint {
        id: "repeats".into(),
        class: ConstraintClass::Blocking,
        scope: Scope::Construct,
        cistron_id: None,
        weight: 1.0,
        enabled: true,
        rule: Rule::KmerRepeat {
            k: 8,
            max_occurrences: 2,
            count_overlapping: true,
            canonical_strand: true,
        },
    });
    set
}
