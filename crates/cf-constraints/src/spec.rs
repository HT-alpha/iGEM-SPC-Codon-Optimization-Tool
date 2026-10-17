//! Constraint-set specification, parsing, and assembly/vendor expansion.
//!
//! Ticket P1-13. Contract: `planning/02-schemas/constraint_set.schema.json`.
//!
//! TWO CLASSES, AND THE DISTINCTION IS OPERATIONAL. A `blocking` constraint
//! makes the construct unbuildable; an `advisory` one is a preference that
//! participates in the objective vector instead. Conflating them is the classic
//! way a codon optimizer returns nothing, so the type keeps them apart and the
//! checker reports them separately.
//!
//! EVERY AUTO-ADDED CONSTRAINT IS REPORTED. Expanding an assembly standard into
//! concrete enzyme sites is convenient and invisible, and an invisible
//! constraint that changes the answer is inexplicable to the user. `expand`
//! therefore returns an [`ExpansionReport`] naming every rule it added and why.

use std::collections::BTreeSet;

use cf_types::error::{code, CfError, Result};
use serde::{Deserialize, Serialize};

/// Which strand a pattern is matched on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Strand {
    Sense,
    #[default]
    Both,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ConstraintClass {
    Blocking,
    Advisory,
}

/// Where a constraint is evaluated.
///
/// The distinction is load-bearing for k-mer repeats: counting them over the
/// coding region alone misses repeats against the immutable backbone, which is
/// what a synthesis vendor actually rejects.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Scope {
    Cds,
    #[default]
    Construct,
    Cistron,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "UPPERCASE")]
pub enum HomopolymerBase {
    A,
    C,
    G,
    T,
    #[serde(rename = "any")]
    #[default]
    Any,
}

fn default_true() -> bool {
    true
}
fn default_zero() -> u32 {
    0
}

/// The rule payload, discriminated by `kind` exactly as in the schema.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Rule {
    ForbiddenMotif {
        motif: String,
        #[serde(default)]
        strand: Strand,
        #[serde(default = "default_zero")]
        max_occurrences: u32,
    },
    RestrictionSite {
        enzyme: String,
        site: String,
        #[serde(default)]
        strand: Strand,
        #[serde(default = "default_zero")]
        max_occurrences: u32,
    },
    Homopolymer {
        #[serde(default)]
        base: HomopolymerBase,
        max_run: u32,
    },
    KmerRepeat {
        k: u32,
        max_occurrences: u32,
        #[serde(default = "default_true")]
        count_overlapping: bool,
        #[serde(default)]
        canonical_strand: bool,
    },
    GcWindow {
        window_nt: u32,
        min: f64,
        max: f64,
    },
    /// Detection lands with the structure engine in Phase 2 (`P2-11`); parsed
    /// here so a Phase 1 constraint set round-trips without loss.
    LocalStructure {
        region: String,
        #[serde(default)]
        span: Option<serde_json::Value>,
        metric: String,
        #[serde(default)]
        min: Option<f64>,
        #[serde(default)]
        max: Option<f64>,
        #[serde(default)]
        engine: Option<String>,
    },
    CodonBan {
        codons: Vec<String>,
        #[serde(default)]
        reason: Option<String>,
    },
    MinCodonFrequency {
        min_relative_adaptiveness: f64,
    },
    CustomRegex {
        pattern: String,
        #[serde(default = "default_zero")]
        max_occurrences: u32,
        #[serde(default)]
        flags: Vec<String>,
    },
}

impl Rule {
    /// Whether Phase 1 can detect violations of this rule.
    ///
    /// Reported rather than silently ignored: a constraint the checker cannot
    /// evaluate must not be presented as satisfied. `custom_regex` is
    /// detectable: a pattern that compiles under the size cap is enforced, and
    /// one that does not compile is refused at `ConstraintChecker::compile`
    /// (`CF-FEAS-006`, `CF-FEAS-007`) rather than treated as satisfied.
    /// `local_structure` and `min_codon_frequency` stay out of this phase.
    pub fn is_detectable_in_phase_1(&self) -> bool {
        !matches!(
            self,
            Rule::LocalStructure { .. } | Rule::MinCodonFrequency { .. }
        )
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Constraint {
    pub id: String,
    pub class: ConstraintClass,
    #[serde(default)]
    pub scope: Scope,
    #[serde(default)]
    pub cistron_id: Option<String>,
    #[serde(default = "one")]
    pub weight: f64,
    #[serde(default = "default_true")]
    pub enabled: bool,
    pub rule: Rule,
}

fn one() -> f64 {
    1.0
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RepairPolicy {
    #[serde(default = "twelve")]
    pub max_window_codons: u32,
    #[serde(default = "two")]
    pub max_escalations: u32,
    #[serde(default = "three")]
    pub max_depth: u32,
    #[serde(default = "thirty_two")]
    pub beam_width: u32,
    #[serde(default = "five")]
    pub tabu_tenure_generations: u32,
    #[serde(default = "default_true")]
    pub objective_aware: bool,
    #[serde(default = "default_true")]
    pub stochastic_tie_break: bool,
    #[serde(default = "default_true")]
    pub lamarckian: bool,
}

fn twelve() -> u32 {
    12
}
fn two() -> u32 {
    2
}
fn three() -> u32 {
    3
}
fn thirty_two() -> u32 {
    32
}
fn five() -> u32 {
    5
}

impl Default for RepairPolicy {
    fn default() -> Self {
        Self {
            max_window_codons: 12,
            max_escalations: 2,
            max_depth: 3,
            beam_width: 32,
            tabu_tenure_generations: 5,
            objective_aware: true,
            stochastic_tie_break: true,
            lamarckian: true,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConstraintSet {
    #[serde(default)]
    pub schema_version: Option<String>,
    pub constraint_set_id: String,
    #[serde(default)]
    pub name: Option<String>,
    #[serde(default)]
    pub derived_from: Vec<String>,
    pub constraints: Vec<Constraint>,
    #[serde(default)]
    pub repair_policy: RepairPolicy,
    #[serde(default)]
    pub digest: Option<String>,
}

impl ConstraintSet {
    pub fn from_json(value: &serde_json::Value) -> Result<Self> {
        serde_json::from_value(value.clone()).map_err(|e| {
            CfError::new(
                code::FEAS_CONSTRAINT_INVALID,
                format!("constraint set failed to parse: {e}"),
            )
        })
    }

    pub fn to_json(&self) -> Result<serde_json::Value> {
        serde_json::to_value(self).map_err(|e| {
            CfError::new(
                code::FEAS_CONSTRAINT_INVALID,
                format!("constraint set failed to serialize: {e}"),
            )
        })
    }

    /// Validate structural properties the JSON Schema cannot express.
    pub fn validate(&self) -> Result<()> {
        let mut seen = BTreeSet::new();
        for constraint in &self.constraints {
            if !seen.insert(constraint.id.as_str()) {
                return Err(CfError::new(
                    code::FEAS_CONSTRAINT_INVALID,
                    format!("duplicate constraint id {:?}", constraint.id),
                ));
            }
            match &constraint.rule {
                Rule::GcWindow { min, max, .. } if min > max => {
                    return Err(CfError::new(
                        code::FEAS_CONSTRAINT_INVALID,
                        format!("{}: GC band min {min} exceeds max {max}", constraint.id),
                    ));
                }
                Rule::KmerRepeat { k, .. } if !(4..=32).contains(k) => {
                    return Err(CfError::new(
                        code::FEAS_CONSTRAINT_INVALID,
                        format!(
                            "{}: k={k} is outside 4..=32; k-mers are packed two bits per base",
                            constraint.id
                        ),
                    ));
                }
                Rule::Homopolymer { max_run, .. } if *max_run < 2 => {
                    return Err(CfError::new(
                        code::FEAS_CONSTRAINT_INVALID,
                        format!("{}: max_run must be at least 2", constraint.id),
                    ));
                }
                _ => {}
            }
        }
        Ok(())
    }

    pub fn blocking(&self) -> impl Iterator<Item = &Constraint> {
        self.constraints
            .iter()
            .filter(|c| c.enabled && c.class == ConstraintClass::Blocking)
    }
}

// --------------------------------------------------------------------------- //
// assembly and vendor expansion
// --------------------------------------------------------------------------- //

/// Type IIS and common type II recognition sequences.
///
/// Recognition sequences only — cut positions are irrelevant here, because the
/// constraint is "this site must not occur", not "cut here". Restricted to
/// enzymes whose sites are unambiguous and widely tabulated; an enzyme absent
/// from this table is reported rather than guessed at.
pub static ENZYME_SITES: &[(&str, &str)] = &[
    ("BsaI", "GGTCTC"),
    ("BsmBI", "CGTCTC"),
    ("Esp3I", "CGTCTC"),
    ("BbsI", "GAAGAC"),
    ("SapI", "GCTCTTC"),
    ("EcoRI", "GAATTC"),
    ("BamHI", "GGATCC"),
    ("HindIII", "AAGCTT"),
    ("NotI", "GCGGCCGC"),
    ("XbaI", "TCTAGA"),
    ("SpeI", "ACTAGT"),
    ("PstI", "CTGCAG"),
    ("XhoI", "CTCGAG"),
    ("SalI", "GTCGAC"),
    ("KpnI", "GGTACC"),
    ("SacI", "GAGCTC"),
    ("NcoI", "CCATGG"),
    ("NdeI", "CATATG"),
    ("BglII", "AGATCT"),
    ("EcoRV", "GATATC"),
    ("PvuII", "CAGCTG"),
    ("SmaI", "CCCGGG"),
    ("ApaI", "GGGCCC"),
];

pub fn enzyme_site(enzyme: &str) -> Option<&'static str> {
    ENZYME_SITES
        .iter()
        .find(|(name, _)| name.eq_ignore_ascii_case(enzyme))
        .map(|(_, site)| *site)
}

/// One line of the report handed back to the user.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpansionNote {
    pub constraint_id: String,
    pub origin: String,
    pub detail: String,
}

#[derive(Clone, Debug, Default)]
pub struct ExpansionReport {
    pub added: Vec<ExpansionNote>,
    pub warnings: Vec<String>,
}

/// Assembly method, mirroring `construct_context.assembly.standard`.
#[derive(Clone, Debug, Default)]
pub struct AssemblySpec {
    pub standard: String,
    pub forbidden_enzymes: Vec<String>,
    pub vendor_profile: Option<String>,
}

/// Expand an assembly standard and vendor profile into concrete rules.
pub fn expand(set: &mut ConstraintSet, assembly: &AssemblySpec) -> Result<ExpansionReport> {
    let mut report = ExpansionReport::default();
    let existing: BTreeSet<String> = set.constraints.iter().map(|c| c.id.clone()).collect();

    let mut enzymes: Vec<String> = match assembly.standard.as_str() {
        "golden_gate_bsai" => vec!["BsaI".into()],
        "golden_gate_bsmbi" => vec!["BsmBI".into()],
        "moclo" => vec!["BsaI".into(), "BsmBI".into()],
        "biobrick_rfc10" => vec![
            "EcoRI".into(),
            "XbaI".into(),
            "SpeI".into(),
            "PstI".into(),
            "NotI".into(),
        ],
        "gibson" | "gene_synthesis_only" | "none" => Vec::new(),
        "restriction_ligation" => Vec::new(),
        other => {
            return Err(CfError::new(
                code::FEAS_CONSTRAINT_INVALID,
                format!("unknown assembly standard {other:?}"),
            ));
        }
    };
    enzymes.extend(assembly.forbidden_enzymes.iter().cloned());
    enzymes.sort();
    enzymes.dedup();

    for enzyme in enzymes {
        let Some(site) = enzyme_site(&enzyme) else {
            report.warnings.push(format!(
                "no recognition sequence is vendored for {enzyme:?}; supply an explicit \
                 restriction_site or forbidden_motif rule instead of relying on expansion"
            ));
            continue;
        };
        let id = format!("c_auto_{}", enzyme.to_ascii_lowercase());
        if existing.contains(&id) {
            continue;
        }
        set.constraints.push(Constraint {
            id: id.clone(),
            class: ConstraintClass::Blocking,
            scope: Scope::Construct,
            cistron_id: None,
            weight: 1.0,
            enabled: true,
            rule: Rule::RestrictionSite {
                enzyme: enzyme.clone(),
                site: site.to_string(),
                strand: Strand::Both,
                max_occurrences: 0,
            },
        });
        report.added.push(ExpansionNote {
            constraint_id: id,
            origin: format!("assembly_standard={}", assembly.standard),
            detail: format!("{enzyme} site {site} forbidden on both strands"),
        });
    }

    // Vendor synthesis limits are deliberately NOT vendored per vendor. Real
    // limits are negotiated, versioned and occasionally confidential, and a
    // plausible-looking wrong number here would silently change every design.
    // A named vendor therefore resolves to a conservative generic profile and
    // says so.
    if let Some(vendor) = &assembly.vendor_profile {
        let id = "c_auto_vendor_homopolymer".to_string();
        if !existing.contains(&id) {
            set.constraints.push(Constraint {
                id: id.clone(),
                class: ConstraintClass::Blocking,
                scope: Scope::Construct,
                cistron_id: None,
                weight: 1.0,
                enabled: true,
                rule: Rule::Homopolymer {
                    base: HomopolymerBase::Any,
                    max_run: 8,
                },
            });
            report.added.push(ExpansionNote {
                constraint_id: id,
                origin: format!("vendor_profile={vendor} (generic)"),
                detail: "homopolymer runs longer than 8 nt forbidden".into(),
            });
        }
        report.warnings.push(format!(
            "vendor-specific synthesis limits are not vendored for {vendor:?}; a conservative \
             generic profile was applied. Supply the vendor's actual limits as explicit \
             constraints before relying on a feasibility result."
        ));
    }

    set.derived_from.push("assembly_standard".into());
    if assembly.vendor_profile.is_some() {
        set.derived_from.push("vendor_profile".into());
    }
    set.derived_from.sort();
    set.derived_from.dedup();
    set.validate()?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn example() -> serde_json::Value {
        let raw = include_str!("../../../planning/02-schemas/examples/constraint_set.example.json");
        serde_json::from_str(raw).unwrap()
    }

    /// contract_*: the worked example parses, round-trips, and validates.
    #[test]
    fn contract_constraint_set_example_round_trips() {
        let value = example();
        let set = ConstraintSet::from_json(&value).unwrap();
        set.validate().unwrap();
        assert_eq!(set.constraint_set_id, "cst_goldengate_twist_default");
        assert_eq!(set.constraints.len(), 9);
        assert_eq!(set.blocking().count(), 7);

        // Every rule kind in the example survives the round trip unchanged.
        let back = set.to_json().unwrap();
        let reparsed = ConstraintSet::from_json(&back).unwrap();
        assert_eq!(set, reparsed);
    }

    #[test]
    fn contract_overlapping_kmer_semantics_are_explicit_in_the_set() {
        // The legacy prototype counts overlapping occurrences; the schema makes
        // that a declared field rather than an inherited accident, and the
        // example must say so.
        let set = ConstraintSet::from_json(&example()).unwrap();
        let kmer_rules: Vec<&Rule> = set
            .constraints
            .iter()
            .map(|c| &c.rule)
            .filter(|r| matches!(r, Rule::KmerRepeat { .. }))
            .collect();
        assert_eq!(kmer_rules.len(), 2);
        for rule in kmer_rules {
            let Rule::KmerRepeat {
                count_overlapping, ..
            } = rule
            else {
                unreachable!()
            };
            assert!(*count_overlapping, "overlapping counting must be explicit");
        }
    }

    #[test]
    fn unit_expansion_reports_every_added_constraint() {
        let mut set = ConstraintSet {
            schema_version: Some("1.0.0".into()),
            constraint_set_id: "cst_test".into(),
            name: None,
            derived_from: vec![],
            constraints: vec![],
            repair_policy: RepairPolicy::default(),
            digest: None,
        };
        let report = expand(
            &mut set,
            &AssemblySpec {
                standard: "moclo".into(),
                forbidden_enzymes: vec!["NotI".into()],
                vendor_profile: Some("twist".into()),
            },
        )
        .unwrap();

        let ids: Vec<&str> = report
            .added
            .iter()
            .map(|n| n.constraint_id.as_str())
            .collect();
        assert!(ids.contains(&"c_auto_bsai"));
        assert!(ids.contains(&"c_auto_bsmbi"));
        assert!(ids.contains(&"c_auto_noti"));
        assert!(ids.contains(&"c_auto_vendor_homopolymer"));
        assert_eq!(set.constraints.len(), report.added.len());
        // The vendor caveat is surfaced, not swallowed.
        assert!(report.warnings.iter().any(|w| w.contains("not vendored")));
        assert!(set.derived_from.contains(&"assembly_standard".to_string()));
    }

    #[test]
    fn unit_unknown_enzyme_warns_rather_than_guessing() {
        let mut set = ConstraintSet {
            schema_version: None,
            constraint_set_id: "cst_test".into(),
            name: None,
            derived_from: vec![],
            constraints: vec![],
            repair_policy: RepairPolicy::default(),
            digest: None,
        };
        let report = expand(
            &mut set,
            &AssemblySpec {
                standard: "none".into(),
                forbidden_enzymes: vec!["NotARealEnzyme".into()],
                vendor_profile: None,
            },
        )
        .unwrap();
        assert!(report.added.is_empty());
        assert_eq!(report.warnings.len(), 1);
        assert!(report.warnings[0].contains("NotARealEnzyme"));
    }

    #[test]
    fn unit_validation_rejects_malformed_sets() {
        let mut set = ConstraintSet::from_json(&example()).unwrap();
        set.constraints.push(set.constraints[0].clone());
        assert_eq!(
            set.validate().unwrap_err().code,
            code::FEAS_CONSTRAINT_INVALID
        );

        let mut set = ConstraintSet::from_json(&example()).unwrap();
        set.constraints[0].rule = Rule::KmerRepeat {
            k: 64,
            max_occurrences: 1,
            count_overlapping: true,
            canonical_strand: false,
        };
        assert!(set.validate().is_err());
    }

    #[test]
    fn unit_undetectable_rules_are_identifiable() {
        let set = ConstraintSet::from_json(&example()).unwrap();
        let deferred: Vec<&str> = set
            .constraints
            .iter()
            .filter(|c| !c.rule.is_detectable_in_phase_1())
            .map(|c| c.id.as_str())
            .collect();
        // The example's local_structure rule needs the Phase 2 structure engine.
        assert_eq!(deferred, vec!["c_initiation_structure"]);
    }
}
