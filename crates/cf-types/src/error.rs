//! The `CF-*` error code registry and RFC 9457 problem details.
//!
//! Ticket P1-04. Contract: `planning/04-api/01-rest-api.md` section 2.6.
//!
//! Every code the system can emit is declared once, in the `codes!` macro
//! below, which generates both the constant and its registry entry. A caller
//! that uses the constants therefore cannot emit an undocumented code, and
//! `TST-API-004` (registry portion) additionally asserts that every code cited
//! anywhere in the planning package exists here — so the plan and the code
//! cannot drift apart in either direction.

use std::borrow::Cow;
use std::fmt;

/// One registry entry: the wire code, its HTTP status, and a short title.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CfCode {
    pub code: &'static str,
    pub status: u16,
    pub title: &'static str,
    /// Where the condition is specified. Kept in the binary because an error a
    /// user cannot trace back to a decision is an error they cannot act on.
    pub doc: &'static str,
}

macro_rules! codes {
    ($( $ident:ident => ($code:literal, $status:literal, $title:literal, $doc:literal) ),* $(,)?) => {
        /// Code constants. Use these rather than string literals.
        pub mod code {
            $( pub const $ident: &str = $code; )*
        }

        /// Every code the system may emit.
        pub static REGISTRY: &[CfCode] = &[
            $( CfCode { code: $code, status: $status, title: $title, doc: $doc } ),*
        ];
    };
}

codes! {
    // ---- CF-SEQ: input sequence and the canonical representation ----------
    SEQ_NOT_IN_FRAME       => ("CF-SEQ-001", 422, "Coding sequence is not a whole number of codons", "planning/02-schemas/sequence.schema.json"),
    SEQ_INVALID_BASE       => ("CF-SEQ-002", 422, "Sequence contains a character that is not A, C, G or T", "planning/02-schemas/sequence.schema.json"),
    SEQ_INTERNAL_STOP      => ("CF-SEQ-003", 422, "Coding sequence contains an internal stop codon", "planning/02-schemas/sequence.schema.json"),
    SEQ_NOT_SYNONYMOUS     => ("CF-SEQ-004", 422, "Substituted codon is not synonymous with the codon it replaces", "planning/02-schemas/landscape_query.schema.json"),
    SEQ_INDEX_OUT_OF_RANGE => ("CF-SEQ-005", 422, "Codon index is outside the sequence", "planning/02-schemas/sequence.schema.json"),
    SEQ_POSITION_LOCKED    => ("CF-SEQ-006", 422, "Codon position is locked and may not be edited", "planning/02-schemas/sequence.schema.json"),
    SEQ_UNKNOWN_TABLE      => ("CF-SEQ-007", 422, "Translation table id is not defined by NCBI", "crates/cf-types/data/genetic_codes.json"),
    SEQ_PROTEIN_MISMATCH   => ("CF-SEQ-008", 422, "Declared protein does not match the translated coding sequence", "planning/02-schemas/sequence.schema.json"),
    SEQ_LOCK_OUT_OF_RANGE  => ("CF-SEQ-009", 422, "Locked span lies outside the sequence", "planning/02-schemas/sequence.schema.json"),
    SEQ_DIGEST_MISMATCH    => ("CF-SEQ-010", 422, "Declared digest does not match the document content", "planning/02-schemas/sequence.schema.json"),
    SEQ_STOP_PLACEMENT     => ("CF-SEQ-011", 422, "Stop codon placement disagrees with includes_stop", "planning/02-schemas/sequence.schema.json"),
    SEQ_EMPTY              => ("CF-SEQ-012", 422, "Coding sequence is empty", "planning/02-schemas/sequence.schema.json"),
    SEQ_BAD_IDENTIFIER     => ("CF-SEQ-013", 422, "Identifier does not match its declared pattern", "planning/02-schemas/common.schema.json"),
    SEQ_NOT_CANONICAL      => ("CF-SEQ-014", 500, "Value cannot be serialized as canonical JSON", "planning/02-schemas/common.schema.json"),

    // ---- CF-HOST: host reference data --------------------------------------
    HOST_UNKNOWN           => ("CF-HOST-001", 404, "Unknown host", "planning/02-schemas/host_profile.schema.json"),
    HOST_PROFILE_INVALID   => ("CF-HOST-002", 422, "Host profile failed validation", "planning/02-schemas/host_profile.schema.json"),
    HOST_TERM_UNSUPPORTED  => ("CF-HOST-003", 422, "Host has no data to support a requested term", "planning/02-schemas/host_profile.schema.json"),

    // ---- CF-FEAS: the five constraint-feasibility checks --------------------
    FEAS_LOCKED_SPAN       => ("CF-FEAS-001", 422, "A blocking violation lies inside immutable sequence", "planning/03-algorithms/06-constraint-repair.md"),
    FEAS_COUNTING_BOUND    => ("CF-FEAS-002", 422, "k-mer constraint is unsatisfiable by the pigeonhole bound", "planning/03-algorithms/06-constraint-repair.md"),
    FEAS_PER_POSITION      => ("CF-FEAS-003", 422, "A position has no admissible synonym", "planning/03-algorithms/06-constraint-repair.md"),
    FEAS_UNREACHABLE       => ("CF-FEAS-004", 422, "Motif constraints jointly close off every continuation", "planning/03-algorithms/06-constraint-repair.md"),
    FEAS_GC_UNATTAINABLE   => ("CF-FEAS-005", 422, "Required GC band is outside the attainable range", "planning/03-algorithms/06-constraint-repair.md"),
    FEAS_CONSTRAINT_INVALID=> ("CF-FEAS-006", 422, "Constraint definition is malformed", "planning/02-schemas/constraint_set.schema.json"),
    FEAS_REGEX_TOO_LARGE   => ("CF-FEAS-007", 422, "Custom pattern exceeds the compiled-size cap", "planning/06-build/06-risk-register.md"),

    // ---- CF-RESOLVE: model resolution --------------------------------------
    RESOLVE_NO_ARTIFACT    => ("CF-RESOLVE-001", 422, "No artifact resolves for a requested term and host", "planning/04-api/04-model-resolution.md"),
    RESOLVE_MAJOR_MISMATCH => ("CF-RESOLVE-002", 422, "A job may not mix landscape contract majors", "planning/04-api/04-model-resolution.md"),
    RESOLVE_INADMISSIBLE   => ("CF-RESOLVE-003", 422, "Nearest-host fallback is outside the admissibility threshold", "planning/04-api/04-model-resolution.md"),
    RESOLVE_FALLBACK_USED  => ("CF-RESOLVE-004", 200, "Model resolved from a related host", "planning/04-api/04-model-resolution.md"),
    RESOLVE_DEPRECATED     => ("CF-RESOLVE-005", 200, "Resolved artifact is deprecated", "planning/04-api/04-model-resolution.md"),
    RESOLVE_DIGEST_MISSING => ("CF-RESOLVE-006", 409, "A pinned artifact digest is unavailable; replay cannot proceed", "planning/04-api/04-model-resolution.md"),

    // ---- CF-BUDGET ---------------------------------------------------------
    BUDGET_UNSATISFIABLE   => ("CF-BUDGET-001", 422, "Requested budget cannot be satisfied", "planning/07-decisions/ADR-0008-compute-budget-and-scaling.md"),
    BUDGET_DOWNGRADED      => ("CF-BUDGET-002", 200, "Budget tier downgraded to fit the available hardware", "planning/07-decisions/ADR-0008-compute-budget-and-scaling.md"),
    BUDGET_EXCEEDED        => ("CF-BUDGET-003", 500, "A phase exceeded its measured compute budget", "planning/03-algorithms/08-complexity-and-performance-budgets.md"),

    // ---- CF-LAND: landscape providers --------------------------------------
    LAND_TERM_UNKNOWN      => ("CF-LAND-001", 422, "No provider is registered for a requested term", "planning/02-schemas/landscape_query.schema.json"),
    LAND_CAPABILITY        => ("CF-LAND-002", 422, "Provider does not offer the capability the request needs", "planning/07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md"),
    LAND_DELTA_INEXACT     => ("CF-LAND-003", 500, "An exact delta was required but the provider returned an estimate", "planning/07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md"),
    LAND_BASE_EVICTED      => ("CF-LAND-010", 409, "Delta query base sequence is no longer cached; resend as a full query", "planning/02-schemas/landscape_query.schema.json"),

    // ---- CF-JOINT ----------------------------------------------------------
    JOINT_NO_CONVERGENCE   => ("CF-JOINT-001", 200, "Price decomposition failed to converge; fell back to fully joint search", "planning/03-algorithms/04-systems-level-joint-search.md"),

    // ---- CF-OBS: wet-lab observation ingestion -----------------------------
    OBS_SCHEMA             => ("CF-OBS-001", 422, "Observation failed schema validation", "planning/02-schemas/observation.schema.json"),
    OBS_AS_BUILT_MISMATCH  => ("CF-OBS-002", 422, "As-built digest does not match the referenced candidate", "planning/04-api/03-active-learning-ingestion.md"),
    OBS_NO_CONTROL         => ("CF-OBS-003", 422, "Submission contains no control", "planning/04-api/03-active-learning-ingestion.md"),
    OBS_CENSORING_MISSING  => ("CF-OBS-004", 422, "Assay dynamic range and censoring must be declared", "planning/04-api/03-active-learning-ingestion.md"),
    OBS_MEASURES_MISSING   => ("CF-OBS-005", 422, "Latent measured quantity must be declared", "planning/04-api/03-active-learning-ingestion.md"),
    OBS_EXPLORATION_MISSING=> ("CF-OBS-006", 422, "Exploration-arm assignment must be declared", "planning/04-api/03-active-learning-ingestion.md"),
    OBS_CONDITION_UNKNOWN  => ("CF-OBS-007", 422, "Host or growth condition could not be resolved", "planning/04-api/03-active-learning-ingestion.md"),
    OBS_REPLICATES         => ("CF-OBS-008", 422, "Replicate count is inconsistent with the reported values", "planning/04-api/03-active-learning-ingestion.md"),

    // ---- CF-AUTH / CF-IDEM / CF-QUOTA --------------------------------------
    AUTH_MISSING           => ("CF-AUTH-001", 401, "Missing or malformed credential", "planning/04-api/01-rest-api.md"),
    AUTH_FORBIDDEN         => ("CF-AUTH-002", 403, "Credential is not permitted to perform this action", "planning/04-api/01-rest-api.md"),
    AUTH_FOREIGN_ARTIFACT  => ("CF-AUTH-003", 403, "Artifact belongs to another tenant", "planning/07-decisions/ADR-0006-active-learning-data-isolation.md"),
    IDEM_KEY_CONFLICT      => ("CF-IDEM-001", 409, "Idempotency key was reused with a different body", "planning/04-api/01-rest-api.md"),
    QUOTA_EXCEEDED         => ("CF-QUOTA-001", 429, "Tenant compute quota exceeded", "planning/04-api/01-rest-api.md"),
}

/// Look up a registry entry by wire code.
pub fn lookup(code: &str) -> Option<&'static CfCode> {
    REGISTRY.iter().find(|entry| entry.code == code)
}

/// An error carrying a registered code.
///
/// `detail` is a `Cow` so that a static message costs no allocation. That
/// matters because `CodonSequence::substitute` is on the per-move path and
/// rejects non-synonymous substitutions by returning one of these.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CfError {
    pub code: &'static str,
    pub detail: Cow<'static, str>,
}

impl CfError {
    pub fn new(code: &'static str, detail: impl Into<Cow<'static, str>>) -> Self {
        debug_assert!(
            lookup(code).is_some(),
            "emitted an unregistered error code: {code}"
        );
        Self {
            code,
            detail: detail.into(),
        }
    }

    pub fn entry(&self) -> &'static CfCode {
        lookup(self.code).expect("registered code")
    }

    /// RFC 9457 problem detail, the single error shape across REST and CLI.
    pub fn to_problem_details(&self) -> serde_json::Value {
        let entry = self.entry();
        serde_json::json!({
            "type": format!("https://codonforge.org/problems/{}", entry.code.to_ascii_lowercase()),
            "title": entry.title,
            "status": entry.status,
            "code": entry.code,
            "detail": self.detail.as_ref(),
            "specification": entry.doc,
        })
    }
}

impl fmt::Display for CfError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {} ({})", self.code, self.entry().title, self.detail)
    }
}

impl std::error::Error for CfError {}

pub type Result<T> = std::result::Result<T, CfError>;

#[cfg(test)]
mod tests {
    use super::*;

    /// Codes are unique, well-formed, and carry a plausible status.
    #[test]
    fn unit_registry_is_well_formed() {
        let mut seen = std::collections::HashSet::new();
        for entry in REGISTRY {
            assert!(seen.insert(entry.code), "duplicate code {}", entry.code);
            let parts: Vec<&str> = entry.code.split('-').collect();
            assert_eq!(parts.len(), 3, "malformed code {}", entry.code);
            assert_eq!(parts[0], "CF");
            assert!(parts[1].chars().all(|c| c.is_ascii_uppercase()));
            assert_eq!(parts[2].len(), 3);
            assert!(parts[2].chars().all(|c| c.is_ascii_digit()));
            assert!(
                (100..600).contains(&entry.status),
                "bad status for {}",
                entry.code
            );
            assert!(!entry.title.is_empty());
            assert!(!entry.doc.is_empty());
        }
        assert!(REGISTRY.len() >= 40, "registry looks truncated");
    }

    #[test]
    fn unit_problem_details_shape() {
        let err = CfError::new(code::SEQ_NOT_SYNONYMOUS, "CGT -> AAA at codon 3");
        let pd = err.to_problem_details();
        assert_eq!(pd["code"], "CF-SEQ-004");
        assert_eq!(pd["status"], 422);
        assert!(pd["type"].as_str().unwrap().ends_with("cf-seq-004"));
        assert!(pd["title"].as_str().unwrap().contains("synonymous"));
    }
}
