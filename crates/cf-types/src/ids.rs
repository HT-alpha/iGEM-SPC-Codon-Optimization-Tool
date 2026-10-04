//! Identifier newtypes, validated against the patterns in `common.schema.json`.
//!
//! Ticket P1-08. Contract: `planning/02-schemas/common.schema.json`.
//!
//! The pattern strings below are copied verbatim from the schema, and
//! `contract_id_patterns_match_common_schema` reads the schema at test time and
//! asserts they are still identical. That makes the schema the single source of
//! truth without paying for regex compilation from JSON on every process start,
//! and it means a schema edit that nobody propagated fails CI rather than
//! producing an identifier the API will later reject.

use std::sync::OnceLock;

use regex::Regex;

use crate::error::{code, CfError, Result};

macro_rules! id_newtype {
    (
        $(#[$meta:meta])*
        $name:ident, $schema_def:literal, $pattern:literal
    ) => {
        $(#[$meta])*
        #[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(String);

        impl $name {
            /// The `$defs` key in `common.schema.json` this type mirrors.
            pub const SCHEMA_DEF: &'static str = $schema_def;
            /// The pattern, copied verbatim from the schema.
            pub const PATTERN: &'static str = $pattern;

            fn regex() -> &'static Regex {
                static RE: OnceLock<Regex> = OnceLock::new();
                RE.get_or_init(|| {
                    Regex::new(Self::PATTERN).expect(concat!(
                        "invalid pattern for ", $schema_def
                    ))
                })
            }

            /// Validate and wrap.
            pub fn new(value: impl Into<String>) -> Result<Self> {
                let value = value.into();
                if Self::regex().is_match(&value) {
                    Ok(Self(value))
                } else {
                    Err(CfError::new(
                        code::SEQ_BAD_IDENTIFIER,
                        format!(
                            "{:?} does not match the {} pattern {}",
                            value, $schema_def, $pattern
                        ),
                    ))
                }
            }

            pub fn as_str(&self) -> &str {
                &self.0
            }

            pub fn into_string(self) -> String {
                self.0
            }
        }

        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str(&self.0)
            }
        }

        impl std::fmt::Debug for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "{}({})", $schema_def, self.0)
            }
        }

        impl std::str::FromStr for $name {
            type Err = CfError;
            fn from_str(s: &str) -> Result<Self> {
                Self::new(s)
            }
        }
    };
}

id_newtype!(
    /// Stable host identifier: snake-case binomial, NCBI taxid, optional strain.
    HostId,
    "HostId",
    r"^[a-z0-9]+(?:_[a-z0-9]+)*:[0-9]+(?:/[A-Za-z0-9._-]+)?$"
);

id_newtype!(
    /// Opaque tenant identifier. The isolation invariant of ADR-0006 is stated
    /// over this value, and it is never read from a request body.
    TenantId,
    "TenantId",
    r"^t_[A-Za-z0-9]{8,32}$"
);

id_newtype!(
    /// ULID-based job id, lexicographically sortable by submission time.
    JobId,
    "JobId",
    r"^job_[0-9A-HJKMNP-TV-Z]{26}$"
);

id_newtype!(
    /// Human-readable artifact coordinate.
    ///
    /// Never sufficient on its own for reproducibility: a channel can be
    /// re-pointed, so this must always travel with a [`crate::Digest`]
    /// (ADR-0009).
    ArtifactId,
    "ArtifactId",
    r"^cf:(epistasis_potts|epistasis_dwell|burden_capacity|lm_naturalness|structure_params|kinetics_domainmap|calibrator|host_profile):(?:[a-z0-9_]+:[0-9]+(?:/[A-Za-z0-9._-]+)?|host_agnostic):[0-9]+\.[0-9]+\.[0-9]+$"
);

id_newtype!(
    /// Namespaced objective or diagnostic term identifier.
    ///
    /// Adding a term is additive and bumps no version, which is why
    /// `requested_terms` has no wildcard.
    TermId,
    "TermId",
    r"^(epistasis|structure|kinetics|lm|burden|usage|constraint|meta)\.[a-z0-9_]+$"
);

id_newtype!(
    /// Semantic version.
    SemVer,
    "SemVer",
    r"^[0-9]+\.[0-9]+\.[0-9]+(?:-[0-9A-Za-z.-]+)?(?:\+[0-9A-Za-z.-]+)?$"
);

id_newtype!(
    /// Document schema version (SemVer without pre-release or build metadata).
    SchemaVersion,
    "SchemaVersion",
    r"^[0-9]+\.[0-9]+\.[0-9]+$"
);

impl TermId {
    /// The owning module namespace, e.g. `epistasis` for
    /// `epistasis.dwell_deviation`.
    pub fn namespace(&self) -> &str {
        self.0
            .split('.')
            .next()
            .expect("validated pattern guarantees a dot")
    }
}

/// Every identifier type, paired with its `$defs` key and pattern.
///
/// Drives the contract test; also lets a tool report the whole identifier
/// surface without reflection.
pub static ID_PATTERNS: &[(&str, &str)] = &[
    (HostId::SCHEMA_DEF, HostId::PATTERN),
    (TenantId::SCHEMA_DEF, TenantId::PATTERN),
    (JobId::SCHEMA_DEF, JobId::PATTERN),
    (ArtifactId::SCHEMA_DEF, ArtifactId::PATTERN),
    (TermId::SCHEMA_DEF, TermId::PATTERN),
    (SemVer::SCHEMA_DEF, SemVer::PATTERN),
    (SchemaVersion::SCHEMA_DEF, SchemaVersion::PATTERN),
];

#[cfg(test)]
mod tests {
    use super::*;

    /// contract_*: the Rust patterns are the schema's patterns.
    ///
    /// Reads `schemas/v1/common.schema.json` — the copy the runtime uses — so
    /// this also fails if `P1-03`'s sync between `planning/02-schemas` and
    /// `schemas/v1` has been broken.
    #[test]
    fn contract_id_patterns_match_common_schema() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../schemas/v1/common.schema.json"
        );
        let raw =
            std::fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
        let schema: serde_json::Value = serde_json::from_str(&raw).expect("valid JSON");
        let defs = schema["$defs"].as_object().expect("$defs object");

        for (def_name, pattern) in ID_PATTERNS {
            let schema_pattern = defs
                .get(*def_name)
                .unwrap_or_else(|| panic!("common.schema.json has no $defs/{def_name}"))
                .get("pattern")
                .unwrap_or_else(|| panic!("$defs/{def_name} has no pattern"))
                .as_str()
                .expect("pattern is a string");
            assert_eq!(
                *pattern, schema_pattern,
                "{def_name}: Rust pattern has drifted from common.schema.json"
            );
        }
    }

    #[test]
    fn unit_identifiers_accept_the_schema_examples() {
        assert!(HostId::new("escherichia_coli:562/K-12_MG1655").is_ok());
        assert!(HostId::new("saccharomyces_cerevisiae:4932/S288C").is_ok());
        assert!(HostId::new("homo_sapiens:9606").is_ok());
        assert!(TenantId::new("t_igem2026lab01").is_ok());
        assert!(JobId::new("job_01JBQ8Z4K7M2N3P5R6T7V8W9XY").is_ok());
        assert!(
            ArtifactId::new("cf:epistasis_potts:escherichia_coli:562/K-12_MG1655:1.4.0").is_ok()
        );
        assert!(ArtifactId::new("cf:lm_naturalness:host_agnostic:2.0.1").is_ok());
        for term in [
            "epistasis.dwell_deviation",
            "structure.local_accessibility",
            "burden.growth_cost",
            "usage.cai",
            "meta.confidence",
        ] {
            assert!(
                TermId::new(term).is_ok(),
                "{term} should be a valid term id"
            );
        }
    }

    #[test]
    fn unit_identifiers_reject_near_misses() {
        // Upper case in the binomial.
        assert!(HostId::new("Escherichia_coli:562").is_err());
        // Missing taxid.
        assert!(HostId::new("escherichia_coli").is_err());
        // ULID alphabet excludes I, L, O and U. Y is included (Crockford `V-Z`).
        assert!(JobId::new("job_01JBQ8Z4K7M2N3P5R6T7V8W9XI").is_err());
        // Too short.
        assert!(JobId::new("job_01JBQ8Z4K7M2N3P5R6T7V8W9X").is_err());
        // Unknown namespace.
        assert!(TermId::new("expression.level").is_err());
        // Unknown artifact kind.
        assert!(ArtifactId::new("cf:mystery_model:host_agnostic:1.0.0").is_err());
        // Version must be three components.
        assert!(ArtifactId::new("cf:lm_naturalness:host_agnostic:2.0").is_err());
        assert_eq!(
            TermId::new("nope").unwrap_err().code,
            code::SEQ_BAD_IDENTIFIER
        );
    }

    #[test]
    fn unit_term_namespace() {
        assert_eq!(
            TermId::new("epistasis.dwell_deviation")
                .unwrap()
                .namespace(),
            "epistasis"
        );
        assert_eq!(TermId::new("usage.cai").unwrap().namespace(), "usage");
    }
}
