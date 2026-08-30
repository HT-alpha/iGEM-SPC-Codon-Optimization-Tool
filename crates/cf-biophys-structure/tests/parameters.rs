//! `TST-STRU-005` — parameters are input data, with provenance and a digest.
//!
//! Ticket P2-07a. `DEP-STRU-001`.
//!
//! The sequencing audit for this ticket asked the right question: a schema that
//! validates a well-formed but fabricated parameter set is the same failure the
//! ticket exists to prevent. So these tests deliberately do *not* assert that
//! any energy is plausible — no test here could tell a real Turner value from a
//! well-chosen invention, and pretending otherwise would be the fabrication
//! problem with a test around it.
//!
//! What they do assert is the part a machine can actually check: that a set
//! without provenance is refused, that the digest is over the parameter content
//! rather than the wrapper, and that the digest reaches the lockfile so two runs
//! cannot silently use different thermodynamics.
//!
//! Every fixture below is *obviously* synthetic — round numbers, a citation
//! naming itself as fake. That is deliberate under the no-fabricated-constants
//! rule: a fixture that looked like a real parameter set would be a fabricated
//! parameter set, sitting in the repository waiting to be mistaken for one.

use cf_biophys_structure::{ParameterSource, StructureParameters};
use cf_types::Digest;
use serde_json::{json, Value};

/// A syntactically complete set whose values are transparently invented.
fn fixture() -> Value {
    let mut document = json!({
        "schema_version": "1.0.0",
        "parameter_set_id": "param_synthetic_test_v1",
        "provenance": {
            "citation": "Synthetic values for CodonForge tests; not measured, not usable for design.",
            "doi": null,
            "source_digest": "sha256:0000000000000000000000000000000000000000000000000000000000000000",
            "licence": "CC0-1.0",
            "redistributable": true,
            "transcribed_by": "test fixture"
        },
        "temperature_k": 310.15,
        "stack": { "AU/UA": -1.0, "GC/CG": -2.0, "GU/UG": -0.5 },
        "hairpin_loop": { "3": 5.0, "4": 4.5 },
        "terminal_au_penalty": 0.5
    });
    let digest = Digest::of_document_omitting(&document, &["digest", "created_at"]).unwrap();
    document["digest"] = json!(digest.as_str());
    document
}

#[test]
fn tst_stru_005_a_well_formed_set_with_provenance_is_accepted() {
    let parameters = StructureParameters::from_json(&fixture()).expect("fixture should parse");
    assert_eq!(parameters.parameter_set_id, "param_synthetic_test_v1");
    assert_eq!(parameters.provenance.licence, "CC0-1.0");
    assert_eq!(parameters.stack("GC/CG").unwrap(), -2.0);
}

/// The check the audit asked for: provenance, not well-formedness.
///
/// Each of these documents is *structurally* fine. A schema checking types and
/// required keys on the top level would pass every one of them. They differ
/// only in declining to say where the numbers came from, which is precisely the
/// distinction that matters and precisely what a plausible-looking fabricated
/// set would also fail.
#[test]
fn tst_stru_005_a_set_that_will_not_say_where_it_came_from_is_refused() {
    let cases: [(&str, &str); 5] = [
        ("citation", "a bare identifier is not a reference"),
        ("source_digest", "a transcription must be re-checkable"),
        ("licence", "AMD-009's whole point"),
        ("redistributable", "an unstated claim is not a claim"),
        ("provenance", "no provenance block at all"),
    ];

    for (field, why) in cases {
        let mut document = fixture();
        if field == "provenance" {
            document.as_object_mut().unwrap().remove("provenance");
        } else {
            document["provenance"]
                .as_object_mut()
                .unwrap()
                .remove(field);
        }

        let error = StructureParameters::from_json(&document)
            .expect_err(&format!("a set with no {field} must be refused: {why}"));
        assert_eq!(error.code, "CF-STRU-002", "removing {field}");
    }
}

/// An empty citation is not a citation.
///
/// Worth its own case because "required" in a schema is satisfied by `""`, and
/// the easiest way to get past a provenance requirement without providing
/// provenance is to type nothing into the box.
#[test]
fn tst_stru_005_an_empty_citation_does_not_satisfy_the_requirement() {
    let mut document = fixture();
    document["provenance"]["citation"] = json!("   ");
    let error = StructureParameters::from_json(&document).expect_err("blank citation");
    assert_eq!(error.code, "CF-STRU-002");
}

/// The digest covers the parameters, not just the wrapper.
///
/// Fixture rule: adversarial floats, crossed through text rather than round-
/// tripped in memory. Each perturbation below is the smallest representable
/// change to an energy — a set that differs by one ulp is a different set, and
/// a digest that cannot see the difference is a digest of the wrong thing.
#[test]
fn tst_stru_005_the_digest_is_over_the_parameter_content() {
    let base = fixture();
    let base_digest = StructureParameters::verify_digest(&base).expect("fixture digest verifies");

    let perturbations: [(&str, f64); 3] = [
        ("GC/CG", -2.000_000_000_000_000_4),
        ("AU/UA", -0.999_999_999_999_999_9),
        ("GU/UG", -0.500_000_000_000_000_1),
    ];

    for (pair, value) in perturbations {
        let mut document = base.clone();
        document["stack"][pair] = json!(value);
        // Through text, so the parse path is exercised, not just an in-memory f64.
        let text = serde_json::to_string(&document).unwrap();
        let reparsed: Value = serde_json::from_str(&text).unwrap();
        let moved = Digest::of_document_omitting(&reparsed, &["digest", "created_at"]).unwrap();
        assert_ne!(
            moved.as_str(),
            base_digest.as_str(),
            "a one-ulp change to {pair} must change the digest; if it does not, two runs \
             with measurably different thermodynamics share a lockfile entry"
        );
    }
}

/// A document whose stated digest disagrees with its content is refused.
#[test]
fn tst_stru_005_an_edited_set_does_not_keep_its_old_digest() {
    let mut document = fixture();
    document["temperature_k"] = json!(298.15);
    let error = StructureParameters::verify_digest(&document)
        .expect_err("an edited document must not verify against its old digest");
    assert_eq!(error.code, "CF-STRU-002");
}

/// No parameters means a refusal, not a default.
#[test]
fn tst_stru_005_a_run_without_parameters_is_refused_not_defaulted() {
    let empty = ParameterSource::none();
    let error = empty.require().expect_err("a build with no table must refuse");
    assert_eq!(error.code, "CF-STRU-001");
    assert!(
        empty.lockfile_digest().is_none(),
        "nothing is pinned when nothing was supplied"
    );

    let supplied = ParameterSource::from(StructureParameters::from_json(&fixture()).unwrap());
    assert!(
        supplied.require().is_ok(),
        "a supplied set is usable, so the refusal above is about absence, not a \
         permanently-closed door"
    );
    assert!(
        supplied.lockfile_digest().is_some(),
        "a supplied set reaches the lockfile"
    );
}

/// A missing stack is refused by name, not silently zeroed.
#[test]
fn tst_stru_005_an_absent_stack_is_named_not_zeroed() {
    let parameters = StructureParameters::from_json(&fixture()).unwrap();
    let error = parameters
        .stack("CG/GC")
        .expect_err("an absent pair must be refused");
    assert_eq!(error.code, "CF-STRU-002");
    assert!(
        error.to_string().contains("CG/GC"),
        "the refusal must name the pair the user has to go and find: {error}"
    );
}

/// The repository ships no parameter table.
///
/// This is the test that keeps `DEP-STRU-001` from quietly eroding. The risk is
/// not that someone writes `const TURNER: ...` in a fit of malice; it is that a
/// later ticket needs a number to make a test pass and types a plausible one in.
#[test]
fn tst_stru_005_no_parameter_table_is_compiled_in() {
    let source = include_str!("../src/parameters.rs");
    let lib = include_str!("../src/lib.rs");

    for (name, text) in [("parameters.rs", source), ("lib.rs", lib)] {
        assert!(
            !text.contains("const STACK") && !text.contains("static STACK"),
            "{name} appears to have grown a compiled-in stacking table; parameters are \
             input data per DEP-STRU-001"
        );
    }

    // The reader holds no energies of its own: the only floats it may name are
    // structural defaults, not thermodynamics.
    let energy_like = source
        .lines()
        .filter(|line| !line.trim_start().starts_with("//"))
        .filter(|line| line.contains("-1.") || line.contains("-2.") || line.contains("kcal"))
        .count();
    assert_eq!(
        energy_like, 0,
        "the parameter reader should contain no energies; it reads them"
    );
}
