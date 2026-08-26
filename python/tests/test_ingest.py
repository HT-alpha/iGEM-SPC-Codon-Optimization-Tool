"""Ingest tests: `TST-INGE-001`, `TST-INGE-002`, `TST-INGE-003`.

Tickets P1-09..P1-12.

The property these tests exist to pin down is the one the round's brief calls
out: a sequence that would fail `CodonSequence`'s invariants must fail **at
ingest**, with a typed code, rather than downstream in a ticket that assumed
ingest had already checked.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator
from referencing import Registry, Resource

from codonforge import _core
from codonforge.ingest import (
    IngestError,
    IngestReport,
    build_host_profile,
    codon_usage_from_cds,
    import_seeds,
    parse_fasta,
    read_cds,
    read_construct,
)

ROOT = Path(__file__).resolve().parents[2]
SCHEMAS = ROOT / "schemas" / "v1"

DEMO_CDS = "ATGAAACGTATTCTGTAA"
DEMO_PROTEIN = "MKRIL*"


@pytest.fixture(scope="module")
def registry() -> Registry:
    resources = []
    for path in sorted(SCHEMAS.glob("*.schema.json")):
        doc = json.loads(path.read_text(encoding="utf-8"))
        if "$id" in doc:
            resources.append((doc["$id"], Resource.from_contents(doc)))
    return Registry().with_resources(resources)


def validator_for(name: str, registry: Registry) -> Draft202012Validator:
    schema = json.loads((SCHEMAS / name).read_text(encoding="utf-8"))
    return Draft202012Validator(schema, registry=registry)


# --------------------------------------------------------------------------- #
# P1-09 — FASTA
# --------------------------------------------------------------------------- #


def test_tst_inge_001_reads_a_cds_and_validates_against_the_schema(registry):
    report = IngestReport()
    document = read_cds(
        f">demo test record\n{DEMO_CDS}\n",
        translation_table=11,
        uri="demo.fa",
        report=report,
    )
    errors = list(validator_for("sequence.schema.json", registry).iter_errors(document))
    assert not errors, errors
    assert document["protein"] == DEMO_PROTEIN
    assert document["source"]["kind"] == "fasta"
    assert document["source"]["record_id"] == "demo"
    # The digest covers `source`, so it must have been recomputed after adding it.
    recomputed = _core.digest_json(json.dumps({k: v for k, v in document.items() if k != "digest"}))
    assert document["digest"] == recomputed


def test_tst_inge_001_rna_and_case_are_transcoded_and_reported():
    report = IngestReport()
    document = read_cds(">demo\naugaaacguauucuguaa\n", translation_table=11, report=report)
    assert document["cds"] == DEMO_CDS
    # Silently altering the user's input without saying so is the failure mode.
    assert any("U->T" in note for note in report.notes), report.notes
    assert any("upper-cased" in note for note in report.notes), report.notes


def test_tst_inge_001_line_wrapping_and_multiple_records():
    records = parse_fasta(">a d1\nATGAAA\nCGTATT\n>b d2\nCTGTAA\n")
    assert [r.identifier for r in records] == ["a", "b"]
    assert records[0].sequence == "ATGAAACGTATT"
    assert records[0].description == "d1"


@pytest.mark.parametrize(
    ("text", "code", "needle"),
    [
        (">\nATGAAATAA\n", "CF-SEQ-001", "header is empty"),
        ("ATGAAATAA\n", "CF-SEQ-001", "before any '>'"),
        ("", "CF-SEQ-001", "no FASTA records"),
        (">a\nATGAAATAA\n>a\nATGAAATAA\n", "CF-SEQ-013", "duplicate"),
        (">a\nATG-AA-TAA\n", "CF-SEQ-002", "alignment gap"),
        (">a\nATGNNNTAA\n", "CF-SEQ-002", "ambiguity code"),
        (">a\nATGZZZTAA\n", "CF-SEQ-002", "not nucleotides"),
    ],
)
def test_tst_inge_001_malformed_input_is_rejected_with_a_typed_code(text, code, needle):
    with pytest.raises(IngestError) as excinfo:
        read_cds(text, translation_table=11)
    assert excinfo.value.code == code
    assert needle in excinfo.value.detail


def test_tst_inge_001_core_invariants_surface_at_ingest_not_downstream():
    """A sequence that would fail `CodonSequence` fails here, with its code."""
    cases = {
        ">a\nATGAAAC\n": "CF-SEQ-001",  # not a whole number of codons
        ">a\nATGTAAAAATAA\n": "CF-SEQ-003",  # internal stop
        ">a\n\n": "CF-SEQ-012",  # empty record body -> empty sequence
    }
    for text, code in cases.items():
        with pytest.raises(IngestError) as excinfo:
            read_cds(text, translation_table=11)
        assert excinfo.value.code == code, f"{text!r} -> {excinfo.value}"

    # An undefined translation table is the core's error too, not re-derived.
    with pytest.raises(IngestError) as excinfo:
        read_cds(f">a\n{DEMO_CDS}\n", translation_table=7)
    assert excinfo.value.code == "CF-SEQ-007"


def test_ingest_error_carries_location_and_problem_details():
    with pytest.raises(IngestError) as excinfo:
        read_cds(">a\nATGNNNTAA\n", translation_table=11, uri="in.fa")
    problem = excinfo.value.to_problem_details()
    assert problem["code"] == "CF-SEQ-002"
    assert problem["status"] == 422
    assert "in.fa" in problem["location"]
    assert "record 'a'" in problem["location"]


# --------------------------------------------------------------------------- #
# P1-10 — GenBank and construct context
# --------------------------------------------------------------------------- #

GENBANK = """LOCUS       demo                      60 bp    DNA     circular SYN 01-JAN-2026
DEFINITION  demo construct.
FEATURES             Location/Qualifiers
     source          1..60
     RBS             1..17
                     /label="rbs1"
     CDS             18..35
                     /label="reporter"
     terminator      36..60
                     /label="term"
ORIGIN
        1 aaagaggaga aagatctatg aaacgtattc tgtaaccagg catcaaataa aacgaaaggc
//
"""


def test_p1_10_reads_a_construct_and_validates(registry):
    report = IngestReport()
    construct = read_construct(
        GENBANK,
        host_id="escherichia_coli:562/K-12_MG1655",
        translation_table=11,
        report=report,
    )
    errors = list(validator_for("construct_context.schema.json", registry).iter_errors(construct))
    assert not errors, errors
    assert construct["topology"] == "circular"
    kinds = {e["kind"] for e in construct["elements"]}
    assert {"rbs", "cds", "terminator"} <= kinds
    assert len(construct["cistrons"]) == 1
    # The nearest upstream RBS is paired, not guessed at.
    cistron = construct["cistrons"][0]
    rbs = next(e for e in construct["elements"] if e["kind"] == "rbs")
    assert cistron["rbs_element_id"] == rbs["element_id"]
    # Only the CDS is editable.
    assert [e["editable"] for e in construct["elements"]] == [
        e["kind"] == "cds" for e in construct["elements"]
    ]


def test_p1_10_inconsistent_cds_span_is_rejected_at_ingest():
    """The one place a frame error is possible checks for it."""
    broken = GENBANK.replace("CDS             18..35", "CDS             18..36")
    with pytest.raises(IngestError) as excinfo:
        read_construct(broken, host_id="escherichia_coli:562", translation_table=11)
    assert excinfo.value.code == "CF-SEQ-001"
    assert "whole number of codons" in excinfo.value.detail


def test_p1_10_feature_running_off_the_record_is_rejected():
    broken = GENBANK.replace("terminator      36..60", "terminator      36..600")
    with pytest.raises(IngestError) as excinfo:
        read_construct(broken, host_id="escherichia_coli:562", translation_table=11)
    assert excinfo.value.code == "CF-SEQ-005"


def test_p1_10_missing_rbs_is_reported_not_invented():
    without_rbs = GENBANK.replace(
        '     RBS             1..17\n                     /label="rbs1"\n', ""
    )
    report = IngestReport()
    construct = read_construct(
        without_rbs,
        host_id="escherichia_coli:562",
        translation_table=11,
        report=report,
    )
    assert construct["cistrons"][0]["rbs_element_id"] is None
    assert any("no upstream RBS" in note for note in report.notes), report.notes


def test_p1_10_record_without_a_cds_is_rejected():
    no_cds = GENBANK.replace(
        '     CDS             18..35\n                     /label="reporter"\n', ""
    )
    with pytest.raises(IngestError) as excinfo:
        read_construct(no_cds, host_id="escherichia_coli:562", translation_table=11)
    assert "no CDS feature" in excinfo.value.detail


# --------------------------------------------------------------------------- #
# P1-11 — host profile
# --------------------------------------------------------------------------- #


def test_tst_inge_002_host_profile_validates_and_declares_its_tier(registry):
    usage = codon_usage_from_cds(
        [DEMO_CDS, "ATGAAGCGCATCCTGTAA", "ATGAAACGTATCCTGTAA"], translation_table=11
    )
    report = IngestReport()
    profile = build_host_profile(
        host_id="escherichia_coli:562/K-12_MG1655",
        taxid=511145,
        translation_table=11,
        codon_usage=usage,
        trna_gene_copy_number={"TTT": 6, "ACG": 4, "CAT": 8},
        elongation_rate_codons_per_s=16.5,
        report=report,
    )
    errors = list(validator_for("host_profile.schema.json", registry).iter_errors(profile))
    assert not errors, errors

    # ADR-0004: a genome-only host is capacity tier 3, and its unit says the
    # burden term is a mismatch score rather than a growth-rate prediction.
    assert profile["capacity"]["estimator_tier"] == 3
    assert profile["capacity"]["theta_unit"] == "fraction_of_native_demand"
    assert profile["capacity"]["uncertainty"]["calibrated"] is False

    # No tAI s-values were supplied, so no weights were invented and the burden
    # term is not offered.
    assert "tai" not in profile
    assert "burden.growth_cost" not in profile["data_completeness"]["supported_terms"]
    assert any("no tAI s-values" in w for w in report.warnings), report.warnings


def test_tst_inge_002_elongation_rate_has_no_default():
    """A measured physiological parameter is never inherited from another host."""
    usage = codon_usage_from_cds([DEMO_CDS], translation_table=11)
    with pytest.raises(TypeError):
        build_host_profile(  # type: ignore[call-arg]
            host_id="x:1",
            taxid=1,
            translation_table=11,
            codon_usage=usage,
            trna_gene_copy_number={"TTT": 1},
        )
    with pytest.raises(IngestError) as excinfo:
        build_host_profile(
            host_id="x:1",
            taxid=1,
            translation_table=11,
            codon_usage=usage,
            trna_gene_copy_number={"TTT": 1},
            elongation_rate_codons_per_s=0.0,
        )
    assert excinfo.value.code == "CF-HOST-002"
    assert "no default" in excinfo.value.detail


def test_tst_inge_002_relative_adaptiveness_is_derived_correctly():
    # Arginine: two CGT and one CGC, so w(CGT) = 1.0 and w(CGC) = 0.5.
    usage = codon_usage_from_cds(["ATGCGTTAA", "ATGCGTTAA", "ATGCGCTAA"], translation_table=11)
    table = usage["table"]
    assert table["CGT"]["count"] == 2
    assert table["CGC"]["count"] == 1
    assert table["CGT"]["relative_adaptiveness"] == pytest.approx(1.0)
    assert table["CGC"]["relative_adaptiveness"] == pytest.approx(0.5)


def test_tst_inge_002_out_of_frame_cds_is_rejected():
    with pytest.raises(IngestError) as excinfo:
        codon_usage_from_cds(["ATGAAAC"], translation_table=11)
    assert excinfo.value.code == "CF-SEQ-001"


# --------------------------------------------------------------------------- #
# P1-12 — seed importers
# --------------------------------------------------------------------------- #


def test_tst_inge_003_accepts_matching_seeds_and_validates(registry):
    seed_set = import_seeds(
        DEMO_PROTEIN,
        11,
        [
            {
                "tool": "lineardesign",
                "payload": ">ld\nATGAAACGCATCCTGTAA\n",
                "parameters": {"lambda": 0.0},
            },
            {"tool": "codontransformer", "payload": "ATGAAACGTATTTTGTAA"},
        ],
    )
    errors = list(validator_for("seed_set.schema.json", registry).iter_errors(seed_set))
    assert not errors, errors
    assert all(s["acceptance"]["accepted"] for s in seed_set["seeds"]), seed_set["seeds"]
    assert seed_set["seeds"][0]["source"]["parameters"] == {"lambda": 0.0}


def test_tst_inge_003_protein_mismatch_reports_a_positional_diff():
    """The most likely user error, diagnosed rather than merely refused."""
    report = IngestReport()
    seed_set = import_seeds(
        DEMO_PROTEIN,
        11,
        [{"tool": "manual", "payload": "ATGAAACATATTCTGTAA"}],  # R -> H at position 2
        report=report,
    )
    seed = seed_set["seeds"][0]
    assert seed["acceptance"]["accepted"] is False
    assert seed["acceptance"]["reason"] == "protein_mismatch"
    assert seed["acceptance"]["protein_diff"] == [{"position": 2, "expected": "R", "found": "H"}]
    assert any("position(s) [2]" in w for w in report.warnings), report.warnings


def test_tst_inge_003_length_mismatch_is_distinguished_from_substitution():
    seed_set = import_seeds(DEMO_PROTEIN, 11, [{"tool": "manual", "payload": "ATGAAACGTTAA"}])
    seed = seed_set["seeds"][0]
    assert seed["acceptance"]["reason"] == "length_mismatch"
    # The diff shows where it ran out, using '-' for the missing residues.
    assert seed["acceptance"]["protein_diff"][0]["found"] in {"*", "-"}


def test_tst_inge_003_one_bad_seed_does_not_lose_the_others():
    seed_set = import_seeds(
        DEMO_PROTEIN,
        11,
        [
            {"tool": "manual", "payload": "ATGAAACGCATCCTGTAA"},
            {"tool": "manual", "payload": "ATGAAACATATTCTGTAA"},
            {"tool": "manual", "payload": "ATGAAACGTATCCTGTAA"},
        ],
    )
    accepted = [s for s in seed_set["seeds"] if s["acceptance"]["accepted"]]
    assert len(accepted) == 2, "a rejected seed must not fail the batch"


def test_tst_inge_003_duplicates_and_bad_input_are_classified():
    seed_set = import_seeds(
        DEMO_PROTEIN,
        11,
        [
            {"tool": "manual", "payload": DEMO_CDS},
            {"tool": "manual", "payload": DEMO_CDS},
            {"tool": "manual", "payload": "ATGAAACG"},
            {"tool": "manual", "payload": "ATGNNNCGTATTCTGTAA"},
        ],
    )
    reasons = [s["acceptance"]["reason"] for s in seed_set["seeds"]]
    assert reasons[0] == "ok"
    assert reasons[1] == "duplicate"
    assert reasons[2] == "not_whole_codons"
    assert reasons[3] == "invalid_alphabet"


def test_tst_inge_003_unknown_tool_is_refused():
    seed_set = import_seeds(DEMO_PROTEIN, 11, [{"tool": "not_a_real_tool", "payload": DEMO_CDS}])
    assert seed_set["seeds"][0]["acceptance"]["accepted"] is False


def test_tst_inge_003_target_protein_digest_matches_the_documented_convention():
    seed_set = import_seeds(DEMO_PROTEIN, 11, [{"tool": "manual", "payload": DEMO_CDS}])
    # seed_set.schema.json: sha256 over the UTF-8 bytes of the protein, via the
    # canonical JSON of the string.
    assert seed_set["target_protein_digest"] == _core.digest_json(json.dumps(DEMO_PROTEIN))
    # And a seed's digest is the CodonSequence DOCUMENT digest, not the cds.
    document = json.loads(_core.sequence_document(DEMO_CDS, 11, None, None, None))
    expected = document["digest"]
    assert seed_set["seeds"][0]["digest"] == expected
