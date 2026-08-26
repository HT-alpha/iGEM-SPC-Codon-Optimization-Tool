"""Cross-runtime contract tests for the canonical sequence representation.

Tickets P1-06, P1-07. Test ids: TST-LAPI-001 (round-trip portion),
TST-CORE-003 (cross-implementation portion).

These are the tests that keep the two runtimes honest. `cf-types` and the JSON
Schemas are separate artifacts maintained by hand, and ADR-0001 accepts that
duplication as the price of the two-language split — on the condition that
drift fails CI rather than surfacing at runtime. That condition is what these
tests implement:

  * a document produced by Rust validates against the schema Python reads;
  * a digest computed by Rust matches one computed by `planning/tools/jcs.py`,
    an independent implementation of RFC 8785 written from the spec.

The second one matters more than it looks. A canonicalizer that is
self-consistent is trivially "correct" against itself; the only way to find out
whether it implements the *spec* is to compare it against a different
implementation. `jcs.py` is also what computes the digests stored in the worked
examples, so this test is what ties those fixtures to the implementation.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator
from referencing import Registry, Resource

from codonforge import _core

ROOT = Path(__file__).resolve().parents[2]
SCHEMAS = ROOT / "schemas" / "v1"
EXAMPLES = ROOT / "planning" / "02-schemas" / "examples"

sys.path.insert(0, str(ROOT / "planning" / "tools"))
import jcs  # noqa: E402

DEMO_CDS = "ATGAAACGTATTCTGTAA"
DEMO_FLANK5 = "AAAGAGGAGAAAGATCT"


# --------------------------------------------------------------------------- #
# schema validation
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def registry() -> Registry:
    resources = []
    for path in sorted(SCHEMAS.glob("*.schema.json")):
        doc = json.loads(path.read_text(encoding="utf-8"))
        if "$id" in doc:
            resources.append((doc["$id"], Resource.from_contents(doc)))
    return Registry().with_resources(resources)


@pytest.fixture(scope="module")
def sequence_validator(registry: Registry) -> Draft202012Validator:
    schema = json.loads((SCHEMAS / "sequence.schema.json").read_text(encoding="utf-8"))
    return Draft202012Validator(schema, registry=registry)


def test_contract_rust_document_validates_against_the_schema(sequence_validator):
    """A document the Rust core produces satisfies the schema Python reads."""
    raw = _core.sequence_document(DEMO_CDS, 11, DEMO_FLANK5, None, "seq_demo")
    doc = json.loads(raw)
    errors = sorted(sequence_validator.iter_errors(doc), key=lambda e: list(e.path))
    assert not errors, "\n".join(
        f"{'/'.join(str(p) for p in e.absolute_path) or '<root>'}: {e.message}" for e in errors
    )
    assert doc["protein"] == "MKRIL*"
    assert doc["codon_count"] == 6
    assert len(doc["cds"]) == 3 * doc["codon_count"]
    assert doc["alphabet"] == "DNA"
    reasons = {lock["reason"] for lock in doc["locked_codons"]}
    assert reasons == {"start_codon", "stop_codon"}


def test_contract_document_round_trips_through_rust():
    """Rust -> JSON -> Rust -> JSON is a fixed point."""
    first = _core.sequence_document(DEMO_CDS, 11, DEMO_FLANK5, None, "seq_demo")
    second = _core.sequence_round_trip(first)
    assert json.loads(first) == json.loads(second)


def test_contract_planning_example_parses_with_its_stored_digest(sequence_validator):
    """The worked example is accepted whole, digest included.

    This previously had to strip `digest` before handing the example to the
    core, because the examples carried hand-written placeholders. They no longer
    do (see planning/02-schemas/examples/digest_map.yaml), so the document is
    passed intact and invariant I5 — the stored digest equals a fresh
    canonicalize-then-hash of the payload — is checked by the implementation
    rather than skipped.

    Nothing here strips a digest for any other reason.
    """
    example = json.loads((EXAMPLES / "sequence.example.json").read_text(encoding="utf-8"))
    assert not list(sequence_validator.iter_errors(example)), "example must match its schema"

    round_tripped = json.loads(_core.sequence_round_trip(json.dumps(example)))

    assert round_tripped == example, (
        "the example must survive Rust parse -> re-emit unchanged; a difference "
        "means the stored digest or some derived field disagrees with the core"
    )


def test_contract_invariants_are_enforced_across_the_boundary():
    """The CF-* code survives the crossing, so Python can branch on it."""
    doc = json.loads(_core.sequence_document(DEMO_CDS, 11, None, None, None))

    doc["protein"] = "MKRIQ*"
    with pytest.raises(ValueError) as excinfo:
        _core.sequence_round_trip(json.dumps(doc))
    problem = json.loads(str(excinfo.value))
    assert problem["code"] == "CF-SEQ-008"
    assert problem["status"] == 422

    fresh = json.loads(_core.sequence_document(DEMO_CDS, 11, None, None, None))
    fresh["digest"] = "sha256:" + "0" * 64
    with pytest.raises(ValueError) as excinfo:
        _core.sequence_round_trip(json.dumps(fresh))
    assert json.loads(str(excinfo.value))["code"] == "CF-SEQ-010"


# --------------------------------------------------------------------------- #
# Rust against the from-spec Python reference
# --------------------------------------------------------------------------- #

CANONICAL_CASES = [
    {"b": 1, "a": 2, "c": {"z": 1, "y": 2}},
    {"nested": [1, [2, [3, {"k": "v"}]]], "empty_obj": {}, "empty_arr": []},
    {"floats": [0.1, 1.5, 100.0, 1e20, 1e21, 1e22, 1e-6, 1e-7, 1.2345e-7]},
    {"edge": [0.0, -0.0, -1.0, 5e-324, 1.7976931348623157e308, 333333333.33333329]},
    {"strings": ['a"b\\c', "tab\there", "nl\nhere", "\u0000\u001f", "é☃𝄞"]},
    {"bools": [True, False, None], "int": 42, "neg": -42},
    {"\u1d11e-ish": 1, "\uff3a": 2, "a": 3},
]


@pytest.mark.parametrize("case", CANONICAL_CASES, ids=range(len(CANONICAL_CASES)))
def test_contract_tst_core_003_rust_and_python_canonicalize_identically(case):
    """Two independent RFC 8785 implementations must agree, byte for byte."""
    payload = json.dumps(case)
    rust = _core.canonical_json(payload)
    python = jcs.canonicalize(json.loads(payload))
    assert rust == python, f"canonical forms differ:\n  rust  : {rust}\n  python: {python}"

    expected = "sha256:" + hashlib.sha256(python.encode("utf-8")).hexdigest()
    assert _core.digest_json(payload) == expected
    assert jcs.digest(json.loads(payload)) == expected


def test_contract_codon_sequence_document_agrees_between_runtimes():
    """The document form used for the example fixtures is the core's own.

    `jcs.codon_sequence_document` is what computes the sequence digests stored
    in the Pareto-front and seed-set examples. If it drifted from
    `CodonSequence::from_cds`, those fixtures would be self-consistent and
    wrong — which is precisely the failure this round set out to remove — so the
    two constructions are compared field by field rather than only by digest.
    """
    tables = jcs.load_genetic_codes(ROOT)
    for cds, table in [
        ("ATGAAACGTATTCTGTAA", 11),
        ("ATGAAACGCATCCTGTAA", 11),
        ("ATGGCAGAACTGAAATAA", 11),
        ("ATGAAACGT", 11),  # no terminal stop: only the start codon is locked
        ("ATGAAACGTATTCTGTAA", 1),
    ]:
        python_doc = jcs.codon_sequence_document(cds, table, tables)
        rust_doc = json.loads(_core.sequence_document(cds, table, None, None, None))
        assert python_doc == rust_doc, f"construction differs for {cds} under table {table}"


def test_contract_digest_omits_the_named_field():
    payload = json.dumps({"a": 1, "digest": "sha256:" + "0" * 64})
    stripped = json.dumps({"a": 1})
    assert _core.digest_json(payload, "digest") == _core.digest_json(stripped)


def test_contract_defined_translation_tables_match_the_vendored_data():
    """The 27-table set, and the six retired ids, are visible from Python."""
    defined = _core.defined_translation_tables()
    assert len(defined) == 27
    assert set(range(1, 34)) - set(defined) == {7, 8, 17, 18, 19, 20}
    with pytest.raises(ValueError) as excinfo:
        _core.translate(DEMO_CDS, 7)
    assert json.loads(str(excinfo.value))["code"] == "CF-SEQ-007"
