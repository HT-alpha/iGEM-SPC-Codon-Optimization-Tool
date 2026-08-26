"""Cross-runtime contract tests for the canonical sequence representation.

Tickets P1-06, P1-07. Test ids: TST-LAPI-001 (round-trip portion),
TST-CORE-003 (cross-implementation portion).

These are the tests that keep the two runtimes honest. `cf-types` and the JSON
Schemas are separate artifacts maintained by hand, and ADR-0001 accepts that
duplication as the price of the two-language split — on the condition that
drift fails CI rather than surfacing at runtime. That condition is what these
tests implement:

  * a document produced by Rust validates against the schema Python reads;
  * a digest computed by Rust matches one computed by an independent
    implementation of RFC 8785 written here, in Python, from the spec.

The second one matters more than it looks. A canonicalizer that is
self-consistent is trivially "correct" against itself; the only way to find out
whether it implements the *spec* is to compare it against a different
implementation.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator
from referencing import Registry, Resource

from codonforge import _core

ROOT = Path(__file__).resolve().parents[2]
SCHEMAS = ROOT / "schemas" / "v1"
EXAMPLES = ROOT / "planning" / "02-schemas" / "examples"

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


def test_contract_planning_example_parses(sequence_validator):
    """The worked example in the planning package is structurally real.

    The example's `digest` is a hand-written placeholder (as are the digests in
    every worked example — they exist to show the field's shape, and several of
    them are cross-referenced between examples, so they cannot all be genuine).
    It is therefore stripped before the invariant check; everything else — frame,
    translation, lock ranges, flanks — is checked against the implementation.
    """
    example = json.loads((EXAMPLES / "sequence.example.json").read_text(encoding="utf-8"))
    assert not list(sequence_validator.iter_errors(example)), "example must match its schema"

    without_digest = {k: v for k, v in example.items() if k != "digest"}
    round_tripped = json.loads(_core.sequence_round_trip(json.dumps(without_digest)))

    assert round_tripped["cds"] == example["cds"]
    assert round_tripped["protein"] == example["protein"]
    assert round_tripped["codon_count"] == example["codon_count"]
    assert round_tripped["five_prime_flank"] == example["five_prime_flank"]
    assert round_tripped["locked_codons"] == example["locked_codons"]
    # And the recomputed digest is a real one.
    assert round_tripped["digest"].startswith("sha256:")
    assert len(round_tripped["digest"]) == len("sha256:") + 64


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
# an independent RFC 8785 implementation, for cross-checking
# --------------------------------------------------------------------------- #


def _es_number(value: float) -> str:
    """ECMAScript Number::toString, written from ECMA-262 6.1.6.1.20.

    Deliberately a second implementation rather than a port of the Rust one:
    the point is to disagree if either is wrong.
    """
    if value != value or value in (float("inf"), float("-inf")):
        raise ValueError("not representable")
    if value == 0:
        return "0"
    if value < 0:
        return "-" + _es_number(-value)

    mantissa, _, exponent = repr(value).partition("e")
    if exponent:
        digits = mantissa.replace(".", "").rstrip("0") or "0"
        n = int(exponent) + 1
    else:
        # repr gives the shortest round-tripping decimal for a Python float.
        intpart, _, frac = mantissa.partition(".")
        frac = frac.rstrip("0")
        if intpart == "0":
            stripped = frac.lstrip("0")
            n = -(len(frac) - len(stripped))
            digits = stripped
        else:
            digits = (intpart + frac).rstrip("0") or "0"
            n = len(intpart)
    digits = digits.rstrip("0") or "0"
    k = len(digits)

    if k <= n <= 21:
        return digits + "0" * (n - k)
    if 0 < n <= 21:
        return digits[:n] + "." + digits[n:]
    if -6 < n <= 0:
        return "0." + "0" * (-n) + digits
    e = n - 1
    sign = "+" if e >= 0 else "-"
    head = digits if k == 1 else digits[0] + "." + digits[1:]
    return f"{head}e{sign}{abs(e)}"


def _canonicalize(value) -> str:
    if value is None:
        return "null"
    if value is True:
        return "true"
    if value is False:
        return "false"
    if isinstance(value, int):
        if abs(value) > 9007199254740991:
            raise ValueError("outside the IEEE-754 safe integer range")
        return str(value)
    if isinstance(value, float):
        return _es_number(value)
    if isinstance(value, str):
        return _escape(value)
    if isinstance(value, list):
        return "[" + ",".join(_canonicalize(v) for v in value) + "]"
    if isinstance(value, dict):
        keys = sorted(value, key=lambda k: k.encode("utf-16-be"))
        return "{" + ",".join(f"{_escape(k)}:{_canonicalize(value[k])}" for k in keys) + "}"
    raise TypeError(type(value))


def _escape(text: str) -> str:
    out = ['"']
    short = {
        '"': '\\"',
        "\\": "\\\\",
        "\b": "\\b",
        "\f": "\\f",
        "\n": "\\n",
        "\r": "\\r",
        "\t": "\\t",
    }
    for ch in text:
        if ch in short:
            out.append(short[ch])
        elif ord(ch) < 0x20:
            out.append(f"\\u{ord(ch):04x}")
        else:
            out.append(ch)
    out.append('"')
    return "".join(out)


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
    python = _canonicalize(json.loads(payload))
    assert rust == python, f"canonical forms differ:\n  rust  : {rust}\n  python: {python}"

    # And therefore so must the digests.
    expected = "sha256:" + hashlib.sha256(python.encode("utf-8")).hexdigest()
    assert _core.digest_json(payload) == expected


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
