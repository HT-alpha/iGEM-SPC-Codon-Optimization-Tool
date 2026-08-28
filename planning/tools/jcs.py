"""Reference implementation of RFC 8785 (JCS) and the CodonSequence document form.

This is a *second* implementation, written from the specifications, of two
things the Rust core also implements:

  * `canonicalize` / `digest`     — RFC 8785 canonical JSON and its sha256
  * `codon_sequence_document`     — the document form of
                                    planning/02-schemas/sequence.schema.json

The duplication is deliberate and is the point. A canonicalizer that agrees
with itself proves nothing; the only way to find out whether an implementation
follows the spec rather than its own habits is to compare it against a
different one. `python/tests/test_contract_sequence_schema.py` asserts that
this module and `codonforge._core` agree byte for byte, so a divergence in
either fails CI.

It is kept dependency-free (stdlib only) so that the planning-package validator
can run in a CI job that has no Rust toolchain and no compiled core.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

__all__ = [
    "canonicalize",
    "digest",
    "digest_of_document",
    "es_number",
    "codon_sequence_document",
    "load_genetic_codes",
]

# Largest integer exactly representable as an IEEE-754 double. RFC 8785 is
# defined over doubles, so an integer beyond this cannot round-trip and must be
# refused rather than silently rendered.
MAX_SAFE_INTEGER = 9_007_199_254_740_991

_SHORT_ESCAPES = {
    '"': '\\"',
    "\\": "\\\\",
    "\b": "\\b",
    "\f": "\\f",
    "\n": "\\n",
    "\r": "\\r",
    "\t": "\\t",
}


def es_number(value: float) -> str:
    """ECMAScript ``Number::toString`` (ECMA-262 6.1.6.1.20).

    RFC 8785 requires this rendering, which is *not* what Python's ``repr`` or
    Rust's ``Display`` produce once exponents are involved: ``1e21`` is
    ``"1e+21"`` here and ``"1000000000000000000000"`` in both languages'
    default formatting.
    """
    if value != value or value in (float("inf"), float("-inf")):
        raise ValueError("NaN and infinity have no JSON representation")
    if value == 0:
        return "0"  # also covers -0.0, which ECMAScript renders as "0"
    if value < 0:
        return "-" + es_number(-value)

    mantissa, _, exponent = repr(value).partition("e")
    if exponent:
        digits = mantissa.replace(".", "").rstrip("0") or "0"
        n = int(exponent) + 1
    else:
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


def _escape(text: str) -> str:
    out = ['"']
    for ch in text:
        if ch in _SHORT_ESCAPES:
            out.append(_SHORT_ESCAPES[ch])
        elif ord(ch) < 0x20:
            out.append(f"\\u{ord(ch):04x}")
        else:
            out.append(ch)
    out.append('"')
    return "".join(out)


def canonicalize(value) -> str:
    """Serialize to the RFC 8785 canonical form."""
    if value is None:
        return "null"
    if value is True:
        return "true"
    if value is False:
        return "false"
    if isinstance(value, int):
        if abs(value) > MAX_SAFE_INTEGER:
            raise ValueError(
                f"integer {value} exceeds the IEEE-754 safe integer range required by RFC 8785"
            )
        return str(value)
    if isinstance(value, float):
        return es_number(value)
    if isinstance(value, str):
        return _escape(value)
    if isinstance(value, list):
        return "[" + ",".join(canonicalize(v) for v in value) + "]"
    if isinstance(value, dict):
        # Sorted by UTF-16 code units, not by bytes: the spec is written in
        # UTF-16, and a key outside the BMP orders differently under the two.
        keys = sorted(value, key=lambda k: k.encode("utf-16-be"))
        return "{" + ",".join(f"{_escape(k)}:{canonicalize(value[k])}" for k in keys) + "}"
    raise TypeError(f"not JSON-representable: {type(value)}")


def digest(value) -> str:
    """sha256 over the canonical form, with its algorithm prefix."""
    return "sha256:" + hashlib.sha256(canonicalize(value).encode("utf-8")).hexdigest()


#: Excluded from every content digest, alongside the self-reference.
#:
#: ``DEP-VER-003``. A document's identity is what it says, not when it was
#: written. Kept in step with ``cf_types::VOLATILE_DIGEST_FIELDS``; the two
#: implementations are independent on purpose, so the constant is duplicated
#: and the *values* must match.
VOLATILE_FIELDS = ("created_at",)


def digest_of_document(document: dict, omit_field: str = "digest") -> str:
    """Digest of a document with the self-reference and volatile fields removed.

    Every ``digest`` field in the schemas is defined this way: the digest a
    document carries cannot include itself, and cannot include when it was
    written -- otherwise two identical resolutions are different lockfiles and
    ``replay --assert-identical`` can never succeed.
    """
    dropped = {omit_field, *VOLATILE_FIELDS}
    return digest({k: v for k, v in document.items() if k not in dropped})


def digest_of_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


# --------------------------------------------------------------------------- #
# CodonSequence document form
# --------------------------------------------------------------------------- #

_BASES = "ACGT"


def load_genetic_codes(root: Path) -> dict:
    """Load the vendored NCBI tables.

    Read from ``crates/cf-types/data/genetic_codes.json`` rather than
    re-deriving them, so that this module and the Rust core are demonstrably
    using the same assignments.
    """
    path = root / "crates" / "cf-types" / "data" / "genetic_codes.json"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing; the digest check needs the vendored NCBI tables. "
            "Regenerate with tools/generate_genetic_codes.py."
        )
    return json.loads(path.read_text(encoding="utf-8"))["tables"]


def _codon_index(codon: str) -> int:
    return 16 * _BASES.index(codon[0]) + 4 * _BASES.index(codon[1]) + _BASES.index(codon[2])


def codon_sequence_document(cds: str, translation_table: int, tables: dict) -> dict:
    """Build the canonical CodonSequence document for a bare coding region.

    Mirrors `CodonSequence::from_cds` in `crates/cf-types/src/seqrec.rs`,
    including its default locks (the start codon, and the terminal stop when one
    is present). No flanks, no source, no sequence id — this is the document
    form of a coding region quoted on its own, which is what a Pareto-front
    candidate or an imported seed carries.
    """
    cds = cds.upper().replace("U", "T")
    if not cds or len(cds) % 3:
        raise ValueError(f"coding sequence is not a whole number of codons: {len(cds)} nt")
    table = tables.get(str(translation_table))
    if table is None:
        raise ValueError(f"NCBI translation table {translation_table} is not defined")

    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    protein = "".join(table["aa"][_codon_index(c)] for c in codons)

    includes_stop = protein.endswith("*")
    interior = protein[:-1] if includes_stop else protein
    if "*" in interior:
        raise ValueError(f"internal stop codon in {cds}")

    locks = [{"span": {"start": 0, "end": 1}, "reason": "start_codon"}]
    if includes_stop:
        locks.append(
            {"span": {"start": len(codons) - 1, "end": len(codons)}, "reason": "stop_codon"}
        )

    document = {
        "schema_version": "1.0.0",
        "alphabet": "DNA",
        "translation_table": translation_table,
        "cds": cds,
        "codon_count": len(codons),
        "protein": protein,
        "includes_stop": includes_stop,
        "five_prime_flank": "",
        "three_prime_flank": "",
        "locked_codons": locks,
    }
    document["digest"] = digest_of_document(document)
    return document
