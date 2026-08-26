"""FASTA reading, and the conversion to the canonical sequence representation.

Ticket P1-09. Test: `TST-INGE-001`.
Contract: `planning/02-schemas/sequence.schema.json`.

Hand-rolled rather than delegating to Biopython, for one reason: Biopython is
lenient about malformed headers and stray characters, and this is the boundary
where leniency turns into a silent wrong answer. A record with an empty header,
a stray sequence block before any header, or a character outside the accepted
alphabet is rejected here with a code the caller can branch on.

WHAT IS ACCEPTED SILENTLY, AND WHAT IS NOT:

  lower case      accepted, upper-cased          (universally used for masking)
  U / u           accepted, transcoded to T      (RNA input), and REPORTED
  whitespace      accepted, stripped             (line wrapping)
  '-' or '.'      REJECTED                       (alignment gaps are not a CDS)
  N and IUPAC     REJECTED for a coding sequence (CF-SEQ-002)

The last one is the important one. An ambiguity code in a sequence about to be
codon-optimized is not something to guess about: `N` in a CDS means the input is
a consensus or a draft assembly, and translating it would invent a residue.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from codonforge.ingest.errors import IngestError, IngestReport, SourceLocation

_ALLOWED = set("ACGTUacgtu")
_GAP = set("-.")
_AMBIGUITY = set("NnRrYySsWwKkMmBbDdHhVv")


@dataclass
class FastaRecord:
    """One record, before any interpretation as a coding sequence."""

    identifier: str
    description: str
    sequence: str
    line: int

    @property
    def header(self) -> str:
        return f"{self.identifier} {self.description}".strip()


def parse_fasta(text: str, *, uri: str | None = None) -> list[FastaRecord]:
    """Parse FASTA text into records, rejecting malformed input."""
    records: list[FastaRecord] = []
    identifier: str | None = None
    description = ""
    chunks: list[str] = []
    start_line = 0

    def flush() -> None:
        if identifier is not None:
            records.append(
                FastaRecord(
                    identifier=identifier,
                    description=description,
                    sequence="".join(chunks),
                    line=start_line,
                )
            )

    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            if not header:
                raise IngestError(
                    "CF-SEQ-001",
                    "FASTA header is empty; a record needs an identifier",
                    location=SourceLocation(uri=uri, line=number),
                )
            identifier, _, description = header.partition(" ")
            description = description.strip()
            chunks = []
            start_line = number
        else:
            if identifier is None:
                raise IngestError(
                    "CF-SEQ-001",
                    "sequence data appears before any '>' header",
                    location=SourceLocation(uri=uri, line=number),
                )
            chunks.append(re.sub(r"\s+", "", line))
    flush()

    if not records:
        raise IngestError(
            "CF-SEQ-001",
            "no FASTA records found",
            location=SourceLocation(uri=uri),
        )
    duplicates = _duplicates([r.identifier for r in records])
    if duplicates:
        raise IngestError(
            "CF-SEQ-013",
            f"duplicate FASTA identifiers: {sorted(duplicates)}",
            location=SourceLocation(uri=uri),
        )
    return records


def _duplicates(values: list[str]) -> set[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return duplicates


def normalize_nucleotides(
    sequence: str,
    *,
    location: SourceLocation,
    report: IngestReport | None = None,
) -> str:
    """Upper-case, transcode RNA, and reject anything that is not a base.

    Returns DNA over `ACGT`. Every transformation applied is recorded in
    `report`, because a silently transcoded sequence is one the user cannot
    reconcile with their input file.
    """
    if not sequence:
        raise IngestError("CF-SEQ-012", "sequence is empty", location=location)

    gaps = _GAP & set(sequence)
    if gaps:
        raise IngestError(
            "CF-SEQ-002",
            f"sequence contains alignment gap character(s) {sorted(gaps)}; "
            "an aligned sequence is not a coding sequence",
            location=location,
        )
    ambiguous = _AMBIGUITY & set(sequence)
    if ambiguous:
        raise IngestError(
            "CF-SEQ-002",
            f"sequence contains ambiguity code(s) {sorted(ambiguous)}; a coding "
            "sequence to be optimized must be unambiguous, because translating "
            "an ambiguity code would invent a residue",
            location=location,
        )
    unknown = set(sequence) - _ALLOWED
    if unknown:
        raise IngestError(
            "CF-SEQ-002",
            f"sequence contains character(s) {sorted(unknown)} that are not nucleotides",
            location=location,
        )

    if report is not None:
        if any(c in "Uu" for c in sequence):
            report.note(
                f"{location.describe()}: RNA input transcoded U->T; the canonical "
                "representation is always DNA"
            )
        if any(c.islower() for c in sequence):
            report.note(f"{location.describe()}: lower-case bases upper-cased")

    return sequence.upper().replace("U", "T")


def read_cds(
    text: str,
    *,
    translation_table: int,
    uri: str | None = None,
    record_id: str | None = None,
    report: IngestReport | None = None,
) -> dict:
    """Read a single coding sequence and validate it through the core.

    Returns the `sequence.schema.json` document. Validation is the core's, not
    this module's: frame, stop placement and translation are checked by
    `CodonSequence`, and its `CF-SEQ-*` code is preserved.
    """
    from codonforge import _core

    records = parse_fasta(text, uri=uri)
    if record_id is not None:
        matching = [r for r in records if r.identifier == record_id]
        if not matching:
            raise IngestError(
                "CF-SEQ-001",
                f"no record {record_id!r} in this file; found {[r.identifier for r in records]}",
                location=SourceLocation(uri=uri),
            )
        record = matching[0]
    elif len(records) != 1:
        raise IngestError(
            "CF-SEQ-001",
            f"expected one record, found {len(records)}; name one with record_id",
            location=SourceLocation(uri=uri),
        )
    else:
        record = records[0]

    location = SourceLocation(uri=uri, record_id=record.identifier, line=record.line)
    cds = normalize_nucleotides(record.sequence, location=location, report=report)

    try:
        document = _core.sequence_document(cds, translation_table, None, None, None)
    except ValueError as exc:
        raise IngestError.from_core(exc, location=location) from exc

    import json as _json

    parsed = _json.loads(document)
    parsed["source"] = {
        "kind": "fasta",
        "record_id": record.identifier,
        **({"uri": uri} if uri else {}),
    }
    # `source` is part of the digested payload, so the digest must be recomputed
    # after adding it rather than kept from the core's first pass.
    parsed.pop("digest", None)
    parsed["digest"] = _core.digest_json(_json.dumps(parsed), "digest")
    return parsed


def read_cds_file(path: str | Path, *, translation_table: int, **kwargs) -> dict:
    path = Path(path)
    return read_cds(
        path.read_text(encoding="utf-8"),
        translation_table=translation_table,
        uri=str(path),
        **kwargs,
    )


def iter_records(text: str, *, uri: str | None = None) -> Iterator[FastaRecord]:
    yield from parse_fasta(text, uri=uri)
