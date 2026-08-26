"""GenBank reading and construct-context assembly.

Ticket P1-10. Contract: `planning/02-schemas/construct_context.schema.json`.

THIS MODULE OWNS THE ONLY NUCLEOTIDE-TO-CODON CONVERSION IN THE SYSTEM.
Everything downstream of ingest is codon-indexed; GenBank feature locations are
nucleotide-indexed. So this is the one place a frame error is possible, and it
is therefore the one place that checks for it explicitly: a CDS whose span is
not a whole number of codons is rejected here, with its feature location in the
message, rather than becoming a confusing failure inside `CodonSequence` three
tickets later.

Biopython parses the format; the interpretation is ours. GenBank in the wild
contains inconsistent CDS spans, features running off the end of the record, and
`/translation` qualifiers that disagree with the underlying sequence — all of
which Biopython will happily hand over.
"""

from __future__ import annotations

import io
import json
from typing import Any

from codonforge.ingest.errors import IngestError, IngestReport, SourceLocation

# Feature types mapped onto `construct_context.elements[].kind`.
_FEATURE_KIND = {
    "CDS": "cds",
    "promoter": "promoter",
    "RBS": "rbs",
    "terminator": "terminator",
    "5'UTR": "utr5",
    "3'UTR": "utr3",
    "operator": "operator",
    "rep_origin": "origin",
    "misc_feature": "other",
    "protein_bind": "operator",
    "regulatory": "other",
}


def _parse(text: str, *, uri: str | None):
    try:
        from Bio import SeqIO
    except ImportError as exc:  # pragma: no cover
        raise IngestError(
            "CF-SEQ-001",
            "biopython is required to read GenBank; install codonforge's "
            "dependencies or supply FASTA instead",
            location=SourceLocation(uri=uri),
        ) from exc

    records = list(SeqIO.parse(io.StringIO(text), "genbank"))
    if not records:
        raise IngestError(
            "CF-SEQ-001",
            "no GenBank records found; the file may be truncated or in another format",
            location=SourceLocation(uri=uri),
        )
    return records


def extract_cds_sequences(text: str, *, uri: str | None = None) -> list[tuple[str, str]]:
    """`(identifier, nucleotide sequence)` for every CDS feature.

    Used by the seed importers, which work purely at the nucleotide level.
    """
    out: list[tuple[str, str]] = []
    for record in _parse(text, uri=uri):
        for index, feature in enumerate(record.features):
            if feature.type != "CDS":
                continue
            name = (
                feature.qualifiers.get("label", [None])[0]
                or feature.qualifiers.get("gene", [None])[0]
                or f"{record.id}_cds{index}"
            )
            out.append((name, str(feature.extract(record.seq))))
    if not out:
        raise IngestError(
            "CF-SEQ-001",
            "GenBank record contains no CDS feature",
            location=SourceLocation(uri=uri),
        )
    return out


def read_construct(
    text: str,
    *,
    host_id: str,
    translation_table: int,
    assembly_standard: str = "none",
    construct_id: str = "cons_imported",
    uri: str | None = None,
    report: IngestReport | None = None,
) -> dict[str, Any]:
    """Build a `ConstructContext` from an annotated GenBank record.

    Cistrons are inferred from CDS features plus the nearest upstream RBS, which
    is the convention the schema documents. An inference that cannot be made
    confidently is reported rather than guessed: a CDS with no upstream RBS gets
    a cistron with `rbs_element_id: null` and a note, not a fabricated pairing.
    """
    from codonforge import _core

    report = report if report is not None else IngestReport()
    records = _parse(text, uri=uri)
    if len(records) != 1:
        raise IngestError(
            "CF-SEQ-001",
            f"expected one GenBank record, found {len(records)}",
            location=SourceLocation(uri=uri),
        )
    record = records[0]
    full_sequence = str(record.seq).upper().replace("U", "T")
    location = SourceLocation(uri=uri, record_id=record.id)

    unknown = set(full_sequence) - set("ACGT")
    if unknown:
        raise IngestError(
            "CF-SEQ-002",
            f"construct sequence contains {sorted(unknown)}; ambiguity codes and "
            "gaps are not accepted in a construct to be designed",
            location=location,
        )

    elements: list[dict[str, Any]] = []
    cistrons: list[dict[str, Any]] = []

    for index, feature in enumerate(record.features):
        if feature.type == "source":
            continue
        kind = _FEATURE_KIND.get(feature.type)
        if kind is None:
            report.note(f"{location.describe()}: feature type {feature.type!r} ignored")
            continue

        start = int(feature.location.start)
        end = int(feature.location.end)
        if start < 0 or end > len(full_sequence) or start >= end:
            raise IngestError(
                "CF-SEQ-005",
                f"{feature.type} feature {index} spans {start}..{end}, which is "
                f"outside the {len(full_sequence)} nt record",
                location=location,
            )

        # THE frame check. A CDS span that is not a whole number of codons
        # cannot be codon-indexed, and every downstream module assumes it can.
        if kind == "cds" and (end - start) % 3 != 0:
            raise IngestError(
                "CF-SEQ-001",
                f"CDS feature {index} spans {start}..{end} = {end - start} nt, "
                "which is not a whole number of codons. Nucleotide-to-codon "
                "conversion happens here and nowhere else, so this is rejected "
                "at ingest rather than downstream.",
                location=location,
            )

        name = (
            feature.qualifiers.get("label", [None])[0]
            or feature.qualifiers.get("gene", [None])[0]
            or f"{feature.type.lower()}_{index}"
        )
        element_id = f"el_{_slug(name)}_{index}"
        elements.append(
            {
                "element_id": element_id,
                "kind": kind,
                "span": {
                    "start_nt": start,
                    "end_nt": end,
                    "strand": "+" if feature.location.strand != -1 else "-",
                },
                "editable": kind == "cds",
                "name": str(name),
            }
        )

        if kind == "cds":
            cds_bytes = full_sequence[start:end]
            try:
                _core.translate(cds_bytes, translation_table)
            except ValueError as exc:
                raise IngestError.from_core(exc, location=location) from exc

            declared = feature.qualifiers.get("translation", [None])[0]
            if declared:
                actual = _core.translate(cds_bytes, translation_table).rstrip("*")
                if actual != declared:
                    report.warn(
                        f"{location.describe()}: CDS {name!r} /translation "
                        "disagrees with the underlying sequence; the sequence wins"
                    )

            rbs = _nearest_upstream(elements, start, "rbs")
            if rbs is None:
                report.note(
                    f"{location.describe()}: CDS {name!r} has no upstream RBS "
                    "annotation; its cistron records rbs_element_id = null "
                    "rather than guessing"
                )
            cistrons.append(
                {
                    "cistron_id": f"cis_{_slug(name)}",
                    "cds_element_id": element_id,
                    "rbs_element_id": rbs,
                    "is_design_target": True,
                    "expected_relative_expression": 1.0,
                    "expression_prior_source": "equal_default",
                    "priority_weight": 1.0,
                    "operon_group": None,
                }
            )

    if not cistrons:
        raise IngestError(
            "CF-SEQ-001",
            "GenBank record contains no CDS feature, so there is nothing to design",
            location=location,
        )

    construct: dict[str, Any] = {
        "schema_version": "1.0.0",
        "construct_id": construct_id,
        "name": record.description or record.id,
        "host_id": host_id,
        "topology": "circular" if _is_circular(record) else "linear",
        "full_sequence": full_sequence,
        "elements": elements,
        "cistrons": cistrons,
        "assembly": {
            "standard": assembly_standard,
            "forbidden_enzymes": [],
            "vendor_profile": None,
        },
    }
    construct["digest"] = _core.digest_json(json.dumps(construct), "digest")
    return construct


def _is_circular(record) -> bool:
    topology = record.annotations.get("topology")
    return isinstance(topology, str) and topology.lower() == "circular"


def _nearest_upstream(elements: list[dict[str, Any]], start: int, kind: str) -> str | None:
    """The closest element of `kind` ending at or before `start`."""
    best: tuple[int, str] | None = None
    for element in elements:
        if element["kind"] != kind:
            continue
        end = element["span"]["end_nt"]
        if end <= start and (best is None or end > best[0]):
            best = (end, element["element_id"])
    return best[1] if best else None


def _slug(text: str) -> str:
    return "".join(c if c.isalnum() else "_" for c in str(text)).strip("_").lower()[:32] or "x"
