"""SBOL 3 export of one finished design.

The bytes come from pySBOL3 (the ``sbol3`` package). This module does not
write RDF by hand. ``sbol3.Document.validate()`` runs the library's object
rules and its SHACL shapes. A document that reports an error or a warning
is not returned.

What the file contains
----------------------
* The recommended candidate's coding sequence, with the sequence document's
  flanks when they are present, as an SBOL ``Component`` of type DNA and an
  SBOL ``Sequence`` (IUPAC DNA encoding).
* Sequence features: the CDS (SO:0000316), each flank when non-empty (no UTR
  role is claimed), each ``locked_codons`` entry (start and stop use
  SO:0000318 and SO:0000319; other lock reasons get a name and a description
  only), and each rationale annotation that has a codon span inside the CDS.
* The sequence document's ``protein`` string, when it is present, as a
  separate protein component. The writer copies that field. It does not
  translate the CDS again.
* The host profile document and the constraint set document, as annotations
  in the ``https://dna_design_suite.org/sbol`` namespace. SBOL 3 has no
  codon-usage class. An SBOL ``Constraint`` is a topological relation between
  features, so this export does not create any.

What the file does not contain
-------------------------------
* An import path. Nothing in this package reads SBOL.
* A round trip. Reading the file back does not reconstruct a design job.
* The other candidates on the front. FASTA and CSV still list those.
* A claim that a flank is a UTR, a promoter, or an RBS.

A missing ``sbol3`` install is ``CF-SBOL-001``. A document this exporter
cannot represent, or that ``Document.validate()`` rejects, is ``CF-SBOL-002``.
"""

from __future__ import annotations

import importlib.util
import json
import math
from typing import Any

from dna_design_suite.codon.orchestrator.design import DesignError

NAMESPACE = "https://dna_design_suite.org/sbol"
_PREDICATE = NAMESPACE + "#"

EXPORT_SCOPE_TEXT = (
    "SBOL 3 export only. This file is not an SBOL import and is not a round trip. "
    "It contains the recommended coding sequence, sequence features from the sequence "
    "document and from the recommended candidate's rationale annotations, the host "
    "profile document, and the constraint set document. Codon usage and design "
    "constraints are annotations in the https://dna_design_suite.org/sbol namespace. "
    "SBOL 3 has no codon-usage class. An SBOL Constraint describes the relative "
    "position of features, so this export writes no sbol:Constraint objects."
)

_SO_START_CODON = "https://identifiers.org/SO:0000318"
_SO_STOP_CODON = "https://identifiers.org/SO:0000319"
_LOCK_ROLES = {"start_codon": _SO_START_CODON, "stop_codon": _SO_STOP_CODON}
_DNA = set("ACGT")


class SbolExportError(DesignError):
    """An export that did not produce a file. The design itself may have succeeded."""

    def __init__(self, code: str, detail: str) -> None:
        super().__init__(code, detail, exit_code=2)


def library_installed() -> bool:
    """Whether ``sbol3`` can be imported, without importing it.

    ``find_spec`` does not run the package. A missing optional library must
    stay a refusal, not an import error from inside a dependency.
    """
    try:
        return importlib.util.find_spec("sbol3") is not None
    except (ImportError, ValueError):
        return False


def to_sbol(
    front: dict[str, Any],
    sequence_document: dict[str, Any],
    host_profile: dict[str, Any],
    constraint_set: dict[str, Any],
    design_host: str | None = None,
) -> str:
    """Render the four documents as sorted N-Triples.

    Sorted N-Triples is the ``sbol3`` serializer that orders its lines, so two
    calls on the same documents produce the same bytes.
    """
    if not library_installed():
        raise SbolExportError(
            "CF-SBOL-001",
            "SBOL 3 export needs the optional sbol3 library (pySBOL3). "
            "Install it with pip install 'dna-design-suite[sbol]' "
            "or pip install 'sbol3>=1.2,<2', then retry. "
            "DNA Design Suite does not write SBOL RDF itself.",
        )
    if not isinstance(front, dict) or not isinstance(sequence_document, dict):
        raise SbolExportError("CF-SBOL-002", "SBOL export needs a front and a sequence document")
    if not isinstance(host_profile, dict) or not isinstance(constraint_set, dict):
        raise SbolExportError(
            "CF-SBOL-002",
            "SBOL export needs the host profile and the constraint set as objects",
        )

    import sbol3

    candidate = _recommended(front)
    cds = _cds(candidate)
    document_cds = sequence_document.get("cds")
    if document_cds != cds:
        raise SbolExportError(
            "CF-SBOL-002",
            "the sequence document cds does not match the recommended candidate",
        )
    flank5 = _dna_field(sequence_document, "five_prime_flank")
    flank3 = _dna_field(sequence_document, "three_prime_flank")
    _require_dna(cds, "cds")
    construct = flank5 + cds + flank3
    table = sequence_document.get("translation_table")
    if not isinstance(table, int):
        raise SbolExportError(
            "CF-SBOL-002",
            "the sequence document has no integer translation_table",
        )

    token = _display_token(str(front.get("job_id") or "job"), "job")
    base = f"{NAMESPACE}/{token}"
    sequence = sbol3.Sequence(
        f"{base}/recommended_sequence",
        elements=construct,
        encoding=sbol3.IUPAC_DNA_ENCODING,
        name=str(candidate.get("candidate_id") or "recommended"),
    )
    component = sbol3.Component(
        f"{base}/recommended_cds",
        types=[sbol3.SBO_DNA],
        sequences=[sequence],
        name=str(candidate.get("candidate_id") or "recommended"),
        description=_component_description(candidate, front, design_host, table),
    )
    features, omitted = _features(sbol3, sequence, candidate, sequence_document, flank5, cds)
    component.features = features

    members: list[Any] = [component, sequence]
    protein = sequence_document.get("protein")
    if isinstance(protein, str) and protein:
        protein_sequence = sbol3.Sequence(
            f"{base}/protein_sequence",
            elements=protein,
            encoding=sbol3.IUPAC_PROTEIN_ENCODING,
            name="protein",
        )
        protein_component = sbol3.Component(
            f"{base}/protein",
            types=[sbol3.SBO_PROTEIN],
            sequences=[protein_sequence],
            name="protein",
            description=(
                "Amino-acid string copied from the sequence document's protein field. "
                f"Translation table {table}. This exporter does not translate the CDS."
            ),
        )
        members.extend([protein_component, protein_sequence])

    collection = sbol3.Collection(
        f"{base}/design",
        members=members,
        name=token,
        description=EXPORT_SCOPE_TEXT,
    )
    _annotate(sbol3, collection, _PREDICATE + "export_scope", EXPORT_SCOPE_TEXT)
    _annotate(sbol3, collection, _PREDICATE + "sbol3_version", sbol3.__version__)
    _annotate(sbol3, collection, _PREDICATE + "serialization", "sorted N-Triples")
    _annotate(sbol3, collection, _PREDICATE + "design_host", design_host)
    _annotate(sbol3, collection, _PREDICATE + "translation_table", table)
    _annotate(sbol3, collection, _PREDICATE + "candidate_id", candidate.get("candidate_id"))
    _annotate_document(sbol3, collection, "host_profile_json", host_profile)
    _annotate_document(sbol3, collection, "constraint_set_json", constraint_set)
    _annotate(sbol3, collection, _PREDICATE + "host_profile_host_id", host_profile.get("host_id"))
    _annotate(sbol3, collection, _PREDICATE + "host_display_name", host_profile.get("display_name"))
    _annotate(sbol3, collection, _PREDICATE + "host_taxid", host_profile.get("taxid"))
    _annotate(sbol3, collection, _PREDICATE + "host_strain", host_profile.get("strain"))
    _annotate(
        sbol3,
        collection,
        _PREDICATE + "constraint_set_id",
        constraint_set.get("constraint_set_id"),
    )
    summaries = _constraint_summaries(constraint_set)
    if summaries:
        sbol3.TextProperty(
            collection,
            _PREDICATE + "constraint_summary",
            0,
            math.inf,
            initial_value=summaries,
        )
    if omitted:
        sbol3.TextProperty(
            collection,
            _PREDICATE + "omitted_feature",
            0,
            math.inf,
            initial_value=omitted,
        )

    document = sbol3.Document()
    document.add(members)
    document.add(collection)
    document.bind("dds", _PREDICATE)
    report = document.validate()
    issues = list(report.errors) + list(report.warnings)
    if issues:
        detail = "; ".join(str(issue) for issue in issues)
        raise SbolExportError(
            "CF-SBOL-002",
            "sbol3.Document.validate() reported "
            f"{len(report.errors)} error(s) and {len(report.warnings)} warning(s): {detail}",
        )
    return document.write_string(sbol3.SORTED_NTRIPLES)


def _recommended(front: dict[str, Any]) -> dict[str, Any]:
    candidates = front.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise SbolExportError("CF-SBOL-002", "the front has no candidates to export")
    recommended = front.get("recommended")
    if not isinstance(recommended, dict) or not recommended.get("candidate_id"):
        raise SbolExportError(
            "CF-SBOL-002",
            "the front does not name a recommended candidate",
        )
    wanted = recommended["candidate_id"]
    matches = [item for item in candidates if item.get("candidate_id") == wanted]
    if len(matches) != 1 or not isinstance(matches[0], dict):
        raise SbolExportError(
            "CF-SBOL-002",
            f"recommended candidate {wanted} is not exactly one row of the front",
        )
    return matches[0]


def _cds(candidate: dict[str, Any]) -> str:
    sequences = candidate.get("sequences")
    if not isinstance(sequences, list) or not sequences or not isinstance(sequences[0], dict):
        raise SbolExportError("CF-SBOL-002", "the recommended candidate has no sequence")
    cds = sequences[0].get("cds")
    if not isinstance(cds, str) or not cds:
        raise SbolExportError("CF-SBOL-002", "the recommended candidate has no CDS")
    return cds


def _dna_field(document: dict[str, Any], key: str) -> str:
    value = document.get(key) or ""
    if not isinstance(value, str):
        raise SbolExportError("CF-SBOL-002", f"sequence document {key} is not a string")
    if value:
        _require_dna(value, key)
    return value


def _require_dna(sequence: str, label: str) -> None:
    bad = sorted({base for base in sequence if base not in _DNA})
    if bad:
        raise SbolExportError(
            "CF-SBOL-002",
            f"{label} contains {''.join(bad)}, which this IUPAC DNA export does not rewrite",
        )


def _component_description(
    candidate: dict[str, Any],
    front: dict[str, Any],
    design_host: str | None,
    table: int,
) -> str:
    host = design_host or "not recorded on this export"
    return (
        f"Recommended coding sequence {candidate.get('candidate_id')} "
        f"from job {front.get('job_id')}. Design host {host}. "
        f"Translation table {table}. "
        "The host profile and the constraint set are annotations on the design collection."
    )


def _features(
    sbol3: Any,
    sequence: Any,
    candidate: dict[str, Any],
    document: dict[str, Any],
    flank5: str,
    cds: str,
) -> tuple[list[Any], list[str]]:
    features: list[Any] = []
    omitted: list[str] = []
    offset = len(flank5)
    if flank5:
        features.append(
            _feature(
                sbol3,
                sequence,
                1,
                len(flank5),
                name="five_prime_flank",
                description=(
                    "Immutable 5' flank copied from the sequence document. "
                    "Not annotated as a UTR, promoter, or RBS."
                ),
            )
        )
    features.append(
        _feature(
            sbol3,
            sequence,
            offset + 1,
            offset + len(cds),
            name="cds",
            roles=[sbol3.SO_CDS],
            description=(
                "Coding sequence of the recommended candidate. "
                "Range is 1-based and inclusive on the construct, after the 5' flank."
            ),
        )
    )
    if document.get("three_prime_flank"):
        start = offset + len(cds) + 1
        features.append(
            _feature(
                sbol3,
                sequence,
                start,
                start + len(document["three_prime_flank"]) - 1,
                name="three_prime_flank",
                description=(
                    "Immutable 3' flank copied from the sequence document. "
                    "Not annotated as a UTR or terminator."
                ),
            )
        )
    locks = document.get("locked_codons") or []
    if not isinstance(locks, list):
        raise SbolExportError("CF-SBOL-002", "locked_codons is not a list")
    ordered_locks = sorted(
        (item for item in locks if isinstance(item, dict)),
        key=lambda item: (
            (item.get("span") or {}).get("start", 0),
            (item.get("span") or {}).get("end", 0),
            str(item.get("reason") or ""),
            str(item.get("note") or ""),
        ),
    )
    for item in ordered_locks:
        span = item.get("span") if isinstance(item.get("span"), dict) else {}
        located = _codon_range(span.get("start"), span.get("end"), offset, len(cds))
        reason = str(item.get("reason") or "locked")
        if located is None:
            omitted.append(f"locked_codons {reason} span {span} {item.get('note') or ''}".strip())
            continue
        note = item.get("note")
        description = f"Locked codon span from the sequence document. Reason {reason}."
        if isinstance(note, str) and note:
            description = f"{description} {note}"
        features.append(
            _feature(
                sbol3,
                sequence,
                located[0],
                located[1],
                name=f"locked_{reason}_{located[0]}_{located[1]}",
                roles=[_LOCK_ROLES[reason]] if reason in _LOCK_ROLES else None,
                description=description,
            )
        )
    rationale = candidate.get("rationale") if isinstance(candidate.get("rationale"), dict) else {}
    annotations = rationale.get("annotations") if isinstance(rationale, dict) else None
    if annotations is None:
        annotations = []
    if not isinstance(annotations, list):
        raise SbolExportError("CF-SBOL-002", "rationale annotations are not a list")
    ordered_notes = sorted(
        (item for item in annotations if isinstance(item, dict)),
        key=lambda item: (
            (item.get("span") or {}).get("start", 0),
            (item.get("span") or {}).get("end", 0),
            str(item.get("tag") or ""),
            str(item.get("text") or ""),
        ),
    )
    for item in ordered_notes:
        span = item.get("span") if isinstance(item.get("span"), dict) else {}
        located = _codon_range(span.get("start"), span.get("end"), offset, len(cds))
        tag = str(item.get("tag") or "rationale")
        if located is None:
            omitted.append(f"rationale {tag} span {span}")
            continue
        text = item.get("text") if isinstance(item.get("text"), str) else ""
        features.append(
            _feature(
                sbol3,
                sequence,
                located[0],
                located[1],
                name=f"rationale_{tag}_{located[0]}_{located[1]}",
                description=(
                    f"Rationale annotation copied from the recommended candidate. Tag {tag}. {text}"
                ),
            )
        )
    return features, omitted


def _feature(
    sbol3: Any,
    sequence: Any,
    start: int,
    end: int,
    *,
    name: str,
    description: str,
    roles: list[str] | None = None,
) -> Any:
    if start < 1 or end < start:
        raise SbolExportError("CF-SBOL-002", f"feature {name} has an empty range")
    length = len(sequence.elements or "")
    if end > length:
        raise SbolExportError(
            "CF-SBOL-002",
            f"feature {name} ends at {end}, past the construct length {length}",
        )
    return sbol3.SequenceFeature(
        [sbol3.Range(sequence, start, end)],
        roles=roles,
        orientation=sbol3.SO_FORWARD,
        name=name,
        description=description,
    )


def _codon_range(start: Any, end: Any, flank5: int, cds_length: int) -> tuple[int, int] | None:
    """Half-open codon span on the CDS to a 1-based inclusive nucleotide range."""
    if not isinstance(start, int) or not isinstance(end, int):
        return None
    if end <= start or start < 0 or end * 3 > cds_length:
        return None
    return flank5 + start * 3 + 1, flank5 + end * 3


def _annotate_document(sbol3: Any, owner: Any, name: str, document: dict[str, Any]) -> None:
    text = json.dumps(document, sort_keys=True, separators=(",", ":"))
    sbol3.TextProperty(owner, _PREDICATE + name, 1, 1, initial_value=text)


def _annotate(sbol3: Any, owner: Any, predicate: str, value: Any) -> None:
    if value is None or isinstance(value, (dict, list)):
        return
    if isinstance(value, bool):
        text = "true" if value else "false"
    elif isinstance(value, int) and not isinstance(value, bool) or isinstance(value, float):
        text = str(value)
    elif isinstance(value, str):
        if value == "":
            return
        text = value
    else:
        return
    sbol3.TextProperty(owner, predicate, 0, 1, initial_value=text)


def _constraint_summaries(constraint_set: dict[str, Any]) -> list[str]:
    items = constraint_set.get("constraints") or []
    if not isinstance(items, list):
        raise SbolExportError("CF-SBOL-002", "constraint_set.constraints is not a list")
    lines: list[str] = []
    for item in items:
        if not isinstance(item, dict):
            continue
        identifier = item.get("id")
        if not isinstance(identifier, str) or not identifier:
            continue
        rule = item.get("rule") if isinstance(item.get("rule"), dict) else {}
        kind = str(rule.get("kind") or "")
        extras = [f"{key}={rule[key]}" for key in sorted(rule) if key != "kind"]
        parts = [identifier, str(item.get("class") or ""), kind]
        if extras:
            parts.append(" ".join(extras))
        lines.append(" ".join(part for part in parts if part))
    return lines


def _display_token(value: str, fallback: str) -> str:
    cleaned = []
    for character in value:
        if character.isalnum() or character == "_":
            cleaned.append(character)
        else:
            cleaned.append("_")
    token = "".join(cleaned).strip("_") or fallback
    if token[0].isdigit():
        token = "n_" + token
    return token
