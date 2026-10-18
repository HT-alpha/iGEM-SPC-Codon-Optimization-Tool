"""Codon-usage tables from text a scientist already has.

Kazusa lines, a delimited codon table, or a set of coding sequences. Relative
adaptiveness is the codon count divided by the maximum count in its synonymous
family. An unobserved family stays at 0. No pseudocount is added.

The stored subset is `user_supplied` unless the caller explicitly labels the
counts `highly_expressed`. `genome_wide` is refused: CAI will not score that
subset (`CF-HOST-003`), and relabelling it silently would claim a different
quantity. `DEP-API-004`.
"""

from __future__ import annotations

import json
import re
from typing import Any

from dna_design_suite.ingest.errors import IngestError, SourceLocation
from dna_design_suite.ingest.fasta import parse_fasta

_CODON = re.compile(r"^[ACGTU]{3}$")
_KAZUSA = re.compile(r"([ACGTU]{3})\s+[0-9]+(?:\.[0-9]+)?\(\s*([0-9]+)\s*\)", re.IGNORECASE)
_SPLIT = re.compile(r"[\t,;]+")


def parse_usage_table(
    text: str,
    *,
    translation_table: int = 11,
    source_format: str = "auto",
    organism: str = "",
    label: str = "user_supplied",
) -> dict[str, Any]:
    """Return a host profile whose only block is `codon_usage`."""
    if label not in {"user_supplied", "highly_expressed"}:
        raise IngestError(
            "CF-HOST-003",
            "label must be user_supplied or highly_expressed. genome_wide is not "
            "a CAI reference: usage.cai refuses that subset. Send the same counts "
            "with label user_supplied; the source text records what the file was.",
            location=SourceLocation(),
        )
    cleaned = text.replace("\r\n", "\n").replace("\r", "\n").strip()
    if not cleaned:
        raise IngestError("CF-SEQ-012", "usage text is empty", location=SourceLocation())
    if source_format == "json" or (source_format == "auto" and cleaned[0] in "{["):
        kind = "json"
    else:
        kind = source_format if source_format != "auto" else _detect(cleaned)
    if kind == "kazusa":
        counts, citation = _kazusa_counts(cleaned)
    elif kind == "table":
        counts, citation = _delimited_counts(cleaned)
    elif kind == "cds":
        counts, citation = _cds_counts(cleaned)
    elif kind == "json":
        counts, citation = _json_counts(cleaned)
    else:
        raise IngestError(
            "CF-CLI-001",
            f"format {source_format!r} is not kazusa, table, cds, json, or auto",
            location=SourceLocation(),
        )
    if not counts:
        raise IngestError(
            "CF-SEQ-012",
            "no codon counts were found. Kazusa text needs a triplet and a "
            "parenthesized count. A table needs a codon column and a count, "
            "number, or fraction column. A CDS set needs A/C/G/T sequence.",
            location=SourceLocation(),
        )
    table, unobserved = _weights(counts, translation_table)
    who = f" {organism.strip()}." if organism.strip() else ""
    label_sentence = (
        "The caller labelled these counts as a highly expressed set."
        if label == "highly_expressed"
        else "subset is user_supplied. These weights are not a Sharp and Li "
        "highly expressed reference unless that is what the counts were."
    )
    source = (
        f"{citation}{who} Relative adaptiveness is the count divided by the "
        f"maximum count in the synonymous family. Unobserved families stay at 0. "
        f"No pseudocount was added. {label_sentence}"
    )
    return {
        "codon_usage": {
            "subset": label,
            "source": source,
            "table": table,
        },
        "unobserved_families": unobserved,
        "cai_usable": not unobserved,
    }


def _detect(text: str) -> str:
    if _KAZUSA.search(text):
        return "kazusa"
    if text.lstrip().startswith(">"):
        return "cds"
    lines = [
        line.strip() for line in text.splitlines() if line.strip() and not line.startswith("#")
    ]
    if not lines:
        return "cds"
    cells = _SPLIT.split(lines[0])
    if len(cells) >= 2 and not _CODON.match(cells[0].upper()):
        return "table"
    if len(cells) >= 2 and _CODON.match(cells[0].upper()) and _number(cells[1]) is not None:
        return "table"
    return "cds"


def _kazusa_counts(text: str) -> tuple[dict[str, float], str]:
    counts: dict[str, float] = {}
    for match in _KAZUSA.finditer(text):
        codon = _dna(match.group(1))
        counts[codon] = counts.get(codon, 0.0) + float(match.group(2))
    citation = (
        "Codon counts parsed from Kazusa text (triplet and parenthesized count; "
        "Nakamura, Gojobori and Ikemura, Nucleic Acids Research 2000)."
    )
    return counts, citation


def _json_counts(text: str) -> tuple[dict[str, float], str]:
    try:
        document = json.loads(text)
    except json.JSONDecodeError as exc:
        raise IngestError(
            "CF-HOST-002",
            f"usage JSON did not parse: {exc.msg}",
            location=SourceLocation(),
        ) from exc
    block = document.get("codon_usage") if isinstance(document, dict) else None
    table = block.get("table") if isinstance(block, dict) else None
    if not isinstance(table, dict):
        raise IngestError(
            "CF-HOST-002",
            "usage JSON needs codon_usage.table",
            location=SourceLocation(),
        )
    counts: dict[str, float] = {}
    for codon, entry in table.items():
        if not _CODON.match(str(codon).upper()):
            continue
        if isinstance(entry, dict):
            raw = entry.get("relative_adaptiveness", entry.get("count", entry.get("fraction")))
        else:
            raw = entry
        value = _number(str(raw)) if not isinstance(raw, (int, float)) else float(raw)
        if value is None:
            continue
        counts[_dna(str(codon))] = value
    citation = (
        "Relative adaptiveness read from a JSON codon_usage.table and "
        "renormalized within each synonymous family. subset stays the label "
        "on this call. This is not a downloaded species catalogue."
    )
    return counts, citation


def _delimited_counts(text: str) -> tuple[dict[str, float], str]:
    lines = [
        line.strip() for line in text.splitlines() if line.strip() and not line.startswith("#")
    ]
    header = _SPLIT.split(lines[0])
    names = [cell.strip().lower() for cell in header]
    codon_at = names.index("codon") if "codon" in names else 0
    value_at = _value_column(header)
    start = 1 if value_at is not None or "codon" in names else 0
    if value_at is None:
        value_at = 1 if codon_at == 0 else codon_at + 1
    counts: dict[str, float] = {}
    for line in lines[start:]:
        cells = _SPLIT.split(line)
        if len(cells) <= max(codon_at, value_at) or not _CODON.match(cells[codon_at].upper()):
            continue
        value = _number(cells[value_at])
        if value is None:
            continue
        codon = _dna(cells[codon_at])
        counts[codon] = counts.get(codon, 0.0) + value
    citation = (
        "Codon counts parsed from a delimited table (CoCoPUTs-style codon and "
        "count, number, or fraction column; Alexaki et al., Journal of Molecular "
        "Biology 2019, for the CoCoPUTs compilation)."
    )
    return counts, citation


def _cds_counts(text: str) -> tuple[dict[str, float], str]:
    sequences = _sequences(text)
    counts: dict[str, float] = {}
    for index, cds in enumerate(sequences):
        if len(cds) % 3:
            raise IngestError(
                "CF-SEQ-001",
                f"CDS {index} is {len(cds)} nt, not a whole number of codons",
                location=SourceLocation(record_id=f"cds[{index}]"),
            )
        for offset in range(0, len(cds), 3):
            codon = cds[offset : offset + 3]
            if any(base not in "ACGT" for base in codon):
                raise IngestError(
                    "CF-SEQ-002",
                    f"CDS {index} has a non-base in codon {codon}",
                    location=SourceLocation(record_id=f"cds[{index}]"),
                )
            counts[codon] = counts.get(codon, 0.0) + 1.0
    citation = (
        "Codon counts from the supplied coding sequences. w_i follows the "
        "definition in Sharp and Li, Nucleic Acids Research 1987: each codon "
        "divided by the maximum in its synonymous family."
    )
    return counts, citation


def _sequences(text: str) -> list[str]:
    stripped = text.strip()
    if stripped.startswith(">"):
        return [record.sequence.upper().replace("U", "T") for record in parse_fasta(stripped)]
    sequences: list[str] = []
    current: list[str] = []
    for line in stripped.splitlines():
        if line.startswith(">") or line.startswith(";"):
            if current:
                sequences.append("".join(current))
                current = []
            continue
        current.append(re.sub(r"[^ACGTUacgtu]", "", line))
    if current:
        sequences.append("".join(current).upper().replace("U", "T"))
    return [item for item in sequences if item]


def _weights(counts: dict[str, float], translation_table: int) -> tuple[dict[str, Any], list[str]]:
    from dna_design_suite import _core

    families: dict[str, list[str]] = {}
    amino: dict[str, str] = {}
    for codon in (a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"):
        members = [str(item) for item in _core.synonyms(codon, translation_table)]
        key = "".join(sorted(members))
        families.setdefault(key, members)
        if codon in {"TAA", "TAG", "TGA"}:
            amino[codon] = "*"
        else:
            amino[codon] = str(_core.translate(codon + "AAA", translation_table))[0]
    table: dict[str, Any] = {}
    unobserved: list[str] = []
    for members in families.values():
        observed = max(counts.get(member, 0.0) for member in members)
        if observed <= 0.0:
            name = amino.get(members[0], "?")
            if name not in unobserved:
                unobserved.append(name)
        for member in members:
            count = counts.get(member, 0.0)
            if observed > 0.0 and count == observed:
                weight = 1.0
            elif observed <= 0.0:
                weight = 0.0
            else:
                weight = count / observed
            table[member] = {"relative_adaptiveness": round(weight, 6)}
    unobserved.sort()
    return dict(sorted(table.items())), unobserved


def _value_column(header: list[str]) -> int | None:
    names = [cell.strip().lower() for cell in header]
    for wanted in (
        "number",
        "count",
        "counts",
        "relative_adaptiveness",
        "relative_frequency",
        "frequency_per_thousand",
        "per_thousand",
        "frequency",
        "fraction",
    ):
        if wanted in names:
            return names.index(wanted)
    return None


def _number(text: str) -> float | None:
    try:
        return float(text.strip())
    except ValueError:
        return None


def _dna(codon: str) -> str:
    return codon.upper().replace("U", "T")


def profile_json(profile: dict[str, Any]) -> str:
    return json.dumps(profile, sort_keys=True)
