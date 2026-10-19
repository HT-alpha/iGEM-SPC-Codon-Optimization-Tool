"""Output writers: FASTA, GenBank, CSV, JSON, and the aligned `explain`.

SBOL 3 is `dna_design_suite.export.sbol.to_sbol`. It follows the same rule —
render documents, do not recompute them — and it needs the optional `sbol3`
library, so it does not live in this module's import graph.

Ticket P1-31. Acceptance: emitted GenBank re-ingests to the same digest.

## The one rule these all follow

**A writer renders the emitted documents. It never reaches into the job.**

Every writer here takes a `pareto_front` document and, where relevant, a
`job_status` document, and formats what it finds. None of them re-queries a
provider for a direction, recomputes a score, re-derives which tiers ran, or
reformats a warning into a shape of its own. That is not tidiness: a value
computed twice is a value that can disagree with itself, and this codebase has
already paid for that once — tier 1 and tier 2 each decided what to do with an
L0 term and reached different answers, until `TierPlan` made the decision once.
`objective_sense`, `resolved_objectives[].direction` and `progress.tiers` each
have exactly one producer, and the writers are consumers.

The corollary is that `--output json` emits the documents verbatim rather than
a rendering of them. A pipeline parsing the JSON is reading the same bytes the
schema validates.
"""

from __future__ import annotations

import csv
import io
import json
from typing import Any

LINE_WIDTH = 60
GENBANK_WIDTH = 60


def to_json(document: dict[str, Any]) -> str:
    """The document, verbatim and stably ordered.

    Sorted keys so that two runs of the same job produce byte-identical output
    — `dna-design-suite replay --assert-identical` diffs these, and a diff that fired
    on key order would be useless.
    """
    return json.dumps(document, indent=2, sort_keys=True)


# --------------------------------------------------------------------------- #
# FASTA
# --------------------------------------------------------------------------- #


def to_fasta(front: dict[str, Any]) -> str:
    """One record per candidate.

    The description line carries the objective values in `objective_order`,
    with the sense the front recorded. Not "the direction this writer thinks
    the term has": a header claiming `maximize` for a minimized term would
    misreport the design in the format most likely to be read by a human and
    least likely to be validated.
    """
    order = front["objective_order"]
    sense = front["objective_sense"]
    lines: list[str] = []

    for candidate in front["candidates"]:
        values = " ".join(
            f"{term}={_number(value)}({direction})"
            for term, value, direction in zip(
                order, candidate["objective_vector"], sense, strict=True
            )
        )
        feasible = "feasible" if candidate["constraint_report"]["feasible"] else "INFEASIBLE"
        lines.append(
            f">{candidate['candidate_id']} job={front['job_id']} "
            f"origin={candidate['provenance']['origin']} {feasible} {values}"
        )
        sequence = candidate["sequences"][0]["cds"]
        for start in range(0, len(sequence), LINE_WIDTH):
            lines.append(sequence[start : start + LINE_WIDTH])

    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# GenBank
# --------------------------------------------------------------------------- #


def to_genbank(
    front: dict[str, Any], sequence_document: dict[str, Any], candidate_index: int = 0
) -> str:
    """A GenBank record for one candidate.

    Takes the `sequence.schema.json` document as well as the front, for the
    same reason `explain` takes the status: a writer renders documents, and the
    flanks are not in the front. `pareto_front.candidates[].sequences[].cds` is
    a `CodonString` and `input_sequence` is a closed object carrying a digest
    and an objective vector — by design, since the front is about candidates
    rather than about construct layout.

    That matters for the acceptance criterion, which is a round trip:
    re-ingesting this must give back the same digest. The digest is of the
    sequence *document*, which includes the flanks, so a record carrying only
    the coding region would re-ingest to a different sequence that merely
    looked similar.
    """
    candidate = front["candidates"][candidate_index]
    cds = sequence_document["cds"]
    five = sequence_document.get("five_prime_flank", "")
    three = sequence_document.get("three_prime_flank", "")
    construct = five + cds + three

    locus = candidate["candidate_id"][:16]
    lines = [
        f"LOCUS       {locus:<16} {len(construct)} bp    DNA     linear   SYN",
        f"DEFINITION  DNA Design Suite design {candidate['candidate_id']}.",
        f"ACCESSION   {candidate['candidate_id']}",
        "VERSION     1",
        "KEYWORDS    dna_design_suite.",
        "SOURCE      synthetic construct",
        "  ORGANISM  synthetic construct",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(construct)}",
        '                     /mol_type="other DNA"',
        f"     CDS             {len(five) + 1}..{len(five) + len(cds)}",
        "                     /codon_start=1",
        f"                     /transl_table={sequence_document['translation_table']}",
        f'                     /note="origin={candidate["provenance"]["origin"]}"',
        f'                     /note="front={front["front_id"]}"',
    ]
    for term, value, direction in zip(
        front["objective_order"],
        candidate["objective_vector"],
        front["objective_sense"],
        strict=True,
    ):
        lines.append(f'                     /note="{term}={_number(value)} ({direction})"')

    lines.append("ORIGIN")
    for offset in range(0, len(construct), GENBANK_WIDTH):
        chunk = construct[offset : offset + GENBANK_WIDTH].lower()
        blocks = " ".join(chunk[i : i + 10] for i in range(0, len(chunk), 10))
        lines.append(f"{offset + 1:>9} {blocks}")
    lines.append("//")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# CSV
# --------------------------------------------------------------------------- #


def to_csv(front: dict[str, Any]) -> str:
    """One row per candidate, one column per objective.

    Column headers carry the sense, taken from the front. A CSV of bare numbers
    with no direction is the format most likely to be read by something that
    then optimizes the wrong way.
    """
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")

    header = ["candidate_id", "origin", "feasible", "nondomination_rank"]
    header += [
        f"{term} ({direction})"
        for term, direction in zip(front["objective_order"], front["objective_sense"], strict=True)
    ]
    header.append("cds")
    writer.writerow(header)

    for candidate in front["candidates"]:
        row = [
            candidate["candidate_id"],
            candidate["provenance"]["origin"],
            candidate["constraint_report"]["feasible"],
            candidate["nondomination_rank"],
        ]
        row += [_number(value) for value in candidate["objective_vector"]]
        row.append(candidate["sequences"][0]["cds"])
        writer.writerow(row)

    return buffer.getvalue()


# --------------------------------------------------------------------------- #
# explain
# --------------------------------------------------------------------------- #


def explain(front: dict[str, Any], status: dict[str, Any], candidate_index: int = 0) -> str:
    """The aligned per-candidate rendering.

    Reads `progress.tiers` and `warnings` as they are. The tier list is
    rendered per tier because it *is* per tier — collapsing it to a percentage
    would lose the distinction between a job that ran two tiers and one that
    skipped one, which is the whole reason the structure is not a pair of
    counters. Warnings keep their `(code, stage, subject)` triple for the same
    reason: the triple is what lets a reader attribute one.
    """
    candidate = front["candidates"][candidate_index]
    out: list[str] = []

    out.append(f"candidate  {candidate['candidate_id']}")
    out.append(f"job        {status['job_id']}  {status['state']}")
    out.append(f"lockfile   {front['landscape_lockfile_digest']}")
    out.append(f"origin     {candidate['provenance']['origin']}")
    out.append("")

    recommended = front.get("recommended")
    if recommended:
        out.append(f"recommended  {recommended['candidate_id']}  rule={recommended['rule']}")
        informed = ", ".join(recommended.get("informed_by") or []) or "(none)"
        excluded = ", ".join(recommended.get("not_informed_by") or []) or "(none)"
        out.append(f"  informed_by      {informed}")
        out.append(f"  not_informed_by  {excluded}")
        out.append("")

    out.extend(_tradeoff_lines(front))

    width = max((len(term) for term in front["objective_order"]), default=0)
    out.append("objectives")
    for term, value, direction in zip(
        front["objective_order"],
        candidate["objective_vector"],
        front["objective_sense"],
        strict=True,
    ):
        display = "UNKNOWN" if value is None else _number(value)
        out.append(f"  {term:<{width}}  {display:>14}  {direction}")
    out.append("")

    out.append("tiers")
    tier_width = max((len(t["tier"]) for t in status["progress"]["tiers"]), default=0)
    for tier in status["progress"]["tiers"]:
        detail = ""
        if tier.get("generation") is not None:
            detail = f"  {tier['generation']}/{tier.get('max_generations')}"
        out.append(f"  {tier['tier']:<{tier_width}}  {tier['state']}{detail}")
    out.append("")

    report = candidate["constraint_report"]
    out.append(f"constraints  {'satisfied' if report['feasible'] else 'VIOLATED'}")
    for violation in report.get("violations", []):
        out.append(f"  {violation['constraint_id']}  {violation['detail']}")
    out.append("")

    rationale = candidate.get("rationale") or {}
    annotations = rationale.get("annotations") or []
    if annotations:
        out.append("rationale")
        for annotation in annotations:
            span = annotation.get("span") or {}
            out.append(
                f"  {annotation.get('tag', '')}  {span.get('start', '')}..{span.get('end', '')}"
            )
            out.append(f"    {annotation.get('text', '')}")
        out.append("")

    if status.get("warnings"):
        out.append("warnings")
        code_width = max(len(w["code"]) for w in status["warnings"])
        for warning in status["warnings"]:
            out.append(
                f"  {warning['code']:<{code_width}}  {warning.get('stage', '?')}"
                f"/{warning.get('subject', '?')}"
            )
            out.append(f"  {'':<{code_width}}  {warning['title']}")
        out.append("")

    reduction = front["reduction"]
    out.append(
        f"front      {reduction['presented_count']} of {reduction['raw_nondominated_count']} "
        f"presented, reduction={reduction['method']}"
    )
    out.append(f"summary    termination={front['summary']['termination_reason']}")
    return "\n".join(out) + "\n"


def _tradeoff_lines(front: dict[str, Any]) -> list[str]:
    """Ranges on the emitted front. The numbers are the vectors, not a new score."""
    order = front.get("objective_order") or []
    sense = front.get("objective_sense") or []
    if len(order) < 2 or len(sense) != len(order):
        return []
    complete = []
    for candidate in front.get("candidates") or []:
        vector = candidate.get("objective_vector") or []
        if len(vector) == len(order) and all(value is not None for value in vector):
            complete.append(vector)
    if len(complete) < 2:
        return []
    lines = [f"trade-off  {len(complete)} complete candidates"]
    moved = False
    for index, term in enumerate(order):
        column = [vector[index] for vector in complete]
        low, high = min(column), max(column)
        if low != high:
            moved = True
        lines.append(f"  {term}  {_number(low)} .. {_number(high)}  {sense[index]}")
    if not moved:
        return []
    lines.append("")
    return lines


def _number(value: float | None) -> str:
    """A stable rendering, so two runs of the same job produce the same bytes."""
    if value is None:
        return ""
    return f"{value:.6f}"
