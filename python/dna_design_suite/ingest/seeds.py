"""Seed importers for external tool output.

Ticket P1-12. Test: `TST-INGE-003`.
Contract: `planning/02-schemas/seed_set.schema.json`, ADR-0007.

Every adapter works purely at the nucleotide level: parse, translate, compare
with the target protein. No adapter imports a tokenizer, a model or a scoring
function from the originating tool, which is what keeps `cf-ingest` free of a
dependency on `cf-biophys-lm` and the module graph acyclic — and what makes
adding an adapter about thirty lines with no new dependency.

A rejected seed carries a POSITIONAL PROTEIN DIFF. A user who pastes the wrong
file must be told which residue disagrees, not that "the import failed"; it is
the most likely user error in the whole feature and the cheapest one to
diagnose properly.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

from dna_design_suite.ingest.errors import IngestError, IngestReport, SourceLocation
from dna_design_suite.ingest.fasta import normalize_nucleotides, parse_fasta

# Tools whose output is accepted. Matches `seed_set.schema.json`.
KNOWN_TOOLS = {
    "lineardesign",
    "cdsfold",
    "codontransformer",
    "codonbert",
    "dnachisel",
    "idt_optimizer",
    "genscript",
    "jcat",
    "manual",
    "previous_job",
    "native",
}


@dataclass
class ProteinDifference:
    position: int
    expected: str
    found: str

    def as_dict(self) -> dict[str, Any]:
        return {"position": self.position, "expected": self.expected, "found": self.found}


@dataclass
class SeedAcceptance:
    accepted: bool
    reason: str
    protein_diff: list[ProteinDifference] = field(default_factory=list)
    repaired: bool = False
    n_repair_edits: int = 0

    def as_dict(self) -> dict[str, Any]:
        return {
            "accepted": self.accepted,
            "reason": self.reason,
            "protein_diff": [d.as_dict() for d in self.protein_diff],
            "repaired": self.repaired,
            "n_repair_edits": self.n_repair_edits,
        }


def _protein_diff(expected: str, found: str) -> list[ProteinDifference]:
    """Positional differences, including a length mismatch as trailing gaps."""
    diffs: list[ProteinDifference] = []
    for index in range(max(len(expected), len(found))):
        left = expected[index] if index < len(expected) else "-"
        right = found[index] if index < len(found) else "-"
        if left != right:
            diffs.append(ProteinDifference(index, left, right))
    return diffs


def extract_sequences(tool: str, payload: str, *, uri: str | None = None) -> list[tuple[str, str]]:
    """Pull `(identifier, nucleotide sequence)` pairs out of a tool's output.

    LinearDesign and CDSfold emit FASTA or a bare sequence; CodonTransformer
    emits a bare sequence; DNA Chisel emits GenBank. Everything reduces to
    nucleotides, which is the whole point of adapting at this level.
    """
    if tool not in KNOWN_TOOLS:
        raise IngestError(
            "CF-SEQ-013",
            f"unknown seed source {tool!r}; known tools are {sorted(KNOWN_TOOLS)}",
            location=SourceLocation(uri=uri),
        )

    stripped = payload.strip()
    if not stripped:
        raise IngestError("CF-SEQ-012", f"{tool} output is empty", location=SourceLocation(uri=uri))

    if stripped.startswith(">"):
        return [(r.identifier, r.sequence) for r in parse_fasta(stripped, uri=uri)]

    if tool == "dnachisel" or stripped.startswith("LOCUS"):
        from dna_design_suite.ingest.genbank import extract_cds_sequences

        return extract_cds_sequences(stripped, uri=uri)

    # A bare sequence, possibly wrapped across lines. Reject anything that looks
    # like prose rather than silently treating it as a sequence.
    joined = "".join(stripped.split())
    if not joined:
        raise IngestError(
            "CF-SEQ-012", f"{tool} output has no sequence", location=SourceLocation(uri=uri)
        )
    return [(f"{tool}_seed", joined)]


def import_seeds(
    target_protein: str,
    translation_table: int,
    seeds: list[dict[str, Any]],
    *,
    host_id: str | None = None,
    seed_set_id: str = "seeds_imported",
    report: IngestReport | None = None,
) -> dict[str, Any]:
    """Validate seeds against the target protein, returning a `SeedSet`.

    A seed is admitted only if it translates to exactly the target protein under
    the target table. Rejections are reported per seed, never as a whole-batch
    failure: one bad file among ten should not lose the other nine.
    """
    from dna_design_suite import _core

    report = report if report is not None else IngestReport()
    entries: list[dict[str, Any]] = []
    seen_digests: set[str] = set()

    for index, spec in enumerate(seeds):
        tool = spec.get("tool", "manual")
        payload = spec.get("payload", "")
        uri = spec.get("uri")
        location = SourceLocation(uri=uri, record_id=f"seed[{index}]")

        try:
            extracted = extract_sequences(tool, payload, uri=uri)
        except IngestError as exc:
            entries.append(
                _entry(
                    f"sd_{tool}_{index}",
                    "",
                    None,
                    spec,
                    SeedAcceptance(False, "invalid_alphabet"),
                    detail=exc.detail,
                )
            )
            report.warn(f"seed {index} ({tool}): {exc.detail}")
            continue

        for offset, (_identifier, raw) in enumerate(extracted):
            seed_id = f"sd_{tool}_{index}" + (f"_{offset}" if len(extracted) > 1 else "")
            try:
                cds = normalize_nucleotides(raw, location=location, report=report)
            except IngestError as exc:
                entries.append(
                    _entry(
                        seed_id,
                        "",
                        None,
                        spec,
                        SeedAcceptance(False, "invalid_alphabet"),
                        exc.detail,
                    )
                )
                report.warn(f"{seed_id}: {exc.detail}")
                continue

            if len(cds) % 3 != 0:
                entries.append(
                    _entry(
                        seed_id,
                        cds,
                        None,
                        spec,
                        SeedAcceptance(False, "not_whole_codons"),
                        f"{len(cds)} nt is not a multiple of three",
                    )
                )
                report.warn(f"{seed_id}: {len(cds)} nt is not a whole number of codons")
                continue

            try:
                document = json.loads(
                    _core.sequence_document(cds, translation_table, None, None, None)
                )
            except ValueError as exc:
                translated = IngestError.from_core(exc, location=location)
                reason = "internal_stop" if translated.code == "CF-SEQ-003" else "invalid_alphabet"
                entries.append(
                    _entry(
                        seed_id, cds, None, spec, SeedAcceptance(False, reason), translated.detail
                    )
                )
                report.warn(f"{seed_id}: {translated.detail}")
                continue

            protein = document["protein"]
            if protein != target_protein:
                diffs = _protein_diff(target_protein, protein)
                reason = (
                    "length_mismatch" if len(protein) != len(target_protein) else "protein_mismatch"
                )
                entries.append(
                    _entry(
                        seed_id,
                        cds,
                        document["digest"],
                        spec,
                        SeedAcceptance(False, reason, protein_diff=diffs[:32]),
                        f"translates to {protein!r}, expected {target_protein!r}",
                    )
                )
                report.warn(
                    f"{seed_id}: protein mismatch at position(s) {[d.position for d in diffs[:8]]}"
                )
                continue

            if document["digest"] in seen_digests:
                entries.append(
                    _entry(
                        seed_id, cds, document["digest"], spec, SeedAcceptance(False, "duplicate")
                    )
                )
                continue
            seen_digests.add(document["digest"])
            entries.append(
                _entry(seed_id, cds, document["digest"], spec, SeedAcceptance(True, "ok"))
            )

    target_digest = _core.digest_json(json.dumps(target_protein))
    seed_set: dict[str, Any] = {
        "schema_version": "1.0.0",
        "seed_set_id": seed_set_id,
        "target_protein": target_protein,
        "target_protein_digest": target_digest,
        "translation_table": translation_table,
        "seeds": entries,
        "injection_policy": {
            "max_fraction_of_population": 0.2,
            "deduplicate": True,
            "repair_before_inject": True,
            "also_inject_mutated_copies": True,
            "protect_for_generations": 0,
        },
    }
    if host_id:
        seed_set["host_id"] = host_id
    seed_set["digest"] = _core.digest_json(json.dumps(seed_set), "digest")
    return seed_set


def _entry(
    seed_id: str,
    cds: str,
    digest: str | None,
    spec: dict[str, Any],
    acceptance: SeedAcceptance,
    detail: str | None = None,
) -> dict[str, Any]:
    entry: dict[str, Any] = {
        "seed_id": seed_id,
        "cds": cds,
        "digest": digest or ("sha256:" + "0" * 64),
        "source": {
            "tool": spec.get("tool", "manual"),
            "version": spec.get("version"),
            "parameters": spec.get("parameters", {}),
            "uri": spec.get("uri"),
        },
        "reported_metrics": spec.get("reported_metrics", {}),
        "acceptance": acceptance.as_dict(),
    }
    if detail:
        entry["acceptance"]["detail"] = detail
    return entry
