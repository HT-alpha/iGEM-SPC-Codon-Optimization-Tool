"""A real-protein panel for the harness. `P2-19`, `DEP-BEN-002`.

The sequences are the CDS of `NC_000913.3` whose headers carry that accession
and a `protein_id`. A record without them is refused: a panel that accepts a
bare string would be how an invented sequence gets into a benchmark.

`GATE-BEN-001` is not evaluated. The gate needs at least 50 proteins and a
competitor arm. This module writes that non-claim into `gates.json` so an
empty file cannot be read as a pass.

Native and CAI-greedy run in process. The weights are the relative
adaptiveness of the ribosomal-protein CDS in the same file, counted by
`codon_usage_from_cds`. The geometric mean of those weights is not computed
here and is not reported as `usage.cai`.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

from codonforge.bench.baselines import cai_greedy, native
from codonforge.bench.reference import BenchError
from codonforge.ingest.host import codon_usage_from_cds

_ROOT = Path(__file__).resolve().parents[3] / "benchmarks" / "corpus"
PANEL_PATH = _ROOT / "panel.json"
CDS_PATH = _ROOT / "NC_000913.3_cds.fna.gz"

_ACCESSION = "NC_000913.3"


def load_panel(path: Path | None = None) -> dict[str, Any]:
    document = json.loads((path or PANEL_PATH).read_text(encoding="utf-8"))
    if document.get("accession") != _ACCESSION:
        raise BenchError(f"panel accession is not {_ACCESSION}")
    if not document.get("genes"):
        raise BenchError("panel names no genes")
    evidence = document.get("expression_evidence") or {}
    if not evidence.get("citation") or not evidence.get("geo"):
        raise BenchError(
            "panel has no expression citation; a protein without one is not on this list"
        )
    if evidence.get("value") is not None:
        raise BenchError("panel carries an expression number; this build does not transcribe those")
    return document


def verify_shipped_cds(panel: dict[str, Any], path: Path | None = None) -> None:
    """The gzip on disk is the file the panel names.

    A panel that records a digest and never checks it accepts a replaced
    file. The accession guard does not catch that: another file can still
    carry `NC_000913.3` in a header someone wrote.
    """
    expected = panel.get("cds_gzip_sha256")
    if not isinstance(expected, str) or len(expected) != 64:
        raise BenchError("panel has no cds_gzip_sha256; the file's identity is unchecked")
    digest = hashlib.sha256((path or CDS_PATH).read_bytes()).hexdigest()
    if digest != expected:
        raise BenchError(
            "corpus gzip sha256 does not match the panel; "
            "this file is not the cited NC_000913.3 CDS dump"
        )


def iter_cds(path: Path | None = None):
    """Yield `(header, sequence)` from the gzipped FASTA."""
    header: str | None = None
    chunks: list[str] = []
    with gzip.open(path or CDS_PATH, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def _gene(header: str) -> str | None:
    marker = "[gene="
    if marker not in header:
        return None
    return header.split(marker, 1)[1].split("]", 1)[0]


def _location_length(header: str) -> int | None:
    marker = "[location="
    if marker not in header:
        return None
    raw = header.split(marker, 1)[1].split("]", 1)[0]
    raw = raw.removeprefix("complement(").removesuffix(")")
    if ".." not in raw:
        return None
    start, end = raw.split("..", 1)
    # Join() spans and order() are not in this file's headers. A span this
    # parser does not understand is refused rather than guessed.
    if not start.isdigit() or not end.isdigit():
        return None
    return int(end) - int(start) + 1


def _header_problem(header: str, sequence: str) -> str | None:
    """Why this record cannot enter the panel. None means it can.

    The accession check is the fabrication guard: a header that is not
    `NC_000913.3` with a `protein_id` does not enter.
    """
    if _ACCESSION not in header or "[protein_id=" not in header or "[gene=" not in header:
        return (
            "corpus CDS has no NC_000913.3 protein_id; a sequence without an "
            "accession is not this panel"
        )
    span = _location_length(header)
    if span is None or span != len(sequence):
        return "corpus CDS length does not match its location span"
    if len(sequence) % 3 != 0:
        return "corpus CDS is not a whole number of codons"
    return None


def select_genes(
    panel: dict[str, Any] | None = None,
    cds_path: Path | None = None,
) -> list[dict[str, Any]]:
    """The named genes, each checked against its header."""
    panel = panel if panel is not None else load_panel()
    if cds_path is None:
        verify_shipped_cds(panel)
    wanted = list(panel["genes"])
    found: dict[str, dict[str, Any]] = {}
    for header, sequence in iter_cds(cds_path):
        gene = _gene(header)
        if gene not in wanted or gene in found:
            continue
        problem = _header_problem(header, sequence)
        if problem:
            raise BenchError(problem)
        if not sequence.startswith("ATG"):
            raise BenchError(f"{gene} does not start with ATG; it is not on this panel")
        protein_id = header.split("[protein_id=", 1)[1].split("]", 1)[0]
        product = header.split("[protein=", 1)[1].split("]", 1)[0]
        found[gene] = {
            "gene": gene,
            "product": product,
            "protein_id": protein_id,
            "cds": sequence,
            "codons": len(sequence) // 3,
        }
    missing = [gene for gene in wanted if gene not in found]
    if missing:
        raise BenchError(f"panel genes missing from the CDS file: {missing}")
    return [found[gene] for gene in wanted]


def ribosomal_usage(cds_path: Path | None = None) -> dict[str, float]:
    """Relative adaptiveness from ribosomal-protein CDS. Counts, not a guessed table."""
    sequences = []
    for header, sequence in iter_cds(cds_path):
        gene = _gene(header) or ""
        if not gene.startswith(("rpl", "rps", "rpm")):
            continue
        if _header_problem(header, sequence):
            continue
        sequences.append(sequence)
    if len(sequences) < 30:
        raise BenchError(
            f"only {len(sequences)} ribosomal CDS; "
            "the reference set is not the genome's ribosomal proteins"
        )
    table = codon_usage_from_cds(sequences, translation_table=11, subset="ribosomal_protein_cds")
    return {codon: float(row["relative_adaptiveness"]) for codon, row in table["table"].items()}


def run_baselines(directory: Path, cds_path: Path | None = None) -> dict[str, Any]:
    """Write native and CAI-greedy sequences. No competitor, no gate pass."""
    from codonforge import _core

    panel = load_panel()
    genes = select_genes(panel, cds_path)
    usage = ribosomal_usage(cds_path)
    rows = []
    fasta_records: list[str] = []
    for gene in genes:
        cds = gene["cds"]
        kept = native(cds)
        greedy = cai_greedy(cds, panel["translation_table"], usage)
        if _core.translate(kept, panel["translation_table"]) != _core.translate(
            greedy, panel["translation_table"]
        ):
            raise BenchError(f"{gene['gene']}: CAI-greedy changed the translation")
        changed = sum(
            1
            for index in range(0, len(cds), 3)
            if kept[index : index + 3] != greedy[index : index + 3]
        )
        rows.append(
            {
                "gene": gene["gene"],
                "protein_id": gene["protein_id"],
                "product": gene["product"],
                "codons": gene["codons"],
                "codons_changed": changed,
                "translation_match": True,
            }
        )
        fasta_records.append(
            f">{gene['gene']}|native\n{kept}\n>{gene['gene']}|cai_greedy\n{greedy}\n"
        )
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "baselines.fasta").write_text("".join(fasta_records), encoding="utf-8")
    manifest = {
        "panel": panel["id"],
        "accession": panel["accession"],
        "n_proteins": len(rows),
        "arms": ["native", "cai_greedy"],
        "seeds": 0,
        "seeds_why": "both arms are deterministic",
        "external_tools": [],
    }
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (directory / "metrics.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    gates = [
        {
            "gate": "GATE-BEN-001",
            "claimed": False,
            "n_proteins": len(rows),
            "reason": panel["claim"],
        }
    ]
    (directory / "gates.json").write_text(json.dumps(gates, indent=2), encoding="utf-8")
    return manifest
