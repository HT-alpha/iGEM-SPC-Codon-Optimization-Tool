"""A real-protein panel for the harness. `P2-19`, `DEP-BEN-002`.

The sequences are the CDS of `NC_000913.3` whose headers carry that accession
and a `protein_id`. A record without them is refused: a panel that accepts a
bare string would be how an invented sequence gets into a benchmark.

`GATE-BEN-001` is claimed only as the CAI-only comparison in
`benchmarks/corpus/gates.json`: 50 of 50 proteins on this record.
The competitor the original gate asks for
is a fold on `structure.mfe_global`.
LinearFold returns that global minimum free energy and a dot-bracket.
`structure.local_accessibility` is windowed ViennaRNA accessibility.
Those are not the same number, and `DEP-STRU-002` keeps global MFE off
the engine seam, so the fronts are not compared.

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

from dna_design_suite.codon.bench.baselines import cai_greedy, native
from dna_design_suite.codon.bench.reference import BenchError
from dna_design_suite.ingest.host import codon_usage_from_cds

_ROOT = Path(__file__).resolve().parents[4] / "benchmarks" / "corpus"
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
    expected_selected = panel.get("selected_cds_sha256")
    if not isinstance(expected_selected, str) or len(expected_selected) != 64:
        raise BenchError("panel has no selected_cds_sha256; the expanded gene list is unchecked")
    selected = _selected_cds_sha256(panel, path)
    if selected != expected_selected:
        raise BenchError(
            "selected CDS sha256 does not match the panel; "
            "the expanded gene list is not the cited sequences"
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


def _selected_cds_sha256(panel: dict[str, Any], path: Path | None = None) -> str:
    """Digest of the named genes in panel order, sequences included.

    The gzip digest does not change when the gene list changes. This one does.
    """
    wanted = list(panel["genes"])
    sequences: dict[str, str] = {}
    for header, sequence in iter_cds(path):
        gene = _gene(header)
        if gene not in wanted or gene in sequences:
            continue
        problem = _header_problem(header, sequence)
        if problem:
            raise BenchError(problem)
        sequences[gene] = sequence
    missing = [gene for gene in wanted if gene not in sequences]
    if missing:
        raise BenchError(f"panel genes missing from the CDS file: {missing}")
    payload = "".join(f"{gene}\n{sequences[gene]}\n" for gene in wanted)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


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
    from dna_design_suite import _core

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
            "comparable": False,
            "selected_cds_sha256": panel["selected_cds_sha256"],
            "reason": panel["claim"],
        }
    ]
    (directory / "gates.json").write_text(json.dumps(gates, indent=2), encoding="utf-8")
    return manifest


def cai_gate_document(rows: list[dict[str, Any]], *, selected_cds_sha256: str) -> dict[str, Any]:
    """The `gates.json` object `run_cai_only` writes for these rows.

    `check_cai_gate` accepts it only when every recorded score meets
    `cai_greedy`. A protein below that baseline is listed in
    `degenerate_cases.excluded` and the claim is not made.
    """
    beats_native = all(row["codonforge_cai"] + 1e-12 >= row["native_cai"] for row in rows)
    behind = [row["gene"] for row in rows if row["codonforge_cai"] + 1e-12 < row["cai_greedy_cai"]]
    claimed_n = len(rows) - len(behind)
    definition = (
        "CAI-only. One axis, usage.cai, maximized. Same ribosomal-protein "
        "relative adaptiveness for every arm, including codons whose count in "
        "the 54 ribosomal CDS is 0. Arms: the native CDS, cai_greedy, and one "
        "DNA Design Suite design (budget fast, 4 generations, seed 7, no structure "
        "objective, no forbidden motif). Claimed only when the DNA Design Suite "
        "value is at least cai_greedy on every protein. The MFE comparison "
        "against LinearDesign or CDSfold is not this record."
    )
    if behind:
        exclusion = f" ({', '.join(behind)} excluded, see degenerate_cases)"
        explanation = (
            "Excluded because the recorded DNA Design Suite score is below cai_greedy. "
            "usage.cai is 0 while any counted codon has weight 0."
        )
        details = []
        for row in rows:
            if row["gene"] not in behind:
                continue
            sites = ", ".join(
                f"{site['codon']} at codon {site['index']}" for site in row["zero_weight_codons"]
            )
            details.append(
                f"{row['gene']} ({row['protein_id']}) DNA Design Suite {row['codonforge_cai']} "
                f"vs cai_greedy {row['cai_greedy_cai']}; {sites or 'no zero-weight sense codon'}"
            )
        reason = (
            definition + " Not claimed. " + " ".join(details) + " AGG has count 0 in the 54 "
            "ribosomal CDS, so its weight is 0. usage.cai is 0 while any counted codon has "
            "weight 0, and this 4-generation search did not remove those codons. cai_greedy "
            "replaces each codon with the highest-weight synonym."
        )
    else:
        exclusion = ""
        explanation = (
            "No protein is excluded. aceB and thrA contain AGG, which has weight 0 in the "
            "54 ribosomal CDS, so usage.cai is 0 while any counted AGG remains. DEP-SRCH-020 "
            "linearizes search on the sum of log-weights in that case. On this record both "
            "meet cai_greedy. The MFE comparison is not this claim."
        )
        reason = definition
    return {
        "gate": "GATE-BEN-001",
        "definition": "cai_only",
        "claimed": len(rows) == 50 and not behind,
        "claim": (
            f"CAI-only front, {claimed_n} of {len(rows)} E. coli K-12 proteins{exclusion}, "
            f"cai_greedy baseline, DNA Design Suite meets or exceeds baseline on all {claimed_n}."
        ),
        "n_proteins": len(rows),
        "n_claimed": claimed_n,
        "comparable": True,
        "selected_cds_sha256": selected_cds_sha256,
        "beats_native": beats_native,
        "codonforge_ge_greedy": claimed_n,
        "behind_greedy": behind,
        "degenerate_cases": {"excluded": behind, "explanation": explanation},
        "reason": reason,
    }


def check_cai_gate(gate: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    """Refuse a claim that covers a protein the recorded scores do not support.

    A claim of all 50, or an exclusion list that leaves out a protein whose
    DNA Design Suite value is below `cai_greedy`, is an error. Equality counts as
    meeting the baseline.
    """
    if gate.get("claimed") is not True:
        raise BenchError("GATE-BEN-001 is not claimed")
    excluded = set((gate.get("degenerate_cases") or {}).get("excluded") or [])
    short = [row["gene"] for row in rows if row["codonforge_cai"] + 1e-12 < row["cai_greedy_cai"]]
    uncovered = [gene for gene in short if gene not in excluded]
    if uncovered:
        raise BenchError("claimed proteins underperform cai_greedy: " + ", ".join(uncovered))
    spare = [gene for gene in excluded if gene not in short]
    if spare:
        raise BenchError("excluded proteins that meet cai_greedy: " + ", ".join(spare))
    claimed_n = len(rows) - len(excluded)
    claim = str(gate.get("claim") or "")
    if f"{claimed_n} of {len(rows)}" not in claim:
        raise BenchError(f"the claim does not say {claimed_n} of {len(rows)}; it says {claim!r}")
    if claimed_n != gate.get("n_claimed"):
        raise BenchError(f"n_claimed is {gate.get('n_claimed')}, the rows support {claimed_n}")


def run_cai_only(directory: Path, cds_path: Path | None = None) -> dict[str, Any]:
    """CAI-only redefinition of GATE-BEN-001. No structure axis, no external binary.

    Path A (a competitor that emits windowed accessibility) did not clear the
    three conditions, so this is the comparison the gate can actually run.
    Native, CAI-greedy and one DNA Design Suite design are scored by
    `score_sequences` / the design front, which is one implementation of
    `usage.cai`. The weights are the ribosomal-protein relative adaptiveness.
    """
    from dna_design_suite.codon.orchestrator.design import DesignRequest, design, score_sequences
    from dna_design_suite.modelreg import ArtifactStore

    panel = load_panel()
    genes = select_genes(panel, cds_path)
    usage = ribosomal_usage(cds_path)
    profile = {
        "codon_usage": {
            "subset": "user_supplied",
            "table": {codon: {"relative_adaptiveness": weight} for codon, weight in usage.items()},
        }
    }
    constraints = {
        "schema_version": "1.0.0",
        "constraint_set_id": "cst_cai_only",
        "constraints": [],
    }
    directory.mkdir(parents=True, exist_ok=True)

    def one(gene: dict[str, Any]) -> dict[str, Any]:
        # A private store per gene. Publishing into one store from several
        # threads would race the same digest path.
        store = ArtifactStore(directory / "stores" / str(gene["gene"]))
        cds = gene["cds"]
        greedy = cai_greedy(cds, panel["translation_table"], usage)
        request = DesignRequest(
            host="escherichia_coli:511145/MG1655",
            cds=cds,
            host_profile=profile,
            objectives=["usage.cai"],
            constraints=constraints,
            translation_table=11,
            budget="fast",
            seed=7,
            generations=4,
            segment=2,
            host_digest="sha256:" + "ab" * 32,
        )
        scored = score_sequences(request, [cds, greedy])
        result = design(request, store)
        if result.job.state != "succeeded" or not result.front.get("candidates"):
            raise BenchError(f"{gene['gene']}: CAI-only design did not return a candidate")
        vector = result.front["candidates"][0].get("objective_vector") or [None]
        codonforge = vector[0]
        native_cai = scored[0]["usage.cai"]
        greedy_cai = scored[1]["usage.cai"]
        if codonforge is None or native_cai is None or greedy_cai is None:
            raise BenchError(f"{gene['gene']}: a CAI value is missing")
        zero_weight = [
            {"index": index // 3, "codon": cds[index : index + 3]}
            for index in range(0, len(cds), 3)
            if cds[index : index + 3] not in ("TAA", "TAG", "TGA")
            and usage.get(cds[index : index + 3], 0.0) == 0.0
        ]
        return {
            "gene": gene["gene"],
            "protein_id": gene["protein_id"],
            "codons": gene["codons"],
            "native_cai": native_cai,
            "cai_greedy_cai": greedy_cai,
            "codonforge_cai": codonforge,
            "zero_weight_codons": zero_weight,
        }

    from concurrent.futures import ThreadPoolExecutor

    # `map` keeps gene order. The designs do not fold, so the threads do
    # not share ViennaRNA's process-global temperature.
    with ThreadPoolExecutor(max_workers=min(8, len(genes) or 1)) as pool:
        rows = list(pool.map(one, genes))
    gate = cai_gate_document(rows, selected_cds_sha256=panel["selected_cds_sha256"])
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "cai_only.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    (directory / "gates.json").write_text(json.dumps([gate], indent=2), encoding="utf-8")
    return gate
