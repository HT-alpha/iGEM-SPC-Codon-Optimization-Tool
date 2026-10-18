"""Fit the context model on NC_000913.3, with the panel held out.

The 13 panel genes are inside the 4140 CDS. They are removed before the
fit. Scoring them afterwards is the held-out number, not a second use of
the training set.
"""

from __future__ import annotations

import json
import statistics
from pathlib import Path
from typing import Any

from dna_design_suite import _core
from dna_design_suite.codon.bench.corpus import CDS_PATH, PANEL_PATH, iter_cds

# Numerical stabilizers on the squared parameters. Not expression rates.
LAMBDA_H = 1e-3
LAMBDA_J = 1e-3
MAX_ITERATIONS = 40
WINDOWS = (5, 10)


def held_out_split(
    panel_genes: list[str],
    records: list[tuple[str, str]] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Split CDS records. A panel gene never appears in the training list."""
    if records is None:
        records = [(header, sequence) for header, sequence in iter_cds()]
    held: list[dict[str, str]] = []
    train: list[dict[str, str]] = []
    wanted = set(panel_genes)
    for header, sequence in records:
        gene = _gene(header)
        row = {"gene": gene or "", "cds": sequence, "header": header}
        if gene in wanted:
            held.append(row)
        else:
            train.append(row)
    found = {row["gene"] for row in held}
    missing = sorted(wanted - found)
    if missing:
        raise ValueError(f"held-out genes missing from the CDS file: {missing}")
    leaked = [row["gene"] for row in train if row["gene"] in wanted]
    if leaked:
        raise ValueError(f"panel genes remained in the training set: {leaked}")
    return train, held


def _gene(header: str) -> str | None:
    marker = "[gene="
    if marker not in header:
        return None
    return header.split(marker, 1)[1].split("]", 1)[0]


def _call(name: str, payload: dict[str, Any]) -> dict[str, Any]:
    return json.loads(getattr(_core, name)(json.dumps(payload)))


def fit_window(sequences: list[str], window: int) -> dict[str, Any]:
    return _call(
        "fit_potts_context",
        {
            "sequences": sequences,
            "window": window,
            "lambda_h": LAMBDA_H,
            "lambda_j": LAMBDA_J,
            "max_iterations": MAX_ITERATIONS,
        },
    )


def score_model(sequences: list[str], fitted: dict[str, Any]) -> dict[str, Any]:
    return _call(
        "score_potts_context",
        {
            "sequences": sequences,
            "window": fitted["window"],
            "h": fitted["h"],
            "couplings": fitted["couplings"],
        },
    )


def per_gene_nats(held: list[dict[str, str]], fitted: dict[str, Any]) -> list[float]:
    values = []
    for row in held:
        scored = score_model([row["cds"]], fitted)
        values.append(float(scored["nats_per_codon"]))
    return values


def select_window(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Smallest W within one standard error of the best held-out mean.

    The standard error is the gene-level standard error of the best window.
    Thirteen genes is a small sample; the rule is still the one in ADR-0003.
    """
    best = max(rows, key=lambda row: row["held_out_nats_per_codon"])
    genes = best["per_gene_nats"]
    error = statistics.pstdev(genes) / (len(genes) ** 0.5)
    eligible = [
        row
        for row in rows
        if row["held_out_nats_per_codon"] >= best["held_out_nats_per_codon"] - error
    ]
    chosen = min(eligible, key=lambda row: row["window"])
    chosen = dict(chosen)
    chosen["standard_error"] = error
    chosen["selection"] = "one_standard_error"
    return chosen


def fit_panel(directory: Path | None = None) -> dict[str, Any]:
    """Write the artifact. Returns the summary, including both windows."""
    panel = json.loads(PANEL_PATH.read_text(encoding="utf-8"))
    # The fit held out the original 13. The benchmark list is larger, and
    # those extra genes were in the training CDS. Do not hold them out now.
    train, held = held_out_split(list(panel["potts_held_out"]))
    train_cds = [row["cds"] for row in train]
    held_cds = [row["cds"] for row in held]
    baseline = _call(
        "score_codon_usage_baseline",
        {"train": train_cds, "held_out": held_cds},
    )
    rows = []
    for window in WINDOWS:
        fitted = fit_window(train_cds, window)
        scored = score_model(held_cds, fitted)
        genes = per_gene_nats(held, fitted)
        rows.append(
            {
                "window": window,
                "held_out_nats_per_codon": scored["nats_per_codon"],
                "n_positions": scored["n_positions"],
                "train_nats_per_codon": fitted["nats_per_codon"],
                "iterations": fitted["iterations"],
                "gradient_norm": fitted["gradient_norm"],
                "converged": fitted["converged"],
                "per_gene_nats": genes,
                "h": fitted["h"],
                "couplings": fitted["couplings"],
            }
        )
    chosen = select_window(rows)
    summary = {
        "n_cds_file": len(train) + len(held),
        "n_train": len(train),
        "n_held_out": len(held),
        "held_out_genes": [row["gene"] for row in held],
        "baseline_nats_per_codon": baseline["nats_per_codon"],
        "baseline_n_scored": baseline["n_scored"],
        "baseline_n_unseen": baseline["n_unseen"],
        "windows": [
            {
                "window": row["window"],
                "held_out_nats_per_codon": row["held_out_nats_per_codon"],
                "train_nats_per_codon": row["train_nats_per_codon"],
                "iterations": row["iterations"],
                "gradient_norm": row["gradient_norm"],
                "converged": row["converged"],
            }
            for row in rows
        ],
        "selected_window": chosen["window"],
        "selected_held_out_nats_per_codon": chosen["held_out_nats_per_codon"],
        "standard_error": chosen["standard_error"],
    }
    if directory is not None:
        from dna_design_suite.codon.epistasis.artifact import write_artifact

        write_artifact(directory, panel, summary, chosen, rows)
    return summary


def corpus_digest_matches() -> None:
    """The fit reads the same gzip the panel names."""
    import hashlib

    panel = json.loads(PANEL_PATH.read_text(encoding="utf-8"))
    digest = hashlib.sha256(CDS_PATH.read_bytes()).hexdigest()
    if digest != panel["cds_gzip_sha256"]:
        raise ValueError("CDS gzip does not match the panel digest")


if __name__ == "__main__":
    corpus_digest_matches()
    destination = Path(__file__).resolve().parents[4] / "benchmarks" / "epistasis"
    print(json.dumps(fit_panel(destination), indent=2))
