"""The Potts artifact and the check that its digest covers the parameters.

`digest_of` hashes the document with the digest field omitted. `h` and the
couplings sit in that document, so changing one entry changes the digest.
A digest of the accession and the window alone would not.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from dna_design_suite.modelreg.manifest import digest_of

ARTIFACT_NAME = "manifest.json"


class ArtifactError(ValueError):
    """The file on disk is not the document its digest names."""


def verify_artifact(document: dict[str, Any]) -> None:
    """Recompute the digest over the whole document, parameters included."""
    recorded = document.get("digest")
    recomputed = digest_of(document)
    if recomputed != recorded:
        raise ArtifactError(f"artifact digest {recorded} does not match its content {recomputed}")


def write_artifact(
    directory: Path,
    panel: dict[str, Any],
    summary: dict[str, Any],
    chosen: dict[str, Any],
    rows: list[dict[str, Any]],
) -> dict[str, Any]:
    directory.mkdir(parents=True, exist_ok=True)
    document = _document(panel, summary, chosen, rows)
    document["digest"] = digest_of(document)
    verify_artifact(document)
    path = directory / ARTIFACT_NAME
    path.write_text(json.dumps(document), encoding="utf-8")
    # The summary is what a reader prints. The manifest is what the digest covers.
    (directory / "summary.json").write_text(
        json.dumps(_public_summary(document), indent=2),
        encoding="utf-8",
    )
    return document


def load_calibration(path: Path | None = None) -> dict[str, Any]:
    """What the result page is allowed to say.

    A missing file is absent. It is not a nats value of zero.
    """
    artifact = path or _default_path()
    if not artifact.is_file():
        return {"state": "absent"}
    document = json.loads(artifact.read_text(encoding="utf-8"))
    verify_artifact(document)
    parameters = document["parameters"]
    return {
        "state": "fitted",
        "window": parameters["selected_window"],
        "held_out_nats_per_codon": parameters["selected_held_out_nats_per_codon"],
        "baseline_nats_per_codon": parameters["baseline_nats_per_codon"],
        "iteration_count": parameters["iteration_count"],
        "final_gradient_norm": parameters["final_gradient_norm"],
        "converged": parameters["converged"],
        "n_train": parameters["n_train"],
        "n_held_out": parameters["n_held_out"],
        "digest": document["digest"],
    }


def artifact_path() -> Path:
    """The shipped Potts manifest, whether or not the file is present."""
    return _default_path()


def _default_path() -> Path:
    return Path(__file__).resolve().parents[4] / "benchmarks" / "epistasis" / ARTIFACT_NAME


def _public_summary(document: dict[str, Any]) -> dict[str, Any]:
    parameters = document["parameters"]
    return {
        "state": "fitted",
        "window": parameters["selected_window"],
        "held_out_nats_per_codon": parameters["selected_held_out_nats_per_codon"],
        "baseline_nats_per_codon": parameters["baseline_nats_per_codon"],
        "windows": parameters["comparison"],
        "iteration_count": parameters["iteration_count"],
        "final_gradient_norm": parameters["final_gradient_norm"],
        "converged": parameters["converged"],
        "digest": document["digest"],
        "n_train": parameters["n_train"],
        "n_held_out": parameters["n_held_out"],
    }


def _document(
    panel: dict[str, Any],
    summary: dict[str, Any],
    chosen: dict[str, Any],
    rows: list[dict[str, Any]],
) -> dict[str, Any]:
    import subprocess

    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=False,
        capture_output=True,
        text=True,
    )
    sha = completed.stdout.strip()
    if len(sha) != 40:
        sha = "0" * 40
    delta = chosen["held_out_nats_per_codon"] - summary["baseline_nats_per_codon"]
    return {
        "schema_version": "1.0.0",
        "artifact_id": "cf:epistasis_potts:escherichia_coli:511145/MG1655:0.1.0",
        "kind": "epistasis_potts",
        "host_id": "escherichia_coli:511145/MG1655",
        "version": "0.1.0",
        "contract_major": 1,
        "provides_terms": ["epistasis.context_naturalness"],
        "created_at": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "code_version": {
            "package": "codonforge",
            "version": "0.1.0",
            "git_sha": sha,
            "dirty": True,
        },
        "training": {
            "datasets": [
                {
                    "dataset_id": "ds_nc000913_3_cds",
                    "digest": "sha256:" + panel["cds_gzip_sha256"],
                    "accession": panel["accession"],
                    "n_records": summary["n_cds_file"],
                }
            ],
            "split_strategy": "by_gene",
            "hyperparameters": {
                "windows": [5, 10],
                "lambda_h": 1e-3,
                "lambda_j": 1e-3,
                "mu": 0,
                "optimizer": "lbfgs",
            },
            "n_parameters": 64 + 4096 * int(chosen["window"]),
        },
        "evaluation": {
            "metrics": {
                "held_out_nats_per_codon": chosen["held_out_nats_per_codon"],
                "codon_usage_baseline_nats_per_codon": summary["baseline_nats_per_codon"],
            },
            "baseline_comparison": {
                "baseline_name": "codon_usage_only",
                "metric": "held_out_nats_per_codon",
                "delta": delta,
            },
            "gates": [],
        },
        "isolation": {"scope": "public", "tenant_id": None},
        "lifecycle": {"status": "draft"},
        "citation": (
            "Sequences: NCBI NC_000913.3 (GCF_000005845.2). "
            "Named Ribo-seq dataset, not consumed by this fit: "
            "Li GW, Burkhardt D, Gross C, Weissman JS. Cell. 2014;157:624-635. "
            "GEO GSE53767."
        ),
        "parameters": {
            "assembly": panel["assembly"],
            "accession": panel["accession"],
            "n_train": summary["n_train"],
            "n_held_out": summary["n_held_out"],
            "held_out_genes": summary["held_out_genes"],
            "selected_window": chosen["window"],
            "selected_held_out_nats_per_codon": chosen["held_out_nats_per_codon"],
            "baseline_nats_per_codon": summary["baseline_nats_per_codon"],
            "baseline_n_unseen": summary["baseline_n_unseen"],
            "standard_error": chosen["standard_error"],
            "comparison": summary["windows"],
            "iteration_count": chosen["iterations"],
            "final_gradient_norm": chosen["gradient_norm"],
            "converged": chosen["converged"],
            "models": [
                {
                    "window": row["window"],
                    "held_out_nats_per_codon": row["held_out_nats_per_codon"],
                    "iterations": row["iterations"],
                    "converged": row["converged"],
                    "gradient_norm": row["gradient_norm"],
                    "h": row["h"],
                    "couplings": row["couplings"],
                }
                for row in rows
            ],
        },
    }
