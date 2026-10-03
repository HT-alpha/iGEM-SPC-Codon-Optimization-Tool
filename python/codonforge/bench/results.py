"""Result directory. `DEP-BEN-001`.

The tidy table is JSON, not Parquet: this package does not depend on a
Parquet library, and the bench target is not allowed to grow one just to
match a filename. `report.html` is filled from that JSON. It does not call
the hypervolume function.
"""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

from codonforge.bench.reference import BenchError, hypervolume


def write_result(
    directory: Path,
    *,
    manifest: dict[str, Any],
    front: dict[str, Any],
    lockfile: dict[str, Any],
    senses: list[str],
    points: list[list[float]],
) -> dict[str, Any]:
    """Write the artifacts and the hypervolume computed from them."""
    directory.mkdir(parents=True, exist_ok=True)
    reference = manifest["reference_point"]
    value = hypervolume(points, reference, senses)
    metric = {
        "suite": manifest["suite"],
        "arm": manifest["arm"],
        "protein": manifest["protein"],
        "seed": manifest["seed"],
        "metric": "hypervolume",
        "value": value,
    }
    (directory / "manifest.json").write_text(json.dumps(manifest, sort_keys=True), encoding="utf-8")
    (directory / "front.json").write_text(json.dumps(front, sort_keys=True), encoding="utf-8")
    (directory / "lockfile.json").write_text(json.dumps(lockfile, sort_keys=True), encoding="utf-8")
    (directory / "metrics.json").write_text(json.dumps([metric], sort_keys=True), encoding="utf-8")
    # No threshold was supplied with this result, so no gate is recorded as
    # passed. An empty list is the absence of a claim.
    (directory / "gates.json").write_text("[]\n", encoding="utf-8")
    (directory / "report.html").write_text(render_report([metric]), encoding="utf-8")
    return metric


def check_reproduces(directory: Path, senses: list[str], points: list[list[float]]) -> float:
    """Recompute the stored hypervolume from the published reference point.

    Raises when the stored metric, the front's lockfile digest and the
    lockfile document disagree. `TST-BEN-001`.
    """
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    front = json.loads((directory / "front.json").read_text(encoding="utf-8"))
    lockfile = json.loads((directory / "lockfile.json").read_text(encoding="utf-8"))
    metrics = json.loads((directory / "metrics.json").read_text(encoding="utf-8"))
    if lockfile.get("digest") != manifest.get("lockfile_digest"):
        raise BenchError("manifest lockfile digest does not match the lockfile document")
    if front.get("landscape_lockfile_digest") != manifest.get("lockfile_digest"):
        raise BenchError("front lockfile digest does not match the manifest")
    stored = next(row["value"] for row in metrics if row["metric"] == "hypervolume")
    recomputed = hypervolume(points, manifest["reference_point"], senses)
    if recomputed != stored:
        raise BenchError(
            f"recomputed hypervolume {recomputed} does not match the stored metric {stored}"
        )
    return recomputed


def render_report(metrics: list[dict[str, Any]]) -> str:
    """HTML whose numbers are the metric rows. No indicator is evaluated here."""
    rows = []
    for row in metrics:
        cells = "".join(
            f"<td>{html.escape(str(row[key]))}</td>"
            for key in ("suite", "arm", "protein", "seed", "metric", "value")
        )
        rows.append(f"<tr>{cells}</tr>")
    body = "\n".join(rows)
    return (
        "<!DOCTYPE html><html><body><table>\n"
        "<tr><th>suite</th><th>arm</th><th>protein</th><th>seed</th>"
        "<th>metric</th><th>value</th></tr>\n"
        f"{body}\n</table></body></html>\n"
    )
