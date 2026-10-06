"""The optional GSE53767 dwell head.

The head is a separate axis. A job that does not name it never reads this
file. A job that names it without a fitted head is `CF-LAND-001`. A zero
deviation is not substituted for a missing fit.

The fit is not in this tree. Li, Burkhardt, Gross and Weissman, Cell 2014
(GEO GSE53767) do not state a 5'-end A-site offset. Their footprints were
size-selected from about 15 to 45 nt, so one offset does not place the
A site. The assignment rule in Li, Oh and Weissman, Nature 2012
(doi:10.1038/nature10965, GEO GSE35641, a different series) puts the A site
at the footprint center — residues at least 12 nt from either end — and
does not state +15. The GEO supplementary files were not reachable from the
environment that would have extracted the counts.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from codonforge.epistasis.artifact import ArtifactError
from codonforge.modelreg.manifest import digest_of

DWELL_TERM = "epistasis.dwell_deviation"
ARTIFACT_NAME = "dwell.json"

ABSENCE_REASON = (
    "epistasis.dwell_deviation was requested but the fitted dwell head is absent. "
    "GSE53767 does not state a 5'-end A-site offset (Li et al., Cell 2014; "
    "footprints ~15–45 nt). Li, Oh and Weissman, Nature 2012 assigns footprint "
    "centers on GSE35641, not a +15 nt offset, and that rule was not applied here. "
    "No dwell deviation of zero was substituted."
)


def dwell_path() -> Path:
    return Path(__file__).resolve().parents[3] / "benchmarks" / "epistasis" / ARTIFACT_NAME


def verify_dwell(document: dict[str, Any]) -> None:
    """The digest covers the document. A replaced file fails here."""
    recorded = document.get("digest")
    recomputed = digest_of(document)
    if recomputed != recorded:
        raise ArtifactError(f"dwell digest {recorded} does not match its content {recomputed}")


def dwell_status(path: Path | None = None) -> dict[str, Any]:
    """What calibration may say. Absent is not a fitted deviation of zero."""
    artifact = path or dwell_path()
    if not artifact.is_file():
        return {
            "state": "absent",
            "term": DWELL_TERM,
            "accession": "GSE53767",
            "a_site_offset_nt": None,
            "normalization": None,
            "reason": ABSENCE_REASON,
        }
    document = json.loads(artifact.read_text(encoding="utf-8"))
    verify_dwell(document)
    return {
        "state": document.get("state", "absent"),
        "term": DWELL_TERM,
        "accession": document.get("accession", "GSE53767"),
        "a_site_offset_nt": document.get("a_site_offset_nt"),
        "normalization": document.get("normalization"),
        "reason": document.get("reason", ABSENCE_REASON),
        "digest": document.get("digest"),
    }
