"""The optional GSE53767 dwell head.

The head is a separate axis. A job that does not name it never reads this
file. A job that names it without a fitted head is `CF-LAND-001`. A zero
deviation is not substituted for a missing fit.

The fit is not in this tree. No paper checked below states one 5'-end
offset for GSE53767. Li et al., Cell 2014, map 20–42 nt reads by giving the
same score to center residues at least 10 nt from either end. That is not a
single codon. Li, Oh and Weissman, Nature 2012, do the same kind of
center-weighting on a different series (GSE35641, 25–42 nt, 12 nt from
either end) and do not state +15. The 2014 methods say the size selection
differs from that earlier procedure. Ingolia et al., Science 2009, is yeast.
Mohammad, Green and Buskirk, eLife 2019, assign their own libraries to the
3' end. RiboToolkit and RiboCode compute a P-site from the reads in front of
them; neither states an offset for GSE53767. The GEO files were not
reachable, so the counts were not extracted under any rule.
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
    "No 5'-end A-site offset is stated for GSE53767. Li et al., Cell 2014, "
    "doi:10.1016/j.cell.2014.02.033, size-select about 15-45 nt and score the "
    "center residues at least 10 nt from either end of 20-42 nt reads. "
    "Li, Oh and Weissman, Nature 2012, doi:10.1038/nature10965, GSE35641, "
    "score center residues at least 12 nt from either end of 25-42 nt reads "
    "and do not state +15; the 2014 size selection is not that procedure. "
    "Ingolia et al., Science 2009, doi:10.1126/science.1168978, is yeast: "
    "5' ends start 12-13 nt upstream of the start codon (Fig. 2). "
    "Mohammad, Green and Buskirk, eLife 2019, doi:10.7554/eLife.42591, assign "
    "their own 3' ends, about 12 nt from the A site, on a flash-freeze "
    "high-Mg library, not GSE53767. RiboToolkit (Nucleic Acids Res. 2020, "
    "doi:10.1093/nar/gkaa395) defaults to 26-32 nt and a user or plastid "
    "P-site. RiboCode (Nucleic Acids Res. 2018, doi:10.1093/nar/gky179) "
    "infers a P-site per read length and does not name GSE53767. "
    "The GEO files were not reachable. No dwell deviation of zero was substituted."
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
        "checked": document.get("checked", []),
        "what_would_resolve": document.get("what_would_resolve"),
        "digest": document.get("digest"),
    }
