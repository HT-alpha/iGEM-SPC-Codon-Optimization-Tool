#!/usr/bin/env python3
"""Copy the demo files a wheel has to serve into the package tree.

``examples/`` is the copy the README and a source checkout use. An installed
wheel has no repository around it, so ``dna_design_suite.api`` also looks in
``python/dna_design_suite/api/examples``. Those bytes are this script's output.
``python/tests/test_frontend_entry.py`` refuses a drift between the two.
"""

from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "examples"
DEST = ROOT / "python" / "dna_design_suite" / "api" / "examples"
NAMES = (
    "biobrick_rfc10.json",
    "gibson.json",
    "golden_gate.json",
    "kazusa_e_coli_316407.json",
    "mg1655_ribosomal_usage.json",
    "rpsL.cds",
    "rpsL_constraints.json",
    "synthesis.json",
    "BBa_E0040.fa",
)


def stage() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    for name in NAMES:
        shutil.copyfile(SOURCE / name, DEST / name)


if __name__ == "__main__":
    stage()
