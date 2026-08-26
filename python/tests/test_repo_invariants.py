"""Repository invariants, run as tests so `pytest` alone is a complete gate.

Ticket P1-02. These wrap the three checker tools so that a contributor who runs
`pytest` locally gets the same answer CI does, rather than discovering a
layering or schema-sync failure only after pushing.

The invariants themselves are not new here — they are the mechanical
enforcement of decisions made in the planning package:

  * `validate_planning_package.py` — schemas, cross-references, dependency
    graph acyclicity, traceability (BRIEF-2.1-b, BRIEF-5)
  * `check_module_layering.py`     — the real source tree matches modules.yaml,
    and pyo3 stays confined to one crate (ADR-0001)
  * `check_schema_sync.py`         — schemas/v1 matches the normative planning
    copy, and neither the identifier patterns nor the error registry has
    drifted from it (P1-03)
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

CHECKERS = [
    pytest.param(
        ["planning/tools/validate_planning_package.py"],
        id="planning-package",
    ),
    pytest.param(
        ["tools/check_module_layering.py", "--check"],
        id="module-layering",
    ),
    pytest.param(
        ["tools/check_schema_sync.py", "--check"],
        id="schema-sync",
    ),
]


@pytest.mark.parametrize("argv", CHECKERS)
def test_repository_invariant(argv: list[str]) -> None:
    result = subprocess.run(
        [sys.executable, *argv],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"{' '.join(argv)} failed\n--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )


def test_no_forbidden_structure_dependency() -> None:
    """AMD-009 / RISK-DEP-001 / RISK-DEP-002, enforced rather than documented.

    ViennaRNA, LinearFold and LinearDesign forbid redistribution for a fee, and
    LinearDesign is patent-encumbered, so none of them may appear in the
    dependency tree — not as an optional extra, not behind a feature flag. This
    test reads the manifests directly, because the failure mode it guards
    against is somebody adding one "just for benchmarking".
    """
    forbidden = ("viennarna", "vienna-rna", "linearfold", "lineardesign", "rnastructure")
    manifests = [ROOT / "pyproject.toml", ROOT / "Cargo.toml"]
    manifests += sorted(ROOT.glob("crates/*/Cargo.toml"))

    offenders = []
    for manifest in manifests:
        text = manifest.read_text(encoding="utf-8").lower()
        # Skip comment lines: the policy is *stated* in comments in several
        # manifests, and matching those would make the test unfalsifiable.
        body = "\n".join(line for line in text.splitlines() if not line.strip().startswith("#"))
        for name in forbidden:
            if name in body:
                offenders.append(f"{manifest.relative_to(ROOT)} references {name!r}")

    assert not offenders, (
        "forbidden structure-prediction dependency found; see AMD-009:\n  " + "\n  ".join(offenders)
    )
