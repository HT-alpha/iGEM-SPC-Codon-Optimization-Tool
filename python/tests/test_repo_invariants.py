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


def test_the_compiled_extension_is_not_older_than_its_sources() -> None:
    """`pytest` run directly can silently test a stale `_core`.

    `make test` depends on `build`, so the supported path rebuilds first. Run
    `pytest` on its own after editing Rust and the extension is whatever was
    last installed — the Python tests then exercise old behaviour and fail, or
    worse pass, for reasons nothing in the output explains. That cost a long
    diagnosis once: a test failure was chased into a Python comprehension that
    turned out to be correct, and the real cause was an extension built before
    the change under test.
    """
    import pathlib

    from codonforge import _core

    root = pathlib.Path(__file__).resolve().parents[2]
    extension = pathlib.Path(_core.__file__)
    assert extension.exists(), "the compiled extension is missing; run `make build`"

    newest = max(
        (
            path.stat().st_mtime
            for path in root.joinpath("crates").rglob("*.rs")
            if "target" not in path.parts
        ),
        default=0.0,
    )
    assert extension.stat().st_mtime >= newest, (
        f"{extension.name} is older than the Rust sources it was built from. "
        "Run `make build` (or `make test`, which does it for you) — running "
        "pytest directly after editing Rust tests the previous build."
    )


def test_amd_013_no_grammar_intersection_reaches_the_search_layer():
    """AMD-013: the LinearDesign exclusion is a claim limitation, not a vibe.

    US 2022/0269951 A1 and EP 4 030 438 A1 both claim, as a required element of
    every independent claim, intersecting a context-free grammar onto a DFA of
    synonymous codon choices and searching the intersected grammar. That element
    is what CodonForge must not practise; the general idea of optimizing codons
    against structure is not what is claimed.

    A prohibition that nothing checks is a prohibition that erodes, so this test
    looks for the machinery such a method would need. It is deliberately keyed
    to grammar-parsing vocabulary rather than to the word "LinearDesign": the
    risk AMD-013 names is a reimplementation under a different name, which a
    name-based check would miss entirely.
    """
    import pathlib
    import re

    root = pathlib.Path(__file__).resolve().parents[2]

    # Vocabulary specific to CFG-over-automaton parsing. "context-free" is
    # excluded on purpose - the codebase uses it in the unrelated sense of a
    # window of 1 in the epistasis model, and a check that fired on that would
    # be noise rather than signal.
    forbidden = re.compile(
        r"\b(scfg|stochastic context.free|nonterminal|bifurcation rule|"
        r"inside.outside|cyk\b|grammar intersection|intersected (cfg|grammar))",
        re.IGNORECASE,
    )

    offenders = []
    for path in list(root.joinpath("crates").rglob("*.rs")) + list(
        root.joinpath("python").rglob("*.py")
    ):
        if "target" in path.parts or path.name == "test_repo_invariants.py":
            continue
        for number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            if forbidden.search(line):
                offenders.append(f"{path.relative_to(root)}:{number}: {line.strip()}")

    assert not offenders, (
        "grammar-parsing machinery appeared in the tree. AMD-013 excludes the "
        "LinearDesign method - a CFG intersected onto a codon DFA - and that "
        "exclusion survives a rename or a from-scratch reimplementation, because "
        "patent claims attach to the method rather than to the code's provenance. "
        "If this is a false positive, narrow the pattern deliberately and say why:\n  "
        + "\n  ".join(offenders)
    )
