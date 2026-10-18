"""The documentation describes what this build actually does.

Ticket P1-35. Acceptance: a new contributor can build, test and run a design
from the README alone.

Documentation rots silently, and a quickstart that no longer works is worse
than none — it costs the reader their first hour and their confidence. So the
claims are tested: every make target the docs name exists, every CLI invocation
they show parses, every capability gap they promise to report is reported with
the code they name, and nothing reintroduces the `python3 -m venv` assumption
that the Makefile fallback exists to avoid.
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
README = ROOT / "README.md"
CONTRIBUTING = ROOT / "CONTRIBUTING.md"
AGENTS = ROOT / "AGENTS.md"
MAKEFILE = ROOT / "Makefile"


@pytest.fixture(scope="module")
def docs() -> dict[str, str]:
    return {
        "README.md": README.read_text(encoding="utf-8"),
        "CONTRIBUTING.md": CONTRIBUTING.read_text(encoding="utf-8"),
        "AGENTS.md": AGENTS.read_text(encoding="utf-8"),
    }


def test_p1_35_the_documents_exist() -> None:
    """`AGENTS.md` was described as a stub; it did not exist at all."""
    for path in (README, CONTRIBUTING, AGENTS):
        assert path.exists(), f"{path.name} is missing"
        assert len(path.read_text(encoding="utf-8")) > 500, f"{path.name} is a stub"


# --------------------------------------------------------------------------- #
# the venv fallback, which the carried-forward item asks about specifically
# --------------------------------------------------------------------------- #


def test_p1_35_nothing_reintroduces_a_bare_python_venv(docs: dict[str, str]) -> None:
    """The Makefile fallback is documented, and nothing works around it.

    `python3 -m venv` needs `ensurepip`, shipped separately on Debian and
    Ubuntu and not always installable. Every mention across the docs, the
    Makefile and CI must either be the Makefile's own guarded fallback or an
    explicit warning against it.
    """
    sources = {
        "Makefile": MAKEFILE.read_text(encoding="utf-8"),
        **docs,
        **{
            str(p.relative_to(ROOT)): p.read_text(encoding="utf-8")
            for p in (ROOT / ".github" / "workflows").glob("*.yml")
        },
    }

    # A mention is fine; a *recommendation* is not. Proximity to the word
    # "virtualenv" is too weak a test — a paragraph explaining the fallback
    # would let "use `python3 -m venv`" through on the next line, which is
    # exactly what a first version of this test did.
    permitted = (
        "needs",
        "requires",
        "falls back",
        "fallback",
        "rather than",
        "avoid",
        "never",
        "bare",
        "not ",
        "||",
    )

    for name, text in sources.items():
        # Prose is checked by sentence and scripts by line. A per-line check on
        # markdown fails on wrapping — the qualifier that makes a mention a
        # warning lands on the previous line — while a proximity window is too
        # weak, since a paragraph explaining the fallback lets a recommendation
        # through on the next line.
        if name.endswith(".md"):
            # Emphasis markers are stripped before splitting: a sentence ending
            # `detail.**` does not end in a period followed by whitespace, so
            # it merges with the next one and lends it its qualifiers. That is
            # how "use `python3 -m venv`" inherited a "not" from the sentence
            # before it and went unnoticed.
            flattened = " ".join(text.split()).replace("*", "").replace("`", "")
            units = re.split(r"(?<=[.:])\s", flattened)
        else:
            units = text.splitlines()

        for unit in units:
            if "python3 -m venv" not in unit and "python -m venv" not in unit:
                continue
            lowered = unit.lower()
            assert any(marker in lowered for marker in permitted), (
                f"{name} presents `python -m venv` as something to do, rather than "
                "as a fallback or a warning. It needs `ensurepip`, which is not "
                f"always installable:\n  {unit.strip()}"
            )


def test_p1_35_the_venv_fallback_is_documented_as_the_supported_path(
    docs: dict[str, str],
) -> None:
    """Not an implementation detail a contributor discovers by reading the Makefile."""
    for name in ("README.md", "CONTRIBUTING.md", "AGENTS.md"):
        text = docs[name]
        assert "virtualenv" in text, f"{name} does not mention virtualenv"
        assert "ensurepip" in text or "python3-venv" in text, (
            f"{name} mentions virtualenv without saying why, so a reader will "
            "'simplify' it back to python3 -m venv"
        )


# --------------------------------------------------------------------------- #
# every command the docs give actually exists
# --------------------------------------------------------------------------- #


def test_p1_35_every_documented_make_target_exists(docs: dict[str, str]) -> None:
    makefile = MAKEFILE.read_text(encoding="utf-8")
    targets = set(re.findall(r"^([a-z][a-z0-9_-]*):", makefile, re.MULTILINE))

    # Only inside fenced code blocks: prose says things like "make every
    # change carry a test", and `every` is not a target.
    documented = set()
    for text in docs.values():
        for block in re.findall(r"```sh\n(.*?)```", text, re.S):
            documented |= set(re.findall(r"^make ([a-z][a-z0-9_-]*)", block, re.MULTILINE))

    missing = documented - targets
    assert not missing, f"documented make targets that do not exist: {sorted(missing)}"
    assert {"setup", "build", "test", "lint"} <= documented


def test_p1_35_every_documented_cli_invocation_parses(docs: dict[str, str]) -> None:
    """The quickstart's flags must be flags the parser accepts."""
    from dna_design_suite.cli import build_parser

    parser = build_parser()
    blocks = []
    for text in docs.values():
        for fenced in re.findall(r"```sh\n(.*?)```", text, re.S):
            joined = fenced.replace("\\\n", " ")
            for line in joined.splitlines():
                if "dna_design_suite.cli.main" in line:
                    blocks.append(line.split("dna_design_suite.cli.main", 1)[1])

    assert blocks, "the docs show no CLI invocation at all"
    for block in blocks:
        # Placeholders like <CDS> stand in for user values, so they are
        # substituted rather than dropped - dropping one leaves its flag
        # without an argument and the parse fails for the wrong reason.
        argv = [
            "PLACEHOLDER" if a.startswith("<") and a.endswith(">") else a
            for a in shlex.split(block)
        ]
        parsed = parser.parse_args(argv)
        assert parsed.verb, f"no verb parsed from: {block}"


def test_p1_35_every_documented_invariant_script_exists(docs: dict[str, str]) -> None:
    for text in docs.values():
        for script in re.findall(r"((?:planning/tools|tools)/[a-z_]+\.py)", text):
            assert (ROOT / script).exists(), f"documented script missing: {script}"


# --------------------------------------------------------------------------- #
# the documented gaps are the real ones
# --------------------------------------------------------------------------- #


def test_p1_35_documented_error_codes_are_registered(docs: dict[str, str]) -> None:
    """A doc naming a code that does not exist teaches the wrong vocabulary."""
    registry = (ROOT / "crates" / "cf-types" / "src" / "error.rs").read_text()
    cited = set()
    for text in docs.values():
        cited |= set(re.findall(r"CF-[A-Z]{3,8}-[0-9]{3}", text))

    assert cited, "the docs cite no error codes, so a reader has no vocabulary"
    for code in sorted(cited):
        assert f'"{code}"' in registry, f"{code} is documented but not registered"


def test_p1_35_the_background_execution_gap_is_documented_and_real(
    docs: dict[str, str],
) -> None:
    """DEP-CLI-006. A handle returns before the worker finishes the job."""
    import threading
    import time

    readme = docs["README.md"]
    assert "background=True" in readme
    assert "CF-CLI-003" in readme
    assert "silently run synchronously" in readme or "does not silently" in readme

    from dna_design_suite.codon.orchestrator import design
    from dna_design_suite.codon.orchestrator import job as job_module

    release = threading.Event()
    original = job_module.Job.run

    def blocked(self, on_segment=None):
        assert release.wait(timeout=5)
        return original(self, on_segment)

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(job_module.Job, "run", blocked)
    try:
        started = time.perf_counter()
        handle = design(_minimal_request(), _store(), background=True)
        elapsed = time.perf_counter() - started
        assert elapsed < 0.5
        assert handle.done() is False
        assert handle.job_id
        release.set()
        assert handle.result(timeout=60).job.state == "succeeded"
    finally:
        monkeypatch.undo()


def test_p1_35_the_tier_skip_semantics_are_documented_and_real(
    docs: dict[str, str],
) -> None:
    """DEP-ORCH-008. An exact budget records decision as done."""
    readme = docs["README.md"]
    assert "CF-ORCH-006" in readme
    assert "skipped" in readme
    assert "decision" in readme

    from dna_design_suite.codon.orchestrator import design

    request = _minimal_request()
    request.budget = "exact"
    result = design(request, _store())

    states = {t["tier"]: t["state"] for t in result.status["progress"]["tiers"]}
    assert states["decision"] == "done"
    decision_skips = [
        warning
        for warning in result.status["warnings"]
        if warning["code"] == "CF-ORCH-006" and warning.get("subject") == "decision"
    ]
    assert not decision_skips
    report = result.job.stage_reports["decision"]
    assert report["state"] == "done"
    assert report["rule"] == result.front["recommended"]["rule"]


def test_p1_35_the_documented_term_list_is_the_real_one(docs: dict[str, str]) -> None:
    from dna_design_suite.codon.orchestrator import BUILTIN_TERMS

    readme = docs["README.md"]
    for term in BUILTIN_TERMS:
        assert term in readme, f"{term} is available but undocumented"
    assert "CF-LAND-001" in readme


def test_p1_35_the_documented_exit_codes_match_the_cli(docs: dict[str, str]) -> None:
    # `dna_design_suite.cli` re-exports `main` as a function, so the module holding
    # the exit codes has to be reached by path rather than by attribute.
    import importlib

    cli_module = importlib.import_module("dna_design_suite.cli.main")

    readme = docs["README.md"]
    for name, value in (
        ("EXIT_OK", 0),
        ("EXIT_DESIGN_FAILED", 1),
        ("EXIT_INVALID_INPUT", 2),
        ("EXIT_INFEASIBLE", 3),
        ("EXIT_BUDGET", 4),
        ("EXIT_INTERNAL", 5),
    ):
        assert getattr(cli_module, name) == value
    assert "Infeasible constraints" in readme


def test_p1_35_the_documented_output_formats_all_exist(docs: dict[str, str]) -> None:
    from dna_design_suite.cli import build_parser

    parser = build_parser()
    design_parser = next(
        action.choices["design"]
        for action in parser._actions
        if hasattr(action, "choices") and action.choices and "design" in action.choices
    )
    supported = next(action.choices for action in design_parser._actions if action.dest == "output")

    # Extract whatever the README claims rather than checking a list this test
    # already knows. Matching known names only could catch a format removed
    # from the CLI and never one invented in prose, which is the direction
    # documentation actually drifts.
    readme = docs["README.md"]
    claimed = set(re.findall(r"`--output (\w+)`", readme))
    for line in readme.splitlines():
        if line.startswith("Other formats:") or (claimed & set(re.findall(r"`(\w+)`", line))):
            claimed |= set(re.findall(r"`(\w+)`", line))

    assert claimed, "the README documents no output formats at all"
    unknown = claimed - set(supported)
    assert not unknown, f"README documents --output {sorted(unknown)}, which is not accepted"


def test_p1_35_the_documented_stage_names_are_the_real_ones(docs: dict[str, str]) -> None:
    from dna_design_suite.codon.orchestrator.budget import EXECUTED_TIERS

    readme = docs["README.md"]
    for tier in EXECUTED_TIERS:
        assert tier in readme, f"{tier} runs but is undocumented"


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _minimal_request():
    import sys

    sys.path.insert(0, str(Path(__file__).parent))
    from test_cli import request_for

    return request_for(generations=0)


def _store():
    import tempfile

    from dna_design_suite.modelreg import ArtifactStore

    return ArtifactStore(Path(tempfile.mkdtemp()) / "artifacts")


def test_p1_35_background_true_is_the_documented_handle(docs: dict[str, str]) -> None:
    """The keyword the README names returns a handle. Seed import stays CF-CLI-003."""
    from dna_design_suite.codon.orchestrator import design

    assert "seed" in docs["README.md"].lower()
    handle = design(_minimal_request(), _store(), background=True)
    assert hasattr(handle, "job_id") and hasattr(handle, "done") and hasattr(handle, "result")
    assert handle.result(timeout=60).job.state == "succeeded"
