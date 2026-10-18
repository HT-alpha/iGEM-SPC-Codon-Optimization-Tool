"""CI and wheel workflows spend Actions minutes only where the gate needs them.

The defect this guards: a push, or a docs-only pull request, starts the
expensive jobs and the suite still looks green. The classifier and the
workflow files are both asserted, because a correct function that the
workflow does not call saves nothing.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
CI = ROOT / ".github" / "workflows" / "ci.yml"
WHEELS = ROOT / ".github" / "workflows" / "wheels.yml"

PUBLISH_PLATFORMS = {
    ("ubuntu-latest", "x86_64"),
    ("ubuntu-latest", "aarch64"),
    ("macos-latest", "x86_64"),
    ("macos-latest", "aarch64"),
    ("windows-latest", "x64"),
}


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def paths():
    return _load("ci_path_filter", ROOT / "tools" / "ci_path_filter.py")


@pytest.fixture(scope="module")
def wheels():
    return _load("wheel_targets", ROOT / "tools" / "wheel_targets.py")


def _pairs(platforms: list[dict[str, str]]) -> set[tuple[str, str]]:
    return {(item["os"], item["target"]) for item in platforms}


def test_benchmark_doc_skips_the_heavy_gate(paths) -> None:
    assert not paths.requires_heavy(["docs/benchmarks/v2.9.2-profiles.md"])


def test_root_readme_skips_the_heavy_gate(paths) -> None:
    assert not paths.requires_heavy(["README.md", "CONTRIBUTING.md", "LICENSE"])


def test_python_and_rust_run_the_heavy_gate(paths) -> None:
    assert paths.requires_heavy(["python/dna_design_suite/cli/main.py"])
    assert paths.requires_heavy(["crates/cf-search/src/lib.rs"])


def test_mixed_docs_and_code_run_the_heavy_gate(paths) -> None:
    assert paths.requires_heavy(
        ["docs/benchmarks/v2.9.2-profiles.md", "python/dna_design_suite/api/__init__.py"]
    )


def test_planning_schemas_and_workflows_run_the_heavy_gate(paths) -> None:
    assert paths.requires_heavy(["planning/06-build/06-risk-register.md"])
    assert paths.requires_heavy(["schemas/v1/job_status.schema.json"])
    assert paths.requires_heavy([".github/workflows/ci.yml"])
    assert paths.requires_heavy(["tools/wheel_targets.py"])


def test_docs_parent_escape_cannot_hide_a_code_path(paths) -> None:
    assert paths.requires_heavy(["docs/../python/dna_design_suite/cli/main.py"])
    assert paths.requires_heavy(["docs/benchmarks/../../Cargo.toml"])


def test_empty_diff_runs_the_heavy_gate(paths) -> None:
    assert paths.requires_heavy([])


def test_cli_prints_a_github_output_line(paths) -> None:
    assert paths.github_output_line(["docs/project-status.md"]) == "heavy=false"
    assert paths.github_output_line(["Makefile"]) == "heavy=true"


def test_dispatch_default_is_linux_x86_64_only(wheels) -> None:
    selected = wheels.select_platforms("workflow_dispatch", "linux-x86_64")
    assert _pairs(selected) == {("ubuntu-latest", "x86_64")}
    assert selected[0]["manylinux"] == "auto"


def test_empty_or_unknown_dispatch_stays_on_linux(wheels) -> None:
    assert _pairs(wheels.select_platforms("workflow_dispatch", "")) == {("ubuntu-latest", "x86_64")}
    assert _pairs(wheels.select_platforms("workflow_dispatch", "macos")) == {
        ("ubuntu-latest", "x86_64")
    }


def test_dispatch_all_and_a_version_tag_build_every_platform(wheels) -> None:
    assert _pairs(wheels.select_platforms("workflow_dispatch", "all")) == PUBLISH_PLATFORMS
    assert _pairs(wheels.select_platforms("push", "")) == PUBLISH_PLATFORMS
    linux = [item for item in wheels.ALL_PLATFORMS if item["os"] == "ubuntu-latest"]
    assert {item["manylinux"] for item in linux} == {"auto"}


def test_cli_json_is_one_line() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "tools" / "wheel_targets.py"),
            "--event-name",
            "workflow_dispatch",
            "--targets",
            "linux-x86_64",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "\n" not in result.stdout.strip()
    assert '"target":"x86_64"' in result.stdout
    assert "macos" not in result.stdout


def _workflow(path: Path) -> dict:
    loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    # PyYAML 1.1 reads a bare `on:` key as boolean true.
    if True in loaded and "on" not in loaded:
        loaded["on"] = loaded.pop(True)
    return loaded


def test_wheels_do_not_run_on_branch_pushes_or_a_schedule() -> None:
    """A bare `push:` or a cron would rebuild every platform for a commit."""
    text = WHEELS.read_text(encoding="utf-8")
    assert "cron:" not in text
    document = _workflow(WHEELS)
    trigger = document["on"]
    assert "schedule" not in trigger
    assert "release" not in trigger
    assert "workflow_dispatch" in trigger
    push = trigger["push"]
    assert set(push) == {"tags"}
    assert "v*.*.*" in push["tags"]
    targets = trigger["workflow_dispatch"]["inputs"]["targets"]
    assert targets["default"] == "linux-x86_64"
    assert targets["options"] == ["linux-x86_64", "all"]
    assert "tools/wheel_targets.py" in text
    assert "sccache:" in text
    build = document["jobs"]["build"]
    assert build["needs"] == "legal-clearance"
    assert "needs.legal-clearance.result == 'success'" in build["if"]
    clearance = document["jobs"]["legal-clearance"]["if"]
    assert "github.ref_name" in clearance
    assert "contains(github.ref_name, '-')" in clearance


def test_ci_skips_heavy_jobs_for_docs_and_stays_on_ubuntu() -> None:
    document = _workflow(CI)
    assert document["concurrency"]["cancel-in-progress"] is True
    group = document["concurrency"]["group"]
    assert "github.head_ref" in group
    assert "github.ref_name" in group
    jobs = document["jobs"]
    always = {"changes", "invariants"}
    assert set(jobs) >= always | {
        "fmt-and-lint",
        "rust-crate",
        "pyext-builds",
        "python",
        "generated-assets",
    }
    for name, job in jobs.items():
        assert job["runs-on"] == "ubuntu-latest", name
        if name in always:
            assert "heavy" not in str(job.get("if", ""))
            continue
        assert job["if"] == "needs.changes.outputs.heavy == 'true'"
        needs = job["needs"]
        assert needs == "changes" or (isinstance(needs, list) and "changes" in needs)
    invariants = jobs["invariants"]
    assert "needs" not in invariants
    text = CI.read_text(encoding="utf-8")
    assert "tools/ci_path_filter.py" in text
    assert "macos-latest" not in text
    assert "windows-latest" not in text
    assert "Swatinem/rust-cache@v2" in text
    assert "cache: pip" in text
    fmt = jobs["fmt-and-lint"]
    uses = [step.get("uses", "") for step in fmt["steps"]]
    assert any(item.startswith("Swatinem/rust-cache@") for item in uses)
