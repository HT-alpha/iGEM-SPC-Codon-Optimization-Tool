"""TST-BEN-001 and TST-BEN-002, plus the external-tool argument check.

The reference point is a function of the published ranges only. Two fronts
with different extents and the same ranges get the same point; that is the
property a nadir-of-this-front implementation would fail.
"""

from __future__ import annotations

import stat
import sys
from pathlib import Path

import pytest

from codonforge.bench import (
    BenchError,
    ExternalTool,
    RunRecord,
    Trial,
    assert_comparable,
    cai_greedy,
    check_reproduces,
    decide,
    fixed_reference_point,
    native,
    render_report,
    run_external,
    write_result,
)

RANGES = [(0.0, 1.0), (0.0, 10.0)]
SENSES = ["minimize", "maximize"]


def test_the_reference_point_is_the_published_nadir() -> None:
    assert fixed_reference_point(RANGES, SENSES) == [1.0, 0.0]
    again = fixed_reference_point(RANGES, SENSES)
    assert again == fixed_reference_point(RANGES, SENSES)


def test_a_reversed_range_is_refused() -> None:
    with pytest.raises(BenchError, match="reversed"):
        fixed_reference_point([(1.0, 0.0)], ["minimize"])


def test_comparison_refuses_a_different_reference_point() -> None:
    shared = dict(
        protein="p",
        seed=1,
        budget_tier="standard",
        lockfile_digest="sha256:" + "ab" * 32,
        hypervolume=0.2,
    )
    left = RunRecord(arm="codonforge", reference_point=(1.0, 0.0), **shared)
    right = RunRecord(arm="native", reference_point=(1.0, 0.0), **shared)
    assert_comparable(left, right)
    drifted = RunRecord(arm="native", reference_point=(2.0, 0.0), **shared)
    with pytest.raises(BenchError, match="reference"):
        assert_comparable(left, drifted)


def test_comparison_refuses_a_different_lockfile_and_tier() -> None:
    base = RunRecord(
        arm="codonforge",
        protein="p",
        seed=1,
        budget_tier="standard",
        lockfile_digest="sha256:" + "ab" * 32,
        reference_point=(1.0, 0.0),
        hypervolume=0.2,
    )
    other_tier = RunRecord(
        arm="native",
        protein="p",
        seed=1,
        budget_tier="fast",
        lockfile_digest=base.lockfile_digest,
        reference_point=base.reference_point,
        hypervolume=0.2,
    )
    with pytest.raises(BenchError, match="budget tier"):
        assert_comparable(base, other_tier)
    other_lock = RunRecord(
        arm="native",
        protein="p",
        seed=1,
        budget_tier="standard",
        lockfile_digest="sha256:" + "cd" * 32,
        reference_point=base.reference_point,
        hypervolume=0.2,
    )
    with pytest.raises(BenchError, match="lockfiles"):
        assert_comparable(base, other_lock)


def test_a_stored_result_recomputes_to_the_stored_hypervolume(tmp_path: Path) -> None:
    digest = "sha256:" + "11" * 32
    reference = fixed_reference_point(RANGES, SENSES)
    points = [[0.2, 8.0], [0.4, 6.0]]
    manifest = {
        "suite": "unit",
        "arm": "native",
        "protein": "fixture",
        "seed": 0,
        "budget_tier": "standard",
        "lockfile_digest": digest,
        "reference_point": reference,
    }
    front = {"landscape_lockfile_digest": digest, "candidates": []}
    metric = write_result(
        tmp_path,
        manifest=manifest,
        front=front,
        lockfile={"digest": digest},
        senses=SENSES,
        points=points,
    )
    assert check_reproduces(tmp_path, SENSES, points) == metric["value"]
    page = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert str(metric["value"]) in page
    assert "hypervolume" in page


def test_cai_greedy_picks_the_highest_weight_synonym() -> None:
    # ATG is Met (one codon). CTG/TTA/TTG/CTA/CTC/CTT are Leu.
    # The table prefers TTA, so the Leu codon changes and ATG does not.
    usage = {"ATG": 1.0, "TTA": 0.9, "TTG": 0.1, "CTT": 0.1, "CTC": 0.1, "CTA": 0.1, "CTG": 0.1}
    assert native("ATGCTG") == "ATGCTG"
    assert cai_greedy("ATGCTG", 1, usage) == "ATGTTA"


def test_an_external_tool_must_receive_the_input_and_keep_the_protein(tmp_path: Path) -> None:
    script = tmp_path / "tool.py"
    script.write_text(
        "import sys\n"
        "from pathlib import Path\n"
        "text = Path(sys.argv[1]).read_text()\n"
        "seq = ''.join(line.strip() for line in text.splitlines() if not line.startswith('>'))\n"
        "print('>out')\n"
        "print(seq)\n",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    tool = ExternalTool("identity", sys.executable, (str(script), "{input}"))
    produced = run_external(tool, "ATGCTGCGTTAA", 1, tmp_path / "work")
    assert produced == "ATGCTGCGTTAA"
    bare = ExternalTool("identity", sys.executable, (str(script),))
    with pytest.raises(BenchError, match="input"):
        bare.command(tmp_path / "input.fa")


def test_exp05_keeps_lamarckian_unless_both_clauses_hold() -> None:
    # Twenty seeds. Baldwinian hypervolume is higher by more than the
    # Lamarckian spread, but repair invocations are not under 2×.
    lamarckian = [Trial(hypervolume=1.0, diversity=0.4, repair_invocations=10) for _ in range(20)]
    costly = [Trial(hypervolume=3.0, diversity=0.5, repair_invocations=30) for _ in range(20)]
    decision = decide(lamarckian, costly)
    assert decision.switch_to_baldwinian is False
    assert decision.repair_ratio >= 2.0
    cheap = [Trial(hypervolume=3.0, diversity=0.5, repair_invocations=15) for _ in range(20)]
    switched = decide(lamarckian, cheap)
    assert switched.switch_to_baldwinian is True
    with pytest.raises(BenchError, match="20"):
        decide(lamarckian[:5], cheap[:5])


def test_the_report_copies_the_metric_value() -> None:
    page = render_report(
        [
            {
                "suite": "unit",
                "arm": "native",
                "protein": "p",
                "seed": 1,
                "metric": "hypervolume",
                "value": 0.125,
            }
        ]
    )
    assert "0.125" in page
