"""Drag selection in the Primer design sequence viewer.

Chrome retargets pointermove to the element that called setPointerCapture.
The viewer does that, so event.target during a drag is the viewer and not the
nucleotide under the pointer. nucleotideIndex must still return that
nucleotide, from the element at the pointer's coordinates.

The behaviour is checked by running viewer.js under Node. CF_REQUIRE_NODE=1
turns a missing Node into a failure. A browser check of the same drag is in
the v2.3.1 notes; this test is what CI runs.
"""

from __future__ import annotations

import base64
import json
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path

import pytest

from dna_design_suite.api import serve
from dna_design_suite.modelreg import ArtifactStore

ROOT = Path(__file__).resolve().parents[2]
VIEWER_JS = ROOT / "python" / "dna_design_suite" / "api" / "static" / "viewer.js"
APP_JS = ROOT / "python" / "dna_design_suite" / "api" / "static" / "app.js"
NODE = shutil.which("node")

STUB = r"""
function node(index, parent) {
  return {
    dataset: index === undefined ? {} : { index: String(index) },
    closest(selector) {
      if (selector !== "[data-index]") return null;
      if (this.dataset.index !== undefined && this.dataset.index !== "") return this;
      return parent ? parent.closest(selector) : null;
    },
  };
}
const viewer = node();
const at = (index) => node(index, viewer);
"""


def _need_node() -> None:
    if NODE:
        return
    if os.environ.get("CF_REQUIRE_NODE") == "1":
        pytest.fail("node is required (CF_REQUIRE_NODE=1) and was not found")
    pytest.skip("node is not installed")


def _node(body: str) -> object:
    _need_node()
    source = base64.b64encode(VIEWER_JS.read_bytes()).decode()
    script = (
        'import("data:text/javascript;base64,' + source + '").then(async (m) => {'
        "const out = await (async () => { " + STUB + "\n" + body + " })();"
        "process.stdout.write(JSON.stringify(out));"
        "}).catch((e) => { console.error(e); process.exit(1); });"
    )
    done = subprocess.run([NODE, "-e", script], capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr
    return json.loads(done.stdout)


def test_a_direct_hit_names_that_nucleotide() -> None:
    assert _node("return m.nucleotideIndex(at(300), viewer);") == 300


def test_a_retargeted_move_names_the_nucleotide_under_the_pointer() -> None:
    """The capture target is the viewer. The nucleotide is `under`."""
    assert _node("return m.nucleotideIndex(viewer, at(345));") == 345


def test_a_move_over_empty_viewer_names_nothing() -> None:
    assert _node("return m.nucleotideIndex(viewer, viewer);") is None


def test_a_non_integer_index_is_not_a_nucleotide() -> None:
    assert _node('const bad = node("x", viewer); return m.nucleotideIndex(viewer, bad);') is None


def test_the_page_asks_for_the_element_under_the_pointer_and_keeps_capture() -> None:
    js = APP_JS.read_text(encoding="utf-8")
    assert (
        "nucleotideIndex(event.target, document.elementFromPoint(event.clientX, event.clientY))"
        in js
    )
    assert "viewer.setPointerCapture(event.pointerId)" in js
    assert 'from "/viewer.js"' in js
    assert "nucleotideIndex" in js
    assert "nearestLineIndex" in js
    assert "nearestBase" in js
    assert "renderSequenceLines" in js
    assert "selectionLabel" in js


def test_a_sequence_longer_than_one_line_wraps() -> None:
    """Sixty bases fill one line. The next base starts the next line."""
    report = _node(
        """
        const width = m.BASES_PER_LINE;
        const html = m.renderSequenceLines("A".repeat(width + 1), 0);
        return {
          lines: (html.match(/class="seq-line"/g) || []).length,
          opensSecond: html.includes('data-line-start="' + width + '"'),
          numbersSecond: html.includes('>' + (width + 1) + '</span>'),
          keepsIndex: html.includes('data-index="' + width + '"'),
        };
        """
    )
    assert report == {
        "lines": 2,
        "opensSecond": True,
        "numbersSecond": True,
        "keepsIndex": True,
    }


def test_line_numbers_are_one_based_and_follow_the_window_origin() -> None:
    html = _node('return m.renderSequenceLines("ACGT", 100);')
    assert 'data-index="100"' in html
    assert 'data-index="103"' in html
    assert ">101</span>" in html
    assert ">104</span>" in html


def test_an_illegal_base_is_escaped() -> None:
    html = _node('return m.renderSequenceLines("A<G&T", 0);')
    assert "&lt;" in html
    assert "&amp;" in html
    assert "<G" not in html


def test_the_selection_label_is_one_based_and_inclusive() -> None:
    """The rpsL demo posts {start: 40, end: 300}, which is bases 41–300."""
    assert _node("return m.selectionLabel({start: 40, end: 300});") == {
        "startBase": 41,
        "endBase": 300,
        "length": 260,
    }
    assert _node("return m.selectionLabel(null);") is None
    assert _node("return m.selectionLabel({start: 5, end: 5});") is None


def test_a_pointer_between_lines_names_the_nearer_line() -> None:
    assert _node(
        """
        const lines = [{top: 0, bottom: 20}, {top: 28, bottom: 48}];
        return m.nearestLineIndex(26, lines, 8);
        """
    ) == 1
    assert _node(
        """
        const lines = [{top: 0, bottom: 20}, {top: 28, bottom: 48}];
        return m.nearestLineIndex(200, lines, 8);
        """
    ) is None


def test_a_pointer_in_a_gap_names_the_nearer_base() -> None:
    assert _node(
        """
        const bases = [
          {index: 61, left: 8, right: 16},
          {index: 70, left: 80, right: 88},
        ];
        return [m.nearestBase(40, bases), m.nearestBase(70, bases), m.nearestBase(12, bases)];
        """
    ) == [61, 70, 61]


@pytest.fixture()
def base(tmp_path: Path):
    server = serve(ArtifactStore(tmp_path / "artifacts"))
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()


def test_the_server_serves_the_viewer_module(base: str) -> None:
    with urllib.request.urlopen(f"{base}/viewer.js") as response:
        assert response.status == 200
        assert response.headers.get("Content-Type", "").startswith("text/javascript")
        body = response.read().decode()
    assert "export function nucleotideIndex" in body
