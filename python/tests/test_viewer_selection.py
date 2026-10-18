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
    assert 'import { nucleotideIndex } from "/viewer.js";' in js


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
