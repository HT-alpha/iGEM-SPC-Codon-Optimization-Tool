"""The DNA Design Suite page and `python -m dna_design_suite.api`.

The entry point binds port 8000, opens a browser, and blocks. The page is a
guided form over the same routes. Download text matches the CLI writers.
"""

from __future__ import annotations

import base64
import errno
import json
import os
import shutil
import subprocess
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

import pytest

from dna_design_suite.api import run_page, serve
from dna_design_suite.api.server import _EXAMPLE_FILES
from dna_design_suite.cli import writers
from dna_design_suite.modelreg import ArtifactStore

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "python" / "dna_design_suite" / "api" / "static"
FORMATS = STATIC / "formats.js"
INDEX = STATIC / "index.html"
APP = STATIC / "app.js"
NODE = shutil.which("node")


class _Controls(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.controls: list[tuple[str, dict[str, str], bool]] = []
        self.label_for: set[str] = set()
        self._label = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr = {key: value or "" for key, value in attrs}
        if tag == "label":
            self._label += 1
            if attr.get("for"):
                self.label_for.add(attr["for"])
        if tag in {"input", "select", "textarea"}:
            self.controls.append((tag, attr, self._label > 0))

    def handle_endtag(self, tag: str) -> None:
        if tag == "label" and self._label:
            self._label -= 1


def _front() -> tuple[dict, dict]:
    front = {
        "job_id": "job_demo",
        "front_id": "front_demo",
        "objective_order": ["usage.cai", "usage.gc_content"],
        "objective_sense": ["maximize", "minimize"],
        "candidates": [
            {
                "candidate_id": "cand_abc",
                "provenance": {"origin": "tier1_seed"},
                "constraint_report": {"feasible": True},
                "nondomination_rank": 0,
                "objective_vector": [0.5, None],
                "sequences": [{"cds": "ATG" + "GCT" * 25}],
            },
            {
                "candidate_id": "cand_def",
                "provenance": {"origin": "tier2_block_ca"},
                "constraint_report": {"feasible": False},
                "nondomination_rank": 1,
                "objective_vector": [0.25, 0.5],
                "sequences": [{"cds": "ATGTTT"}],
            },
        ],
    }
    document = {
        "cds": front["candidates"][0]["sequences"][0]["cds"],
        "five_prime_flank": "AA",
        "three_prime_flank": "TT",
        "translation_table": 11,
    }
    return front, document


def _need_node() -> None:
    if NODE:
        return
    if os.environ.get("CF_REQUIRE_NODE") == "1":
        pytest.fail("node is required (CF_REQUIRE_NODE=1) and was not found")
    pytest.skip("node is not installed")


def test_download_text_matches_the_cli_writers() -> None:
    """A byte drift here is a sequence a biologist would file under the wrong name."""
    _need_node()
    front, document = _front()
    payload = base64.b64encode(json.dumps({"front": front, "document": document}).encode()).decode()
    module = base64.b64encode(FORMATS.read_bytes()).decode()
    script = (
        'import("data:text/javascript;base64,' + module + '").then((m) => {'
        f'const input = JSON.parse(Buffer.from("{payload}", "base64").toString());'
        "const out = {"
        "fasta: m.toFasta(input.front),"
        "csv: m.toCsv(input.front),"
        "genbank: m.toGenbank(input.front, input.document, 0),"
        "};"
        "process.stdout.write(JSON.stringify(out));"
        "}).catch((e) => { console.error(e); process.exit(1); });"
    )
    done = subprocess.run([NODE, "-e", script], capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr
    rendered = json.loads(done.stdout)
    assert rendered["fasta"] == writers.to_fasta(front)
    assert rendered["csv"] == writers.to_csv(front)
    assert rendered["genbank"] == writers.to_genbank(front, document, 0)


def test_every_form_control_has_a_label() -> None:
    parser = _Controls()
    parser.feed(INDEX.read_text(encoding="utf-8"))
    missing = []
    for tag, attr, inside in parser.controls:
        named = inside or attr.get("id") in parser.label_for or attr.get("aria-label")
        if not named:
            missing.append(attr.get("id") or tag)
    assert not missing, missing


def test_the_page_explains_the_words_a_first_run_needs() -> None:
    page = INDEX.read_text(encoding="utf-8")
    script = APP.read_text(encoding="utf-8")
    for phrase in (
        "Codon adaptation index",
        "GC window",
        "mg1655_ribosomal_usage.json",
        "rpsL",
        "What this version does not do",
        "name=\"viewport\"",
        "Run design",
        "Find primer pairs",
        "Select the whole sequence",
    ):
        assert phrase in page, phrase
    assert "What to do." in script
    assert "CF-CLI-003" in script
    assert "This version does not provide that feature." in script
    for code in ("CF-AUTH-001", "CF-ORCH-006", "CF-SRCH-010", "CF-PRIM-001", "CF-FEAS-005"):
        assert code in script


def test_the_module_entry_point_asks_for_port_8000_and_a_browser(monkeypatch) -> None:
    captured: dict = {}

    def fake_run_page(
        store,
        host="127.0.0.1",
        port=8000,
        *,
        open_browser=True,
        opener=None,
        block=True,
    ):
        captured.update(
            host=host, port=port, open_browser=open_browser, block=block, store=store
        )
        if open_browser:
            (opener or (lambda url: captured.setdefault("opened", []).append(url)))(
                f"http://{host}:{port}/"
            )

    monkeypatch.setattr("dna_design_suite.api.__main__.run_page", fake_run_page)
    from dna_design_suite.api.__main__ import main

    assert main([]) == 0
    assert captured["host"] == "127.0.0.1"
    assert captured["port"] == 8000
    assert captured["open_browser"] is True
    assert captured["block"] is True
    assert captured["opened"] == ["http://127.0.0.1:8000/"]

    captured.clear()
    assert main(["--no-browser", "--port", "8765"]) == 0
    assert captured["open_browser"] is False
    assert captured["port"] == 8765
    assert "opened" not in captured


def test_run_page_defaults_to_port_8000_and_opens_that_address(tmp_path: Path) -> None:
    opened: list[str] = []
    try:
        server = run_page(
            ArtifactStore(tmp_path / "artifacts"),
            open_browser=True,
            opener=opened.append,
            block=False,
        )
    except OSError as exc:
        assert exc.errno == errno.EADDRINUSE
        return
    try:
        assert server.server_address[1] == 8000
        assert opened == ["http://127.0.0.1:8000/"]
        with urllib.request.urlopen(opened[0]) as response:
            assert response.status == 200
            assert "Run design" in response.read().decode()
    finally:
        server.shutdown()
        server.server_close()


@pytest.fixture()
def base(tmp_path: Path):
    server = serve(ArtifactStore(tmp_path / "artifacts"))
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()
    server.server_close()


def test_the_page_names_the_same_illegal_base_as_the_api() -> None:
    """cleanedTemplate is the page copy of clean_template.

    The page keeps the illegal base on screen so the position is visible.
    The API refuses it. Both number that base from the header-stripped sequence.
    """
    if NODE is None:
        pytest.skip("node is not installed")
    source = APP.read_text(encoding="utf-8")
    start = source.index("function cleanedTemplate")
    end = source.index("function cleanedCds")
    script = source[start:end] + (
        "const samples = ["
        "  '>GFP demo\\natgc',"
        "  '>GFP demo\\nATGNCCC',"
        "  '  atg c  ',"
        "];"
        "const tail = 'ATG' + 'GCT'.repeat(237) + 'TAATAA';"
        "console.log(JSON.stringify({"
        "  cleaned: samples.map((item) => cleanedTemplate(item)),"
        "  trimmed: firstStopTrim(tail).length,"
        "  untrimmed: tail.length,"
        "}));"
    )
    completed = subprocess.run(
        [NODE, "-e", script], check=False, capture_output=True, text=True
    )
    assert completed.returncode == 0, completed.stderr
    reported = json.loads(completed.stdout)
    cleaned = reported["cleaned"]
    assert cleaned[0] == {"sequence": "ATGC", "illegal": []}
    assert cleaned[1]["sequence"] == "ATGNCCC"
    assert cleaned[1]["illegal"] == [{"base": "N", "position": 4}]
    assert cleaned[2] == {"sequence": "ATGC", "illegal": []}
    assert reported["untrimmed"] == 720
    assert reported["trimmed"] == 717
    page = INDEX.read_text(encoding="utf-8")
    assert 'id="trim-stop"' in page
    assert "data-preset=\"biobrick_rfc10.json\"" in page
    assert "data-preset=\"golden_gate.json\"" in page
    assert "data-preset=\"gibson.json\"" in page
    assert "data-preset=\"synthesis.json\"" in page
    assert "checks nothing" in page
    assert "informational, not a filter" in APP.read_text(encoding="utf-8")
    wired = APP.read_text(encoding="utf-8")
    for route in (
        "/v1/hosts",
        "/v1/models",
        "/v1/feasibility",
        "/v1/jobs",
        "objectives_preset",
        "structure_engines",
    ):
        assert route in wired, route
    assert 'kind: "protein"' in wired
    assert 'kind: "construct"' in wired
    assert 'kind: "single_cds"' in wired
    assert "BBa_E0040" in page


def test_product_sequence_is_sliced_from_the_submitted_template() -> None:
    """The page copies the amplicon without waiting for `product_sequence`."""
    _need_node()
    source = APP.read_text(encoding="utf-8")
    start = source.index("function productSequenceOf")
    end = source.index("function detailPanel")
    script = source[start:end] + (
        "const template = 'ACGTACGTAC';"
        "const cases = ["
        "  productSequenceOf({forward_start: 2, reverse_end: 6}, template),"
        "  productSequenceOf({product_sequence: 'ZZ', forward_start: 0, reverse_end: 1}, template),"
        "  productSequenceOf({forward_start: 0}, template),"
        "  productSequenceOf({forward_start: 1.5, reverse_end: 4}, template),"
        "  productSequenceOf({forward_start: 8, reverse_end: 12}, template),"
        "];"
        "console.log(JSON.stringify(cases));"
    )
    completed = subprocess.run(
        [NODE, "-e", script], check=False, capture_output=True, text=True
    )
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout) == ["GTAC", "ZZ", "", "", ""]
    assert "submittedPrimerTemplate = sequence" in source


def test_the_packaged_examples_match_the_repository_copies() -> None:
    packaged = ROOT / "python" / "dna_design_suite" / "api" / "examples"
    for _path, (name, _kind) in _EXAMPLE_FILES.items():
        assert (packaged / name).read_bytes() == (ROOT / "examples" / name).read_bytes()


def test_the_bundled_examples_are_the_files_on_disk(base: str) -> None:
    for path, (name, _kind) in _EXAMPLE_FILES.items():
        with urllib.request.urlopen(base + path) as response:
            assert response.status == 200
            body = response.read()
        assert body == (ROOT / "examples" / name).read_bytes()


def test_formats_module_is_served(base: str) -> None:
    with urllib.request.urlopen(f"{base}/formats.js") as response:
        assert response.status == 200
        assert response.headers.get("Content-Type", "").startswith("text/javascript")
        body = response.read().decode()
    assert "export function toFasta" in body
    assert "export function toGenbank" in body
