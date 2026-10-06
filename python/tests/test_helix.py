"""The DNA helix (python/codonforge/api/static/helix.js and its stylesheet block).

The helix is a frontend-only element, so most of what matters is checked from
three directions. The numbers it draws must be traceable: the radius and the
strand offset come from PDB 1BNA through tools/measure_bform_1bna.py, the
rungs carry a real sequence from the shipped corpus, and the base colours are
a stated transform of IGV's palette. The geometry must be B-form: 36 degrees
and 3.4 angstroms per base pair, complementary and antiparallel strands, a
right-handed turn. And the module must stay inside its container, with
prefers-reduced-motion a stylesheet rule and not a script check.

The geometry and the state machine are exercised by running helix.js under
Node, which every GitHub-hosted runner has. CF_REQUIRE_NODE=1 turns a missing
Node into a failure, so CI cannot pass by not running them. A browser is not a
test dependency; the rendered result was checked in Chrome for the release
notes, and the stylesheet is checked here as text.
"""

from __future__ import annotations

import base64
import colorsys
import gzip
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import urllib.request
from pathlib import Path

import pytest

from codonforge.api import serve
from codonforge.modelreg import ArtifactStore

ROOT = Path(__file__).resolve().parents[2]
STATIC = ROOT / "python" / "codonforge" / "api" / "static"
HELIX_JS = STATIC / "helix.js"
APP_JS = STATIC / "app.js"
CSS = STATIC / "codonforge-ui.css"
INDEX = STATIC / "index.html"
RECORD = ROOT / "docs" / "helix-1bna-measurement.json"
CORPUS = ROOT / "benchmarks" / "corpus"

NODE = shutil.which("node")


def _need_node() -> None:
    if NODE:
        return
    if os.environ.get("CF_REQUIRE_NODE") == "1":
        pytest.fail("node is required (CF_REQUIRE_NODE=1) and was not found")
    pytest.skip("node is not installed")


def _node(body: str) -> object:
    """Run `body` with the helix module bound to `m` and return its JSON output.

    The module is loaded from a data: URL, which Node always treats as an ES
    module, so the check does not depend on a package.json `type` field or on
    the Node version's module detection.
    """
    _need_node()
    source = base64.b64encode(HELIX_JS.read_bytes()).decode()
    script = (
        f'import("data:text/javascript;base64,{source}").then(async (m) => {{'
        f"const out = await (async () => {{ {body} }})();"
        "process.stdout.write(JSON.stringify(out));"
        "}).catch((e) => { console.error(e); process.exit(1); });"
    )
    done = subprocess.run([NODE, "-e", script], capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr
    return json.loads(done.stdout)


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"(^|[^:])//.*", r"\1", text)


def _load_tool():
    spec = importlib.util.spec_from_file_location(
        "measure_bform_1bna", ROOT / "tools" / "measure_bform_1bna.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _synthetic_duplex(
    radius: float = 9.0,
    twist: float = 36.0,
    rise: float = 3.4,
    offset: float = 150.0,
    hand: int = 1,
    pairs: int = 20,
) -> str:
    """A duplex of whole turns with known parameters, in PDB ATOM format.

    Chain A is residues 1..pairs and chain B the next `pairs`, residue i paired
    with 2 * pairs + 1 - i, as in 1BNA. Phosphates sit at `radius` on two
    helical tracks `offset` degrees apart and C1' atoms at 5.6 A. Whole turns
    keep the principal-axis estimate exact; 1.2 turns tilts it by a few
    hundredths of an angstrom. This is a fixture for the measuring script, not
    data the UI shows.
    """
    lines = []
    serial = 1

    def atom(name: str, chain: str, number: int, radius_: float, angle: float, z: float) -> None:
        nonlocal serial
        x = radius_ * math.cos(math.radians(angle))
        y = radius_ * math.sin(math.radians(angle))
        lines.append(
            f"ATOM  {serial:5d} {name:<4s} DN{chain} {chain}{number:4d}"
            f"    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00"
        )
        serial += 1

    for i in range(1, pairs + 1):
        z = rise * (i - 1)
        angle = hand * twist * (i - 1)
        partner = 2 * pairs + 1 - i
        if i > 1:
            atom("P", "A", i, radius, angle, z)
        if partner < 2 * pairs:
            atom("P", "B", partner, radius, angle + hand * offset, z)
        atom("C1'", "A", i, 5.6, angle, z)
        atom("C1'", "B", partner, 5.6, angle + hand * 130.0, z)
    return "\n".join(lines) + "\n"


def test_the_measuring_script_recovers_known_parameters() -> None:
    measure = _load_tool().measure
    got = measure(_synthetic_duplex())
    assert got["handedness"] == "right"
    assert got["phosphate_radius_angstrom"]["mean"] == pytest.approx(9.0, abs=0.05)
    assert got["twist_deg_per_bp"]["mean"] == pytest.approx(36.0, abs=0.1)
    assert got["rise_angstrom_per_bp"]["mean"] == pytest.approx(3.4, abs=0.05)
    assert got["phosphate_track_offset_deg"] == pytest.approx(150.0, abs=1.0)
    assert got["minor_groove_arc_deg"] == pytest.approx(150.0, abs=1.0)
    assert got["major_groove_arc_deg"] == pytest.approx(210.0, abs=1.0)


def test_the_measuring_script_tells_the_hand_of_the_helix() -> None:
    measure = _load_tool().measure
    left = measure(_synthetic_duplex(hand=-1))
    assert left["handedness"] == "left"
    assert left["phosphate_track_offset_deg"] == pytest.approx(150.0, abs=1.0)


def test_the_constants_in_the_module_come_from_the_recorded_measurement() -> None:
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    assert record["pdb_id"] == "1BNA"
    assert record["handedness"] == "right"
    assert len(record["source_sha256"]) == 64
    form = _node("return m.B_FORM;")

    # Measured on 1BNA, so equal to the record after rounding.
    assert form["phosphateRadius"] == round(record["phosphate_radius_angstrom"]["mean"], 1)
    assert form["strandOffset"] == round(record["phosphate_track_offset_deg"])

    # Published idealised values the brief asks for; 1BNA agrees within two
    # standard errors of its own step-to-step spread.
    twist, rise = record["twist_deg_per_bp"], record["rise_angstrom_per_bp"]
    assert form["twistPerBp"] == 36
    assert form["risePerBp"] == 3.4
    assert abs(form["twistPerBp"] - twist["mean"]) <= 2 * twist["sd"] / math.sqrt(twist["n"])
    assert abs(form["risePerBp"] - rise["mean"]) <= 2 * rise["sd"] / math.sqrt(rise["n"])

    # Ten base pairs and 34 angstroms per turn, and a minor groove narrower
    # than the major groove, so the two grooves alternate.
    assert 360 / form["twistPerBp"] == 10
    assert form["risePerBp"] * 10 == pytest.approx(34.0)
    assert 90 < form["strandOffset"] < 180


def test_the_rungs_carry_the_start_of_lacz_from_the_shipped_corpus() -> None:
    panel = json.loads((CORPUS / "panel.json").read_text(encoding="utf-8"))
    raw = (CORPUS / "NC_000913.3_cds.fna.gz").read_bytes()
    assert hashlib.sha256(raw).hexdigest() == panel["cds_gzip_sha256"]
    records = gzip.decompress(raw).decode().split(">")[1:]
    (lacz,) = [r for r in records if "[gene=lacZ]" in r.split("\n", 1)[0]]
    sequence = "".join(lacz.split("\n")[1:])
    default = _node("return m.DEFAULT_SEQUENCE;")
    assert len(default) == 20
    assert sequence.startswith(default)


def _hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))  # type: ignore[return-value]


# IGV.js, js/util/nucleotideColors.js, defaultNucleotideColors (igvteam/igv.js).
IGV_DEFAULT = {"A": (0, 200, 0), "C": (0, 0, 200), "T": (255, 0, 0), "G": (209, 113, 5)}


def _muted(rgb: tuple[int, int, int]) -> tuple[int, int, int]:
    """HSL saturation x 0.4, lightness 0.58, hue kept."""
    h, _l, s = colorsys.rgb_to_hls(*(v / 255 for v in rgb))
    r, g, b = colorsys.hls_to_rgb(h, 0.58, s * 0.4)
    return round(r * 255), round(g * 255), round(b * 255)


def test_base_colours_are_the_muted_igv_defaults() -> None:
    css = CSS.read_text(encoding="utf-8")
    for base, source in IGV_DEFAULT.items():
        found = re.search(rf"--base-{base.lower()}:\s*(#[0-9a-f]{{6}});", css)
        assert found, f"--base-{base.lower()} is not defined"
        got = _hex_to_rgb(found.group(1))
        want = _muted(source)
        assert all(abs(a - b) <= 1 for a, b in zip(got, want, strict=True)), (base, got, want)
    assert "igvteam/igv.js" in css
    assert "G rgb(209,113,5)" in css


def test_the_helix_uses_the_existing_colour_tokens() -> None:
    css = CSS.read_text(encoding="utf-8")
    # The backbone is the tab indicator's token, not a second blue.
    assert re.search(r"--recommended:\s*#1a4a73;", css)
    assert re.search(r"button\.tab\.active\s*\{[^}]*var\(--recommended\)", css)
    assert re.search(r"\.hx-bb path\s*\{[^}]*stroke:\s*var\(--recommended\)", css)
    # The error colour is the one the refused-zone chip uses.
    assert re.search(r"\.zone-red\s*\{[^}]*color:\s*var\(--refusal-ink\)", css)
    assert re.search(
        r'data-helix-state="error"\] \.hx-rung path\s*\{[^}]*var\(--refusal-ink\)', css
    )
    assert re.search(r"--refusal-ink:\s*#5c1515;", css)


def test_the_geometry_is_b_form() -> None:
    g = _node("return m.helixGeometry();")
    form = _node("return m.B_FORM;")
    rungs = g["rungs"]
    assert len(rungs) == 20
    assert g["sequence"] == "".join(r["baseA"] for r in rungs)

    # Rise and radius keep the published ratio 3.4 : 9.2, whatever the scale.
    spacing = [rungs[i]["y"] - rungs[i + 1]["y"] for i in range(len(rungs) - 1)]
    assert all(s == pytest.approx(spacing[0], abs=0.15) for s in spacing)
    assert g["rise"] / g["radius"] == pytest.approx(form["risePerBp"] / form["phosphateRadius"])
    assert g["pitch"] == pytest.approx(g["rise"] * 10)
    assert g["bpPerTurn"] == 10

    # 36 degrees per base pair, so a rung ten places up faces the same way.
    step = [rungs[i + 1]["angleA"] - rungs[i]["angleA"] for i in range(len(rungs) - 1)]
    assert step == [form["twistPerBp"]] * (len(rungs) - 1)

    # Strand B is the complement at every rung, so read 5' to 3' from the top
    # it is the reverse complement of strand A: the strands are antiparallel.
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    assert [r["baseB"] for r in rungs] == [comp[r["baseA"]] for r in rungs]
    strand_b = "".join(r["baseB"] for r in reversed(rungs))
    assert strand_b == "".join(comp[b] for b in reversed(g["sequence"]))


def test_the_backbone_follows_the_helix_and_is_drawn_back_to_front() -> None:
    g = _node("return m.helixGeometry();")
    segments = g["segments"]
    fronts = [s["front"] for s in segments]
    assert fronts == sorted(fronts), "painter's order must hold at every angle"
    assert {s["strand"] for s in segments} == {"a", "b"}

    def bezier(d: str, t: float) -> tuple[float, float, float, float]:
        nums = [float(n) for n in re.findall(r"-?\d+(?:\.\d+)?", d)]
        (x0, y0, x1, y1, x2, y2, x3, y3) = nums
        u = 1 - t
        x = u**3 * x0 + 3 * u * u * t * x1 + 3 * u * t * t * x2 + t**3 * x3
        y = u**3 * y0 + 3 * u * u * t * y1 + 3 * u * t * t * y2 + t**3 * y3
        return x, y, x0 - x3, y0 - y3

    offset = {"a": 0.0, "b": g["strandOffset"]}
    twist = _node("return m.B_FORM.twistPerBp;")
    worst = 0.0
    for s in segments:
        x, y, dx, dy = bezier(s["d"], 0.5)
        j = (g["yBottom"] - y) / g["rise"]
        angle = math.radians(g["startAngle"] + offset[s["strand"]] + twist * j)
        exact = g["centre"] + g["radius"] * math.cos(angle)
        worst = max(worst, abs(x - exact))
        # Right-handed: wherever a strand is clearly at the front it climbs to
        # the right (x grows as y falls).
        if s["front"] >= 0.8:
            assert dx * dy < 0, s
    assert worst < 0.35, f"a backbone segment strays {worst:.2f} units from the helix"


def test_sequences_are_validated() -> None:
    result = _node(
        'try { m.helixGeometry({ sequence: "ACGX" }); return "no error"; }'
        " catch (e) { return e.name; }"
    )
    assert result == "RangeError"
    ok = _node('return m.helixGeometry({ sequence: "acgtacgtac" }).rungs.length;')
    assert ok == 10


STUB = r"""
const log = [];
const made = [];
function node(name) {
  const style = { props: {}, setProperty(k, v) { this.props[k] = v; } };
  const n = { name, attrs: {}, children: [], listeners: {}, style };
  n.setAttribute = (k, v) => { n.attrs[k] = v; };
  n.append = (...c) => { n.children.push(...c); };
  n.replaceChildren = (...c) => { n.children = c; };
  n.dataset = {};
  n.addEventListener = (t, f) => { (n.listeners[t] = n.listeners[t] || []).push(f); };
  n.removeEventListener = (t, f) => {
    n.listeners[t] = (n.listeners[t] || []).filter((x) => x !== f);
  };
  n.getAnimations = () => (n.isSvg ? animations : []);
  made.push(n);
  return n;
}
const animations = [];
const doc = {
  createElementNS: (ns, name) => {
    const n = node(name);
    n.isSvg = name === "svg";
    return n;
  },
};
const container = node("div");
container.ownerDocument = doc;
"""


def test_the_state_machine_waits_for_a_whole_turn_before_docking() -> None:
    got = _node(
        STUB
        + r"""
for (const kind of ["hx-move", "hx-fade", "hx-slide"]) {
  const a = { animationName: kind, rates: [], updatePlaybackRate(r) { this.rates.push(r); } };
  animations.push(a);
}
animations.push({
  animationName: undefined,
  updatePlaybackRate() { throw new Error("touched a transition"); },
});
const helix = m.createHelix(container);
const snap = () => ({ state: container.dataset.helixState, pose: container.dataset.helixPose });
const backbone = made.find((n) => n.name === "g" && n.attrs.class === "hx-bb");
const fire = (event) => [...(backbone.listeners.animationiteration || [])].forEach((f) => f(event));
const steps = {};
steps.start = snap();
helix.setState("running");
steps.running = snap();
steps.runningRates = animations.slice(0, 3).map((a) => a.rates.at(-1));
helix.setState("result");
steps.resultBeforeTurn = snap();
steps.ratesBeforeTurn = animations.slice(0, 3).map((a) => a.rates.at(-1));
// A pulse's iteration, or another target, must not count as the helix's turn.
fire({ target: {}, animationName: "hx-slide" });
fire({ target: backbone, animationName: "hx-glow" });
steps.afterWrongEvents = snap();
fire({ target: backbone, animationName: "hx-slide" });
steps.afterTurn = snap();
steps.ratesAfterTurn = animations.slice(0, 3).map((a) => a.rates.at(-1));
steps.listenersLeft = (backbone.listeners.animationiteration || []).length;
helix.setState("running");
helix.setState("error");
steps.error = snap();
let bad = "none";
try { helix.setState("done"); } catch (e) { bad = e.name; }
steps.bad = bad;
return steps;
"""
    )
    assert got["start"] == {"state": "idle", "pose": "centre"}
    assert got["running"] == {"state": "running", "pose": "corner"}
    assert got["runningRates"] == [1.8, 1.8, 1.8]  # 9 s idle turn -> 5 s running turn
    assert got["resultBeforeTurn"] == {"state": "result", "pose": "corner"}
    assert got["ratesBeforeTurn"] == [1.8, 1.8, 1.8], "the speed holds until the turn ends"
    assert got["afterWrongEvents"]["pose"] == "corner"
    assert got["afterTurn"] == {"state": "result", "pose": "badge"}
    assert got["ratesAfterTurn"] == [1, 1, 1]
    assert got["listenersLeft"] == 0
    assert got["error"] == {"state": "error", "pose": "corner"}, "an error stays where it is"
    assert got["bad"] == "RangeError"


def test_without_a_running_animation_the_result_pose_applies_at_once() -> None:
    got = _node(
        STUB
        + r"""
const helix = m.createHelix(container);
helix.setState("running");
helix.setState("result");
return { state: container.dataset.helixState, pose: container.dataset.helixPose };
"""
    )
    assert got == {"state": "result", "pose": "badge"}


def test_a_new_state_cancels_a_pending_dock() -> None:
    got = _node(
        STUB
        + r"""
animations.push({ animationName: "hx-slide", updatePlaybackRate() {} });
const helix = m.createHelix(container);
const backbone = made.find((n) => n.name === "g" && n.attrs.class === "hx-bb");
helix.setState("running");
helix.setState("result");
helix.setState("running");
const handlersAfterRerun = (backbone.listeners.animationiteration || []).length;
helix.setState("error");
const { helixState: state, helixPose: pose } = container.dataset;
return { handlersAfterRerun, state, pose };
"""
    )
    assert got == {"handlersAfterRerun": 0, "state": "error", "pose": "corner"}


def test_the_module_touches_only_its_container_and_no_script_media_check() -> None:
    code = _strip_comments(HELIX_JS.read_text(encoding="utf-8"))
    for forbidden in (
        r"\bdocument\b",
        r"\bwindow\b",
        r"\bglobalThis\b",
        r"querySelector",
        r"getElementById",
        r"requestAnimationFrame",
        r"\bmatchMedia\b",
        r"setInterval",
        r"visibilityState",
        r"\blocalStorage\b",
        r"\bfetch\b",
    ):
        assert not re.search(forbidden, code), f"helix.js must not use {forbidden}"
    assert "ownerDocument" in code
    assert "export function createHelix" in code
    assert "export function helixGeometry" in code


def test_the_helix_file_and_the_svg_it_builds_fit_in_thirty_kilobytes() -> None:
    size = HELIX_JS.stat().st_size
    markup = _node(
        STUB
        + r"""
const esc = (s) => String(s).replace(/&/g, "&amp;").replace(/"/g, "&quot;");
function html(n) {
  const attrs = Object.entries(n.attrs).map(([k, v]) => ` ${k}="${esc(v)}"`).join("");
  const css = Object.entries(n.style.props).map(([k, v]) => `${k}: ${v};`).join(" ");
  const extra = [
    n.style.transform && `transform: ${n.style.transform};`,
    n.style.opacity && `opacity: ${n.style.opacity};`,
  ].filter(Boolean).join(" ");
  const style = [css, extra].filter(Boolean).join(" ");
  const open = `<${n.name}${attrs}${style ? ` style="${esc(style)}"` : ""}`;
  return n.children.length ? `${open}>${n.children.map(html).join("")}</${n.name}>` : `${open}/>`;
}
m.createHelix(container);
return html(container.children[0]).length;
"""
    )
    assert size < 30_000
    assert size + markup < 30_000, f"JS {size} B plus SVG {markup} B"


def _block(css: str, opener: str) -> str:
    start = css.index(opener)
    depth = 0
    for i in range(css.index("{", start), len(css)):
        depth += {"{": 1, "}": -1}.get(css[i], 0)
        if depth == 0:
            return css[start : i + 1]
    raise AssertionError(f"unbalanced block after {opener}")


def test_reduced_motion_is_a_stylesheet_rule() -> None:
    css = CSS.read_text(encoding="utf-8")
    reduce = _block(css, "@media (prefers-reduced-motion: reduce) {\n  .hx-bb")
    for selector in (".hx-bb", ".hx-rung .hx-a", ".hx-rung .hx-b", '[data-helix-state="running"]'):
        assert selector in reduce
    assert "animation: none" in reduce
    assert "transition: none" in reduce

    # The moves and resizes between poses exist only when motion is allowed.
    allowed = _block(css, "@media (prefers-reduced-motion: no-preference) {\n  .helix-stage")
    assert 'data-helix-pose="corner"' in allowed
    assert 'data-helix-pose="badge"' in allowed
    base = css.replace(allowed, "")
    assert 'data-helix-pose="badge"' not in base

    # Every helix animation is declared with a plain keyframe name the reduce
    # block can switch off, and none of them is started from script.
    code = _strip_comments(HELIX_JS.read_text(encoding="utf-8"))
    assert ".animate(" not in code
    assert "animation =" not in code


def test_the_error_message_waits_for_the_helix() -> None:
    css = CSS.read_text(encoding="utf-8")
    assert re.search(r"\.reveal\.reveal-late\s*\{[^}]*animation-delay:\s*1500ms", css)
    js = APP_JS.read_text(encoding="utf-8")
    # No error callout is still wrapped in the plain, undelayed reveal.
    assert not re.search(r'class="reveal">\$\{[^}]*(problemCallout|callout refusal|body)', js)
    assert js.count("lateReveal(") >= 12


def test_the_polling_code_is_the_only_caller_and_covers_all_four_states() -> None:
    js = APP_JS.read_text(encoding="utf-8")
    assert 'import { createHelix } from "/helix.js";' in js
    assert js.count("createHelix(") == 2, "one instance per tab"
    for tab in ("codonHelix", "primerHelix"):
        states = set(re.findall(rf'{tab}\.setState\("(\w+)"\)', js))
        assert {"running", "error"} <= states, tab
        assert len(re.findall(rf"{tab}\.setState\(", js)) >= 4, tab
    assert 'codonHelix.setState(renderFront(front.payload) ? "result" : "error")' in js
    primer_result = "renderPrimerResult(status.payload.primer_result)"
    assert f'primerHelix.setState({primer_result} ? "result" : "error")' in js
    # Nothing else in the page reaches into the module.
    assert len(re.findall(r"\bsetState\(", js)) == len(
        re.findall(r"(?:codon|primer)Helix\.setState\(", js)
    )
    assert "helixGeometry" not in js


@pytest.fixture()
def base(tmp_path: Path):
    server = serve(ArtifactStore(tmp_path / "artifacts"))
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()


def test_the_server_serves_the_module_and_the_page_loads_it_in_both_tabs(base: str) -> None:
    with urllib.request.urlopen(f"{base}/helix.js") as response:
        assert response.status == 200
        assert response.headers.get("Content-Type", "").startswith("text/javascript")
        body = response.read().decode()
    assert "export function createHelix" in body
    with urllib.request.urlopen(f"{base}/") as response:
        page = response.read().decode()
    assert '<script type="module" src="/app.js"></script>' in page
    for stage in ("codon-helix", "primer-helix"):
        assert f'id="{stage}" class="helix-stage" aria-hidden="true"' in page
    # The stage is the first thing in each result region and displaces none of
    # the existing content: it sits before the placeholders, which stay.
    assert page.index('id="codon-helix"') < page.index('<h2 id="status">')
    assert page.index('id="primer-helix"') < page.index('id="primer-error"')
    assert "No job yet." in page
    assert "No front yet." in page
    assert "NCBI BLAST is not part of this tool." in page
