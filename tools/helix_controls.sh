#!/usr/bin/env bash
# Negative controls for the DNA helix (python/tests/test_helix.py).
#
# Each injection breaks one property the tests claim to guard, runs the test
# that names it, and reverts. A control that does not fire means the test does
# not check what its name says. Run from a clean tree.
#
# CF_REQUIRE_NODE=1 matters: without it the Node-based tests skip, a skip exits
# zero, and every control would read SILENT for a reason that has nothing to do
# with the tests.
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

export CF_REQUIRE_NODE=1
PY=".venv/bin/python -m pytest -q -p no:cacheprovider python/tests/test_helix.py"
HELIX="python/codonforge/api/static/helix.js"
CSS="python/codonforge/api/static/codonforge-ui.css"
APP="python/codonforge/api/static/app.js"
SERVER="python/codonforge/api/server.py"
TOOL="tools/measure_bform_1bna.py"

# Anchor-verified edit: fails loudly if the anchor is absent or not unique.
inject() {
  python3 - "$@" <<'PYEOF' || exit 1
import sys
sys.path.insert(0, "tools")
import patch
path, old, new = sys.argv[1:4]
patch.sub(path, old, new)
PYEOF
}

# --- geometry ---------------------------------------------------------------
inject "$HELIX" 'G: "C", C: "G" });' 'G: "C", C: "C" });'
control "H/1: C pairs with C, so strand B is not the complement" \
  "$PY -k geometry_is_b_form"

inject "$HELIX" '  twistPerBp: 36,' '  twistPerBp: 35,'
control "H/2: twist per base pair is 35, not 36" \
  "$PY -k 'geometry_is_b_form or constants_in_the_module'"

inject "$HELIX" '  phosphateRadius: 9.3,' '  phosphateRadius: 10,'
control "H/3: the radius no longer matches the recorded 1BNA measurement" \
  "$PY -k constants_in_the_module"

inject "$HELIX" '  strandOffset: 145,' '  strandOffset: 180,'
control "H/4: strand offset 180 degrees, so the grooves are equal" \
  "$PY -k constants_in_the_module"

inject "$HELIX" '"ATGACCATGATTACGGATTC"' '"ATGACCATGATTACGGATTG"'
control "H/5: one base of the rung sequence is not lacZ" \
  "$PY -k lacz"

inject "$HELIX" 'segments.sort((p, q) => p.front - q.front);' 'segments.sort((p, q) => q.front - p.front);'
control "H/6: backbone drawn front to back" \
  "$PY -k backbone_follows_the_helix"

# --- colours ----------------------------------------------------------------
inject "$CSS" '--base-g: #bd966b;' '--base-g: #ffff00;'
control "H/7: G is a saturated yellow, not the muted IGV colour" \
  "$PY -k base_colours"

inject "$CSS" '  stroke: var(--recommended);
  stroke-linecap: butt;
}

.hx-rung path {' '  stroke: #1a4a74;
  stroke-linecap: butt;
}

.hx-rung path {'
control "H/8: the backbone is a second blue, one digit off the tab token" \
  "$PY -k existing_colour_tokens"

inject "$CSS" '.helix-stage[data-helix-state="error"] .hx-rung path {
  stroke: var(--refusal-ink);' '.helix-stage[data-helix-state="error"] .hx-rung path {
  stroke: #ff0000;'
control "H/9: the error red is not the refused-zone chip red" \
  "$PY -k existing_colour_tokens"

# --- state machine ----------------------------------------------------------
inject "$HELIX" '    if (!turning().length) {' '    if (true) {'
control "H/10: the result pose is applied without waiting for a turn" \
  "$PY -k state_machine_waits"

inject "$HELIX" '        setRate(1);
        container.dataset.helixPose = "badge";' '        container.dataset.helixPose = "badge";'
control "H/11: the speed is not restored when the helix docks" \
  "$PY -k state_machine_waits"

inject "$HELIX" '        container.dataset.helixPose = "badge";
      });
    }
  }' '        container.dataset.helixPose = "badge";
      });
    } else {
      container.dataset.helixPose = "centre";
    }
  }'
control "H/12: an error snaps the helix back to the centre" \
  "$PY -k state_machine_waits"

inject "$HELIX" '    cancelPending();
    state = next;' '    state = next;'
control "H/13: a new state leaves the old dock pending" \
  "$PY -k new_state_cancels"

inject "$HELIX" 'if (event.target === backbone && event.animationName === "hx-slide") finish();' 'if (event.target === backbone) finish();'
control "H/14: any animation iteration, even a pulse, counts as the helix's turn" \
  "$PY -k state_machine_waits"

# --- independence and reduced motion ---------------------------------------
inject "$HELIX" '  let state = "idle";' '  let state = window.matchMedia ? "idle" : "idle";'
control "H/15: the module consults matchMedia" \
  "$PY -k touches_only_its_container"

inject "$HELIX" '  const { svg, backbone } = buildSvg(container.ownerDocument, helixGeometry(options));' '  const { svg, backbone } = buildSvg(document, helixGeometry(options));'
control "H/16: the module builds from the global document, not its container" \
  "$PY -k touches_only_its_container"

inject "$CSS" '  .hx-bb,
  .hx-rung .hx-a,
  .hx-rung .hx-b,
  .helix-stage[data-helix-state="running"] .hx-rung {
    animation: none;
  }' '  .hx-bb,
  .hx-rung .hx-a,
  .hx-rung .hx-b,
  .helix-stage[data-helix-state="running"] .hx-rung {
    animation-iteration-count: infinite;
  }'
control "H/17: the reduced-motion rule no longer switches the animations off" \
  "$PY -k reduced_motion_is_a_stylesheet_rule"

python3 - <<'PYEOF' || exit 1
import sys
sys.path.insert(0, "tools")
import patch
patch.append("python/codonforge/api/static/codonforge-ui.css", '\n.helix-stage[data-helix-pose="badge"] { height: 1px; }\n')
PYEOF
control "H/18: the docked pose applies outside the motion-allowed block" \
  "$PY -k reduced_motion_is_a_stylesheet_rule"

# --- wiring and serving -----------------------------------------------------
inject "$SERVER" '    "/helix.js": ("helix.js", "text/javascript; charset=utf-8"),
' ''
control "H/19: the server does not list /helix.js" \
  "$PY -k server_serves_the_module"

inject "$APP" '  frontRoot.innerHTML = lateReveal(body);' '  frontRoot.innerHTML = `<div class="reveal">${body}</div>`;'
control "H/20: a front failure is shown without the delayed reveal" \
  "$PY -k error_message_waits"

inject "$APP" '    codonHelix.setState(renderFront(front.payload) ? "result" : "error");' '    renderFront(front.payload);'
control "H/21: a rendered front no longer sets the result state" \
  "$PY -k polling_code_is_the_only_caller"

python3 - <<'PYEOF' || exit 1
import sys
sys.path.insert(0, "tools")
import patch
patch.append("python/codonforge/api/static/helix.js", "/* " + "x" * 20000 + " */\n")
PYEOF
control "H/22: 20 kB of padding in helix.js" \
  "$PY -k thirty_kilobytes"

# --- the measuring script ---------------------------------------------------
inject "$TOOL" '    if sign < 0:
        offset = (-offset) % 360.0' '    if False:
        offset = (-offset) % 360.0'
control "H/23: the measuring script does not mirror the offset of a left-handed helix" \
  "$PY -k hand_of_the_helix"

inject "$TOOL" '    best = nelder_mead(objective, [0.0, 0.0, 0.0, 0.0], [0.5, 0.5, 0.05, 0.05])' '    best = [0.0, 0.0, 0.0, 0.0]'
control "H/24: the axis is not refined, so the offset carries the end-effect bias" \
  "$PY -k measuring_script_recovers"

control_summary
