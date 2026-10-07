#!/bin/sh
# Negative control for the PrimerForge drag selection.
#
# Ignoring the element under the pointer is the Chrome bug: setPointerCapture
# retargets pointermove to the viewer, and a selection that only reads
# event.target stays on the nucleotide where the drag started. The retargeted
# test must fail. Reverts only the file it edits.
set -eu
cd "$(dirname "$0")/.."
fail() { echo "CONTROL FAILED: $1" >&2; exit 1; }

python3 - <<'PY'
import sys
sys.path.insert(0, "tools")
import patch

patch.sub(
    "python/codonforge/api/static/viewer.js",
    """export function nucleotideIndex(target, under) {
  const direct = closestIndex(target);
  if (direct !== null) return direct;
  return closestIndex(under);
}""",
    """export function nucleotideIndex(target, under) {
  void under;
  return closestIndex(target);
}""",
)
PY

set +e
CF_REQUIRE_NODE=1 .venv/bin/pytest -q -p no:cacheprovider \
  python/tests/test_viewer_selection.py -k retargeted \
  >/tmp/viewer-control.log 2>&1
code=$?
set -e
git checkout -- python/codonforge/api/static/viewer.js
if [ "$code" -eq 0 ]; then
  fail "retarget injection was silent"
fi
echo "CAUGHT retargeted drag (exit $code)"
grep -E 'AssertionError|FAILED' /tmp/viewer-control.log | head -4
