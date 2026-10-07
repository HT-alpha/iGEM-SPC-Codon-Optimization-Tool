#!/bin/sh
# Negative control for the CAI-only gate claim.
#
# A gates.json that claims all 50 proteins, with aceB and thrA no longer
# excluded, must be refused by check_cai_gate: those two recorded scores are
# below cai_greedy. Reverts only the file it edits.
set -eu
cd "$(dirname "$0")/.."
fail() { echo "CONTROL FAILED: $1" >&2; exit 1; }
trap 'status=$?; git checkout -- benchmarks/corpus/gates.json; exit "$status"' EXIT

python3 - <<'PY'
import sys
sys.path.insert(0, "tools")
import patch

patch.sub(
    "benchmarks/corpus/gates.json",
    '"n_claimed": 48,',
    '"n_claimed": 50,',
)
patch.sub(
    "benchmarks/corpus/gates.json",
    '"claim": "CAI-only front, 48 of 50 E. coli K-12 proteins (aceB and thrA excluded, see degenerate_cases), cai_greedy baseline, CodonForge meets or exceeds baseline on all 48."',
    '"claim": "CAI-only front, 50 of 50 E. coli K-12 proteins, cai_greedy baseline, CodonForge meets or exceeds baseline on all 50."',
)
patch.sub(
    "benchmarks/corpus/gates.json",
    '"excluded": [\n        "aceB",\n        "thrA"\n      ]',
    '"excluded": []',
)
PY

set +e
.venv/bin/python - <<'PY' >/tmp/gate-control.log 2>&1
import json
from codonforge.bench.corpus import BenchError, check_cai_gate

gate = json.loads(open("benchmarks/corpus/gates.json", encoding="utf-8").read())[0]
rows = json.loads(open("benchmarks/corpus/cai_only.json", encoding="utf-8").read())
try:
    check_cai_gate(gate, rows)
except BenchError as exc:
    text = str(exc)
    if "underperform" not in text:
        raise SystemExit(f"refused for the wrong reason: {text}")
    print(f"CAUGHT {text}")
    raise SystemExit(0)
raise SystemExit("all-50 claim was silent")
PY
code=$?
set -e
if [ "$code" -ne 0 ]; then
  cat /tmp/gate-control.log >&2
  fail "all-50 claim was not refused (exit $code)"
fi
echo "CAUGHT all-50 claim"
grep underperform /tmp/gate-control.log
