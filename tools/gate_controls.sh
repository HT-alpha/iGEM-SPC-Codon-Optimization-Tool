#!/bin/sh
# Negative control for the CAI-only gate claim.
#
# The committed claim is 50 of 50. Setting aceB's recorded CodonForge score
# back to 0, below cai_greedy, must be refused by check_cai_gate. Reverts
# only the file it edits.
set -eu
cd "$(dirname "$0")/.."
fail() { echo "CONTROL FAILED: $1" >&2; exit 1; }
cp benchmarks/corpus/cai_only.json /tmp/cai_only.json.gate-control
trap 'status=$?; cp /tmp/cai_only.json.gate-control benchmarks/corpus/cai_only.json; exit "$status"' EXIT

python3 - <<'PY'
import sys
sys.path.insert(0, "tools")
import patch

patch.sub(
    "benchmarks/corpus/cai_only.json",
    '''    "gene": "aceB",
    "protein_id": "NP_418438.1",
    "codons": 534,
    "native_cai": 0.0,
    "cai_greedy_cai": 1.0,
    "codonforge_cai": 1.0,''',
    '''    "gene": "aceB",
    "protein_id": "NP_418438.1",
    "codons": 534,
    "native_cai": 0.0,
    "cai_greedy_cai": 1.0,
    "codonforge_cai": 0.0,''',
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
    if "underperform" not in text or "aceB" not in text:
        raise SystemExit(f"refused for the wrong reason: {text}")
    print(f"CAUGHT {text}")
    raise SystemExit(0)
raise SystemExit("a score of 0 on aceB was silent")
PY
code=$?
set -e
if [ "$code" -ne 0 ]; then
  cat /tmp/gate-control.log >&2
  fail "below-greedy aceB score was not refused (exit $code)"
fi
echo "CAUGHT below-greedy aceB score"
grep underperform /tmp/gate-control.log
