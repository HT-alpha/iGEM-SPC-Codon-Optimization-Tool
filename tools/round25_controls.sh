#!/usr/bin/env bash
# Negative controls for front reduction, comparison refusal, idempotency,
# and the external-tool argv guard. Each injection is the branch the check
# actually takes.
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

# 1. EXTREME UNION. Presentation starts from the protected extremes. Clearing
#    that set in the reference-direction path lets the cap of two drop an axis
#    extreme. Clearing a different method, or a comment, leaves this test green.
python3 - <<'PYEOF' || exit 1
p = "python/codonforge/decision/reduce.py"
s = open(p).read()
old = """    chosen = set(protected)
    if not usable:
        return chosen
    rows = _normalized([_row(c, senses) for c in usable])
    objectives = len(senses)
"""
new = """    chosen = set()
    if not usable:
        return chosen
    rows = _normalized([_row(c, senses) for c in usable])
    objectives = len(senses)
"""
assert old in s, "anchor 1"
assert s.count(old) == 1, "direction union is not unique"
open(p, "w").write(s.replace(old, new, 1))
PYEOF
control "R/1: reference-direction reduction drops an axis extreme" \
  "$PY python/tests/test_front_reduction.py -q -k reduction_keeps_every_axis_extreme"

# 2. REFERENCE POINT. The third comparison check is the reference point.
#    Deleting it lets two hypervolumes with different nadirs compare.
python3 - <<'PYEOF' || exit 1
p = "python/codonforge/bench/compare.py"
s = open(p).read()
old = """    if left.reference_point != right.reference_point:
        raise BenchError(
            "the two hypervolumes use different reference points and cannot be compared"
        )
"""
assert old in s, "anchor 2"
assert s.count(old) == 1, "reference check is not unique"
open(p, "w").write(s.replace(old, "", 1))
PYEOF
control "R/2: comparison accepts a different reference point" \
  "$PY python/tests/test_bench_harness.py -q -k different_reference_point"

# 3. IDEMPOTENCY HASH. Same key and a different body must be 409. Skipping
#    the hash compare replays the original job.
python3 - <<'PYEOF' || exit 1
p = "python/codonforge/api/server.py"
s = open(p).read()
old = """                    if previous_hash != digest:
                        return 409, problem(
"""
new = """                    if False and previous_hash != digest:
                        return 409, problem(
"""
assert old in s, "anchor 3"
assert s.count(old) == 1, "idempotency hash check is not unique"
open(p, "w").write(s.replace(old, new, 1))
PYEOF
control "R/3: idempotency ignores a different body" \
  "$PY python/tests/test_api.py -q -k idempotency_replays"

# 4. ARGV FIDELITY. command() refuses a template with no {input} before the
#    process starts. Deleting that guard makes bare.command return.
python3 - <<'PYEOF' || exit 1
p = "python/codonforge/bench/baselines.py"
s = open(p).read()
old = """        if not any("{input}" in arg for arg in self.argv):
            raise BenchError(
                f"{self.name}: the argv template does not contain {{input}}, so the "
                "sequence would never be passed to the tool"
            )
"""
assert old in s, "anchor 4"
assert s.count(old) == 1, "argv guard is not unique"
open(p, "w").write(s.replace(old, "", 1))
PYEOF
control "R/4: an external tool runs without receiving the input path" \
  "$PY python/tests/test_bench_harness.py -q -k external_tool_must_receive_the_input"

control_summary
