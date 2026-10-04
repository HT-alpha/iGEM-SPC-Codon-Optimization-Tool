#!/usr/bin/env bash
# Negative controls for the Potts artifact digest and for EXP-05 decide().
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

# 1. DIGEST. Changing a coupling must change the digest. Disabling the
#    comparison lets a mutated document pass verify_artifact.
python3 - <<'PYEOF' || exit 1
p = "python/codonforge/epistasis/artifact.py"
s = open(p).read()
old = """    if recomputed != recorded:
"""
new = """    if False and recomputed != recorded:
"""
assert old in s, "anchor 1"
assert s.count(old) == 1, "digest comparison is not unique"
open(p, "w").write(s.replace(old, new, 1))
PYEOF
control "R/1: a mutated coupling still matches the recorded digest" \
  "$PY python/tests/test_epistasis.py -q -k artifact_digest_covers_a_coupling"

# 2. SEED COUNT. Fewer than 20 seeds is not a decision. Disabling the
#    length check lets a 5-seed call return a decision.
python3 - <<'PYEOF' || exit 1
p = "python/codonforge/bench/exp05.py"
s = open(p).read()
old = """    if len(lamarckian) < SEED_COUNT or len(baldwinian) < SEED_COUNT:
"""
new = """    if False and (len(lamarckian) < SEED_COUNT or len(baldwinian) < SEED_COUNT):
"""
assert old in s, "anchor 2"
assert s.count(old) == 1, "seed-count gate is not unique"
open(p, "w").write(s.replace(old, new, 1))
PYEOF
control "R/2: five seeds are accepted as an EXP-05 decision" \
  "$PY python/tests/test_bench_harness.py -q -k exp05_keeps_lamarckian"

# 3. THE RULE. A clear Baldwinian win switches. Reversing the margin test
#    keeps Lamarckian when the written rule would switch.
python3 - <<'PYEOF' || exit 1
p = "python/codonforge/bench/exp05.py"
s = open(p).read()
old = """    switch = margin > deviation and ratio < 2.0
"""
new = """    switch = margin < deviation and ratio < 2.0
"""
assert old in s, "anchor 3"
assert s.count(old) == 1, "decision inequality is not unique"
open(p, "w").write(s.replace(old, new, 1))
PYEOF
control "R/3: a clear hypervolume win does not switch" \
  "$PY python/tests/test_bench_harness.py -q -k exp05_keeps_lamarckian"

control_summary
