#!/usr/bin/env bash
# Negative controls for round 20 (recommended-axis disclosure).
#
# P1-26 was re-read and not activated, so there is no overlapping-block
# fallback to inject against. The controls here are for Part 2.
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

# 1. THE GAP: an out-of-core axis on the front is silently dropped from
#    the recommendation. Omitting not_informed_by, or emitting it empty
#    while the structure term is on the front, is the failure TST-CLI-012
#    exists to catch.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/design.py'; s=open(p).read()
old='"not_informed_by": list(scalarization["not_informed_by"]),'
new='"not_informed_by": [],'
assert old in s, "anchor 1"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/1: recommendation silently drops an out-of-core axis" \
  "$PY python/tests/test_cli.py -q -k cli_012"

# 2. A second, independent axes list that could drift from the scalarization.
#    Building informed_by from plan.scored minus a hardcoded L0 name, rather
#    than from the score loop, is the duplicated-rule pattern.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/design.py'; s=open(p).read()
old='"informed_by": list(scalarization["informed_by"]),'
new='"informed_by": [t for t in job.plan.scored if t != "usage.cai"],'
assert old in s, "anchor 2"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/2: informed_by is rebuilt from plan.scored, not the score loop" \
  "$PY python/tests/test_cli.py -q -k 'cli_012 or cli_009'"

# 3. The writer re-derives exclusion instead of rendering the document.
python3 - <<'PYEOF'
p='python/codonforge/cli/writers.py'; s=open(p).read()
old='        excluded = ", ".join(recommended.get("not_informed_by") or []) or "(none)"'
new='        excluded = ", ".join(recommended.get("informed_by") or []) or "(none)"'
assert old in s, "anchor 3"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/3: explain swaps not_informed_by for informed_by" \
  "$PY python/tests/test_writers.py -q -k recommendation_axes"

control_summary
exit 0
