#!/usr/bin/env bash
# Negative controls for the stash/sentinel invariant and the extract() arm.
# Documentation round: the checks already exist. These injections hit the
# branches a future edit would actually take (R/4: the nearby branch is silent).
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

# 1. THE TAKE. search_tier3_segment returns the sentinel CfError as published
#    by to_py_err. The noise test must then record CF-LAND-012 instead of
#    raising. Deleting a comment, or forging on the JSON-parse branch, does
#    not hit this.
python3 - <<'PYEOF'
p='crates/cf-pyext/src/lib.rs'; s=open(p).read()
old='''    if let Some(escaped) = pipeline::take_escaped_scorer_error() {
        return Err(escaped);
    }
    result'''
new='    result'
assert old in s, "anchor 1"
open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/1: search_tier3_segment publishes the sentinel as CF-LAND-012" \
  "$PY python/tests/test_orchestrator.py -q -k non_problem_scorer_error"

# 2. THE EXTRACT ARM. A shape mismatch is classified like a raised scorer,
#    so a returned object whose `code` is CF-STRU-003 stops being CF-LAND-012.
#    The short-batch fixture extracts successfully and never reaches this line.
python3 - <<'PYEOF'
p='crates/cf-pyext/src/pipeline.rs'; s=open(p).read()
old='            .map_err(|e| CfError::new(code::LAND_BATCH_UNEXPECTED, e.to_string()))?;'
new='            .map_err(scorer_raised)?;'
assert old in s, "anchor 2"
assert s.count(old) == 1, "extract arm is not unique"
open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/2: extract() failures are classified by a code attribute" \
  "$PY python/tests/test_orchestrator.py -q -k ill_typed_l0_return"

control_summary
