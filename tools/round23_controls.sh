#!/usr/bin/env bash
# Negative controls for DEP-ORCH-004
# (one Job.run lift for every CF-* problem; a non-document scorer error still crashes).
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

# 1. THE SINGLE CATCH. Job.run no longer lifts ValueError into problem.
#    A short L0 batch must then escape as a raw ValueError, which is the
#    gap this round closes.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/job.py'; s=open(p).read()
old='            refusal = _engine_refusal(exc, self._engine_stage)\n            if refusal is exc:\n                raise'
new='            raise exc'
assert old in s, "anchor 1"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/1: Job.run no longer lifts a tier-3 CF-* ValueError" \
  "$PY python/tests/test_orchestrator.py -q -k short_l0_batch"

# 2. THE BATCH BOUNDARY FORGES 012. Every scorer exception becomes
#    CF-LAND-012, so a programming error is absorbed as a problem.
python3 - <<'PYEOF'
p='crates/cf-pyext/src/pipeline.rs'; s=open(p).read()
old='        let returned = self.callable.call1((sequences,)).map_err(scorer_raised)?;'
new='        let returned = self.callable.call1((sequences,)).map_err(|e| CfError::new(code::LAND_BATCH_UNEXPECTED, e.to_string()))?;'
assert old in s, "anchor 2"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/2: every scorer exception is relabeled CF-LAND-012" \
  "$PY python/tests/test_orchestrator.py -q -k non_problem_scorer_error"

# 3. A coded engine error is also relabeled, so CF-STRU-003 never arrives.
python3 - <<'PYEOF'
p='crates/cf-pyext/src/pipeline.rs'; s=open(p).read()
old='        let returned = self.callable.call1((sequences,)).map_err(scorer_raised)?;'
new='        let returned = self.callable.call1((sequences,)).map_err(|e| CfError::new(code::LAND_BATCH_UNEXPECTED, e.to_string()))?;'
assert old in s, "anchor 3"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/3: a missing engine is reported as CF-LAND-012" \
  "$PY python/tests/test_orchestrator.py -q -k missing_structure_engine"

# 4. THE MAPPER FORGES A DOCUMENT. Noise becomes TierPlanError.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/job.py'; s=open(p).read()
old='    except json.JSONDecodeError:\n        return exc'
new='    except json.JSONDecodeError:\n        return TierPlanError(code="CF-LAND-012", detail=str(exc), subject="batch", stage=stage, title="forged", status=500)'
assert old in s, "anchor 4"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/4: the mapper forges a problem for a non-document ValueError" \
  "$PY python/tests/test_orchestrator.py -q -k non_problem_scorer_error"

# 5. Planning validator: a cited DEP that is not in the index.
python3 - <<'PYEOF'
p='planning/03-algorithms/02-block-local-coordinate-ascent.md'; s=open(p).read()
old='`DEP-ORCH-004` is that mapping'
new='`DEP-ORCH-099` is that mapping'
assert old in s, "anchor 5"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/5: planning validator misses an unindexed DEP citation" \
  ".venv/bin/python planning/tools/validate_planning_package.py"

# 6. Release-clearance remains BLOCKED.
control "R/6: the release-clearance gate still refuses" \
  ".venv/bin/python tools/check_release_clearance.py"

# 7. Schema sync: an unregistered code cited from the registry comment is
#    not the check. Cite a fake code from a schema-adjacent source the
#    checker reads: the Rust registry.
python3 - <<'PYEOF'
p='crates/cf-types/src/error.rs'; s=open(p).read()
old='    LAND_BATCH_UNEXPECTED  => ("CF-LAND-012", 500, "Batch results do not match the outstanding batch", "planning/03-algorithms/03-nsga-loop-and-operators.md"),'
new='    LAND_BATCH_UNEXPECTED  => ("CF-LAND-099", 500, "Batch results do not match the outstanding batch", "planning/03-algorithms/03-nsga-loop-and-operators.md"),'
assert old in s, "anchor 7"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/7: schema sync misses an unregistered error code" \
  ".venv/bin/python tools/check_schema_sync.py --check"

# 8. ruff sees an unknown import.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/job.py'; s=open(p).read()
old='from __future__ import annotations'
new='from __future__ import annotations\nimport this_module_does_not_exist'
assert old in s, "anchor 8"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/8: ruff misses an unknown import" \
  ".venv/bin/ruff check python/codonforge/orchestrator/job.py"

# 9. clippy 1.98 sees an unused binding. The toolchain file must be in
#    force; a 1.83 default would still catch `unused_variables` under -D warnings.
python3 - <<'PYEOF'
p='crates/cf-pyext/src/pipeline.rs'; s=open(p).read()
old='fn scorer_raised(err: pyo3::PyErr) -> CfError {'
new='fn scorer_raised(err: pyo3::PyErr) -> CfError {\n    let unused_binding = 1;'
assert old in s, "anchor 9"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/9: clippy misses an unused binding" \
  "cargo clippy -p cf-pyext --all-targets -- -D warnings"

control_summary
