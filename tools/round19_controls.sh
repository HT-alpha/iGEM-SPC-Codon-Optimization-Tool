#!/usr/bin/env bash
# Negative controls for round 19 (Parts A, B, C).
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

# Part B: the harness refuses to measure a stale build.
control_require_fresh_build

PY=".venv/bin/python -m pytest"

# ---- Part A: the three schema fields ----

# 1. `recommended` recomputed at the boundary instead of read from best().
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/design.py'; s=open(p).read()
old='"candidate_id": f"cand_{best.digest.split(\':\')[1][:16]}",'
new='"candidate_id": front_module.candidate_id_of(ordered[-1]) if False else f"cand_{ordered[-1].digest.split(\':\')[1][:16]}",'
assert old in s, "anchor 1"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "A/1: recommended names a candidate other than best()" \
  "$PY python/tests/test_cli.py -q -k cli_009"

# 2. role reports the ceiling instead of the effective registration.
python3 - <<'PYEOF'
p='crates/cf-pyext/src/pipeline.rs'; s=open(p).read()
old="out.insert(term.clone(), json!(provider.role.as_str()));"
new="out.insert(term.clone(), json!(provider.provider.max_role().as_str()));"
assert old in s, "anchor 2"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "A/2: role reports max_role rather than the effective role" \
  "$PY python/tests/test_cli.py -q -k 'cli_009 or p2_11'"

# 3. unknown_terms omitted when empty, so absent and complete look alike.
python3 - <<'PYEOF'
p='python/codonforge/decision/front.py'; s=open(p).read()
old='    entry["unknown_terms"] = list(unknown_terms or [])'
new='    if unknown_terms:\n        entry["unknown_terms"] = list(unknown_terms)'
assert old in s, "anchor 3"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "A/3: unknown_terms disappears when empty" \
  "$PY python/tests/test_cli.py -q -k cli_009"

# 4. feasibility made tri-state, re-conflating what DEP-CONS-004 separated.
python3 - <<'PYEOF'
p='python/codonforge/decision/front.py'; s=open(p).read()
old='            "feasible": feasible,'
new='            "feasible": feasible if not unknown_terms else "unknown",'
assert old in s, "anchor 4"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "A/4: feasible becomes tri-state" \
  "$PY python/tests/test_cli.py -q -k cli_009"

# ---- Part C: P2-11 ----

# 5. THE FINDING: a term that declares it cannot be ranked on is ranked anyway.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/budget.py'; s=open(p).read()
old='            if role != "pareto_axis":'
new='            if False:'
assert old in s, "anchor 5"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "C/5: the not-rankable guard is removed" \
  "$PY python/tests/test_orchestrator.py -q"

# 6. the structure term is scored in core by a cheap proxy.
python3 - <<'PYEOF'
p='crates/cf-landscape/src/terms/structure.rs'; s=open(p).read()
old='    fn score(&self, _ctx: &EvalContext<\'_>) -> Result<f64> {\n        Err(CfError::new('
new='    fn score(&self, _ctx: &EvalContext<\'_>) -> Result<f64> {\n        return Ok(0.5);\n        #[allow(unreachable_code)]\n        Err(CfError::new('
assert old in s, "anchor 6"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "C/6: an in-core proxy stands in for a real fold" \
  "cargo test -p cf-landscape 2>&1 && $PY python/tests/test_cli.py -q -k p2_11"

# 7. capabilities asserted rather than read from the registry.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/design.py'; s=open(p).read()
old='        capabilities=json.loads(_core.term_capabilities(_directions_request(request))),'
new='        capabilities={term: "L1" for term in request.objectives},'
assert old in s, "anchor 7"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "C/7: capabilities are asserted, not read" \
  "$PY python/tests/test_cli.py -q -k p2_11"

# 8. the L0 term routed into tier 1, which has no delta table for it.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/budget.py'; s=open(p).read()
old='                plan.scored.append(term)\n                plan.directions[term] = directions.get(term, "minimize")'
new='                plan.scored.append(term)\n                plan.tier1.append(term)\n                plan.directions[term] = directions.get(term, "minimize")'
assert old in s, "anchor 8"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "C/8: an out-of-core term is handed to tier 1" \
  "$PY python/tests/test_cli.py -q -k p2_11"

# 9. the diagnostic caveat survives the promotion to axis.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/design.py'; s=open(p).read()
old='    if structure_diag.TERM_ID in job.plan.scored:'
new='    if False:'
assert old in s, "anchor 9"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "C/9: CF-LAND-013 still fires for a term that is now an axis" \
  "$PY python/tests/test_cli.py -q -k p2_11"

# 10. the two term whitelists drift apart again.
python3 - <<'PYEOF'
p='python/codonforge/modelreg/resolve.py'; s=open(p).read()
old='    "structure.local_accessibility",\n)'
new=')'
assert old in s, "anchor 10"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "C/10: the CLI advertises a term the resolver rejects" \
  "$PY python/tests/test_cli.py -q -k cli_010"

# ---- validators, unchanged obligations ----

control "V/11: the release-clearance gate still refuses" \
  ".venv/bin/python tools/check_release_clearance.py"

python3 - <<'PYEOF'
p='python/codonforge/predictors/cache.py'; s=open(p).read()
old="class AccessibilityCache:"
new="# Intersected CFG with bifurcation rules over the codon lattice.\nclass AccessibilityCache:"
assert old in s, "anchor 12"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "V/12: the AMD-013 grammar guard still bites" \
  "$PY python/tests/test_repo_invariants.py -q -k amd_013"

control_summary
exit 0
