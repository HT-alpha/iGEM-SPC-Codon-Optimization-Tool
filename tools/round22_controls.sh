#!/usr/bin/env bash
# Negative controls for DEP-SRCH-014 / DEP-ORCH-003
# (B≥R is CF-SRCH-008 alone; CF-SRCH-009 is config; mapper lifts every CF-*).
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

# 1. THE REMAINING B≥R GATE: Tier2Plan::check no longer refuses R > B.
#    The 009 copy in tier2_block_ca is gone. If this is silent, either the
#    production-path test is not hitting Rust, or a second gate still exists.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier2.rs'; s=open(p).read()
old='                    if codons > config.block {'
new='                    if false && codons > config.block {'
assert old in s, "anchor 1"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/1: Tier2Plan::check no longer refuses a term wider than the block" \
  "$PY python/tests/test_orchestrator.py -q -k wider_than_the_block"

# 2. THE DISTINCT 009 PATH: check no longer refuses a zero block.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier2.rs'; s=open(p).read()
old='        if config.block == 0 {'
new='        if false && config.block == 0 {'
assert old in s, "anchor 2"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/2: Tier2Plan::check no longer refuses a zero block" \
  "$PY python/tests/test_orchestrator.py -q -k zero_block"

# 3. MAPPER: only CF-SRCH-008 is lifted. A zero-block 009 must then
#    escape as a raw ValueError, which is the bug this round closes.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/job.py'; s=open(p).read()
old='    if not isinstance(code, str) or not code.startswith("CF-"):'
new='    if not isinstance(code, str) or code != "CF-SRCH-008":'
assert old in s, "anchor 3"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/3: mapper lifts only CF-SRCH-008, so a zero block is unstructured" \
  "$PY python/tests/test_orchestrator.py -q -k zero_block"

# 4. TITLE: hardcoded 008 title on every TierPlanError. A 009 must not
#    be documented as an unbounded term.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/budget.py'; s=open(p).read()
old='        self.title = title or "Tier 2 refused the job"'
new='        self.title = "A term couples positions more widely than the tier 2 block size"'
assert old in s, "anchor 4"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/4: TierPlanError hardcodes the CF-SRCH-008 title" \
  "$PY python/tests/test_orchestrator.py -q -k zero_block"

# 5. Planning validator: a cited DEP that is not in the index.
python3 - <<'PYEOF'
p='planning/03-algorithms/02-block-local-coordinate-ascent.md'; s=open(p).read()
old='`DEP-SRCH-014` — that is the only `B ≥ R` check.'
new='`DEP-SRCH-099` — that is the only `B ≥ R` check.'
assert old in s, "anchor 5"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/5: planning validator misses an unindexed DEP citation" \
  ".venv/bin/python planning/tools/validate_planning_package.py"

# 6. Release-clearance remains BLOCKED.
control "R/6: the release-clearance gate still refuses" \
  ".venv/bin/python tools/check_release_clearance.py"

control_summary
