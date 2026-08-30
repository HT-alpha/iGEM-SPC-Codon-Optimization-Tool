#!/usr/bin/env bash
# Negative controls for round 17 (P2-10, the diagnostic wiring, the sixth rule).
set -u
cd "$(dirname "$0")/.."
pass=0; fail=0
control() {
  local name="$1"; shift; local cmd="$1"; shift
  echo "=================================================================="
  echo "CONTROL: $name"
  if eval "$cmd" >/tmp/c17.out 2>&1; then
    echo "  SILENT — the check did not fire"; fail=$((fail+1))
  else
    echo "  CAUGHT:"; grep -E 'AssertionError|FAILED|Error|assert' /tmp/c17.out | head -3 | sed 's/^/    /'
    pass=$((pass+1))
  fi
  git checkout -- . 2>/dev/null
}
PYT=".venv/bin/python -m pytest"
CACHE="python/tests/test_structure_cache.py"

# ---- the three classes the sixth rule names, against the cache ----

# 1. STALE: the cache ignores the edit range and reuses everything.
#    Caught only by agreement — this is the class TST-STRU-002 does cover.
python3 - <<'PY'
p='python/codonforge/predictors/cache.py'; s=open(p).read()
old="        for index in range(low, high):"
new="        for index in range(low, low):"
assert old in s, "anchor 1"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-10/1 [stale]: the edit range is ignored and values are reused" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 $PYT $CACHE -q"

# 2. NOT-A-CACHE: correct on every input, recomputes everything, every time.
#    Passes TST-STRU-002 perfectly. This is the class the sixth rule adds.
python3 - <<'PY'
p='python/codonforge/predictors/cache.py'; s=open(p).read()
old="        edited = self._sequence[:start] + replacement + self._sequence[stop:]"
new="        edited = self._sequence[:start] + replacement + self._sequence[stop:]\n        self._sequence = edited\n        return self.prime(edited)"
assert old in s, "anchor 2"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-10/2 [not-a-cache]: correct, but refolds everything every time" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 $PYT $CACHE -q"

# 3. ARGUMENT FIDELITY: the affected range stops depending on where the edit was.
python3 - <<'PY'
p='python/codonforge/predictors/cache.py'; s=open(p).read()
old="        first_end = max(self.window_nt - 1, start - self.span_nt)\n        last_end = min(len(self._sequence) - 1, stop - 1 + self.span_nt)"
new="        first_end = self.window_nt - 1\n        last_end = len(self._sequence) - 1"
assert old in s, "anchor 3"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-10/3 [argument fidelity]: the edit position stops mattering" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 $PYT $CACHE -q"

# 4. PLAUSIBILITY (range): a value outside [0,1] is served.
python3 - <<'PY'
p='python/codonforge/predictors/cache.py'; s=open(p).read()
old="                values[index] = folded[offset]"
new="                values[index] = folded[offset] * 1.5"
assert old in s, "anchor 4"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-10/4 [plausibility]: accessibility escapes [0, 1]" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 $PYT $CACHE -q"

# 5. PLAUSIBILITY (scaling): recomputation starts growing with sequence length.
python3 - <<'PY'
p='python/codonforge/predictors/cache.py'; s=open(p).read()
old="        first_end = max(self.window_nt - 1, start - self.span_nt)"
new="        first_end = max(self.window_nt - 1, start - self.span_nt - len(self._sequence) // 4)"
assert old in s, "anchor 5"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-10/5 [plausibility]: work grows with the construct" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 $PYT $CACHE -q"

# 6. A one-window offset in the reuse loop — the subtle index bug.
python3 - <<'PY'
p='python/codonforge/predictors/cache.py'; s=open(p).read()
old="            offset = end_in_sequence - region_start - (self.window_nt - 1)"
new="            offset = end_in_sequence - region_start - self.window_nt"
assert old in s, "anchor 6"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-10/6: a one-window offset in the reuse loop" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 $PYT $CACHE -q"

# 7. The unprimed refusal becomes a plausible empty answer.
python3 - <<'PY'
p='python/codonforge/predictors/cache.py'; s=open(p).read()
old='''    def current(self) -> list[float]:
        if self._sequence is None:'''
new='''    def current(self) -> list[float]:
        if self._sequence is None:
            return []
        if False:'''
assert old in s, "anchor 7"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-10/7: an unprimed cache returns [] instead of refusing" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 $PYT $CACHE -q -k unprimed"

# ---- the diagnostic wiring ----

# 8. The caveat stops being emitted; the numbers still ship.
python3 - <<'PY'
p='python/codonforge/orchestrator/structure_diag.py'; s=open(p).read()
old="    return [detail(value) for value in values], diagnostic_warning(engine_name)"
new="    return [detail(value) for value in values], {'code': 'CF-NONE', 'stage': 'x', 'subject': 'x', 'detail': '', 'title': ''}"
assert old in s, "anchor 8"; open(p,'w').write(s.replace(old,new,1))
PY
control "wiring/8: the front carries the numbers without the caveat" \
  ".venv/bin/python -m pytest python/tests/test_cli.py -q -k dep_stru_003"

# 9. The diagnostic is promoted into the ranked objective order.
python3 - <<'PY'
p='python/codonforge/orchestrator/design.py'; s=open(p).read()
old="        objective_order=list(job.plan.scored),"
new="        objective_order=list(job.plan.scored) + ['structure.local_accessibility'],"
assert old in s, "anchor 9"; open(p,'w').write(s.replace(old,new,1))
PY
control "wiring/9: a term nothing optimized appears in the objective order" \
  ".venv/bin/python -m pytest python/tests/test_cli.py -q -k dep_stru_003"

# 10. The diagnostic reaches only the first candidate.
python3 - <<'PY'
p='python/codonforge/orchestrator/design.py'; s=open(p).read()
old="                    structure_details[position] if structure_details else None"
new="                    structure_details[0] if structure_details and position == 0 else None"
assert old in s, "anchor 10"; open(p,'w').write(s.replace(old,new,1))
PY
control "wiring/10: only the first candidate gets the diagnostic" \
  ".venv/bin/python -m pytest python/tests/test_cli.py -q -k dep_stru_003"

# ---- planning ----

# 11. The sixth rule is removed from the test plan.
python3 - <<'PY'
p='planning/05-testing/01-test-plan.md'; s=open(p).read()
i=s.index('| **A binding to an external tool needs more than reference agreement.**')
j=s.index('\n', i)
open(p,'w').write(s[:i] + s[j+1:])
PY
control "planning/11: the sixth fixture rule disappears from §1a" \
  "grep -q 'binding to an external tool needs more' planning/05-testing/01-test-plan.md"

# 12. DEP-STRU-003 cited but unregistered.
python3 - <<'PY'
p='planning/00-scope/04-departure-index.md'; s=open(p).read()
i=s.index('| `DEP-STRU-003` |'); j=s.index('\n', i)
open(p,'w').write(s[:i] + s[j+1:])
PY
control "planning/12: a cited departure loses its definition" \
  ".venv/bin/python planning/tools/validate_planning_package.py"

# 13. The release gate must still block — confirm it was not weakened this round.
control "gate/13: the release-clearance gate still refuses" \
  ".venv/bin/python tools/check_release_clearance.py"

# 14. The AMD-013 grammar guard against this round's new code.
python3 - <<'PY'
p='python/codonforge/predictors/cache.py'; s=open(p).read()
old="class AccessibilityCache:"
new="# Uses an intersected CFG with bifurcation rules over the codon lattice.\nclass AccessibilityCache:"
assert old in s, "anchor 14"; open(p,'w').write(s.replace(old,new,1))
PY
control "AMD-013/14: the excluded pattern appears in this round's new code" \
  ".venv/bin/python -m pytest python/tests/test_repo_invariants.py -q -k amd_013"

# 15. The dependency boundary still holds after this round's Python additions.
python3 - <<'PY'
p='pyproject.toml'; s=open(p).read()
old='dependencies = [\n  "biopython>=1.83",'
new='dependencies = [\n  "ViennaRNA>=2.6",\n  "biopython>=1.83",'
assert old in s, "anchor 15"; open(p,'w').write(s.replace(old,new,1))
PY
control "dep/15: the cache work quietly makes ViennaRNA a runtime dependency" \
  ".venv/bin/python -m pytest python/tests/test_repo_invariants.py -q -k forbidden"

echo "=================================================================="
echo "controls fired: $pass    silent: $fail"
git checkout -- . 2>/dev/null
exit 0
