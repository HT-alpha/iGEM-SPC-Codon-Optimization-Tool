#!/usr/bin/env bash
# Negative controls for round 16 (P2-09, P2-08, the release gate, DEP-SRCH-011).
set -u
cd "$(dirname "$0")/.."
pass=0; fail=0
control() {
  local name="$1"; shift; local cmd="$1"; shift
  echo "=================================================================="
  echo "CONTROL: $name"
  if eval "$cmd" >/tmp/c.out 2>&1; then
    echo "  SILENT — the check did not fire"; fail=$((fail+1))
  else
    echo "  CAUGHT:"; grep -E 'FAILED|panicked|assertion|AssertionError|BLOCKED|Error|error' /tmp/c.out | head -3 | sed 's/^/    /'
    pass=$((pass+1))
  fi
  git checkout -- . 2>/dev/null; rm -f docs/legal-clearance.yaml
}

# ---- P2-08: the bugs that produce plausible numbers ----

# 1. Off-by-one in the 1-based extraction.
python3 - <<'PY'
p='python/codonforge/predictors/viennarna.py'; s=open(p).read()
old="        for end in range(window_nt, len(rna) + 1):"
new="        for end in range(window_nt - 1, len(rna)):"
assert old in s, "anchor 1"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-08/1: a one-nucleotide index shift" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 .venv/bin/python -m pytest python/tests/test_structure_engine.py -q"

# 2. Temperature ignored — the binding folds at whatever the library last saw.
python3 - <<'PY'
p='python/codonforge/predictors/viennarna.py'; s=open(p).read()
old="        self._rna.cvar.temperature = float(temperature_c)"
new="        pass  # temperature dropped"
assert old in s, "anchor 2"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-08/2: the temperature argument is ignored" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 .venv/bin/python -m pytest python/tests/test_structure_engine.py -q"

# 3. Impossible windows passed through as 0.0 (reads as maximal inaccessibility).
python3 - <<'PY'
p='python/codonforge/predictors/viennarna.py'; s=open(p).read()
old="        for end in range(window_nt, len(rna) + 1):\n            values.append(float(matrix[end][window_nt]))"
new="        for end in range(1, len(rna) + 1):\n            values.append(float(matrix[end][window_nt]))"
assert old in s, "anchor 3"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-08/3: impossible windows reported as zero rather than absent" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 .venv/bin/python -m pytest python/tests/test_structure_engine.py -q"

# 4. The engine returns a plausible constant instead of folding.
python3 - <<'PY'
p='python/codonforge/predictors/viennarna.py'; s=open(p).read()
old="        values = []"
new="        return [0.5] * (len(rna) - window_nt + 1)\n        values = []"
assert old in s, "anchor 4"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-08/4: a constant that looks like a probability" \
  "CF_REQUIRE_STRUCTURE_ENGINE=1 .venv/bin/python -m pytest python/tests/test_structure_engine.py -q"

# 5. AMD-001: global MFE reappears on the engine the search loop can reach.
python3 - <<'PY'
p='python/codonforge/predictors/viennarna.py'; s=open(p).read()
old="    def local_accessibility("
new="    def global_mfe(self, sequence, temperature_c):\n        return self._rna.fold(sequence)\n\n    def local_accessibility("
assert old in s, "anchor 5"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-08/5: global MFE creeps back onto the per-move path" \
  ".venv/bin/python -m pytest python/tests/test_structure_engine.py -q -k amd_001"

# ---- P2-09: the seam ----

# 6. No engine installed degrades to a proxy instead of refusing.
python3 - <<'PY'
p='python/codonforge/predictors/structure.py'; s=open(p).read()
old='        raise StructureEngineError(\n            "CF-STRU-003",\n            "no structure engine is installed'
new='        return ViennaRNAEngine()  # silent fallback\n        raise StructureEngineError(\n            "CF-STRU-003",\n            "no structure engine is installed'
assert old in s, "anchor 6"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-09/6: absence degrades instead of refusing" \
  ".venv/bin/python -m pytest python/tests/test_structure_engine.py -q -k stru_003"

# 7. Detection imports the package instead of using find_spec.
python3 - <<'PY'
p='python/codonforge/predictors/structure.py'; s=open(p).read()
old="            if importlib.util.find_spec(module) is not None:\n                found.append(engine)"
new="            import RNA  # noqa: F401\n            found.append(engine)"
assert old in s, "anchor 7"; open(p,'w').write(s.replace(old,new,1))
PY
control "P2-09/7: engine detection imports the engine" \
  ".venv/bin/python -m pytest python/tests/test_structure_engine.py -q -k no_import"

# ---- release gate ----

# 8. The gate passes with an unsigned clearance file.
cat > docs/legal-clearance.yaml <<'YAML'
amendment: "AMD-013"
scope: "looks fine to me"
YAML
control "gate/8: an incomplete clearance record is accepted" \
  ".venv/bin/python tools/check_release_clearance.py"

# 9. The gate accepts a sign-off on a different analysis.
cat > docs/legal-clearance.yaml <<'YAML'
signed_off_by: "Someone"
role: "faculty advisor"
date: "2026-06-01"
amendment: "AMD-009"
scope: "reviewed the old amendment"
YAML
control "gate/9: a sign-off on the wrong amendment clears this gate" \
  ".venv/bin/python tools/check_release_clearance.py"

# 10. The gate accepts a stale signature.
cat > docs/legal-clearance.yaml <<'YAML'
signed_off_by: "Someone"
role: "faculty advisor"
date: "2020-01-01"
amendment: "AMD-013"
scope: "reviewed a long time ago"
YAML
control "gate/10: a five-year-old sign-off still clears" \
  ".venv/bin/python tools/check_release_clearance.py"

# 11. The gate is wired into the release workflow, not just runnable.
python3 - <<'PY'
import yaml
p='.github/workflows/wheels.yml'; d=yaml.safe_load(open(p))
assert d['jobs']['build'].get('needs')=='legal-clearance', "build must depend on the gate"
assert 'legal-clearance' in d['jobs'], "gate job must exist"
print("wired")
PY
python3 - <<'PY'
import yaml,io
p='.github/workflows/wheels.yml'; s=open(p).read()
s=s.replace("    needs: legal-clearance\n","",1)
open(p,'w').write(s)
PY
control "gate/11: the release build stops depending on the gate" \
  "python3 -c \"
import yaml,sys
d=yaml.safe_load(open('.github/workflows/wheels.yml'))
sys.exit(0 if d['jobs']['build'].get('needs')=='legal-clearance' else 1)\""

# ---- dependency policy ----

# 12. A structure engine becomes a runtime dependency.
python3 - <<'PY'
p='pyproject.toml'; s=open(p).read()
old='dependencies = [\n  "biopython>=1.83",'
new='dependencies = [\n  "ViennaRNA>=2.6",\n  "biopython>=1.83",'
assert old in s, "anchor 12"; open(p,'w').write(s.replace(old,new,1))
PY
control "dep/12: ViennaRNA becomes a shipped runtime dependency" \
  ".venv/bin/python -m pytest python/tests/test_repo_invariants.py -q -k forbidden"

# 13. LinearDesign appears in an extra (no extra clears a patent).
python3 - <<'PY'
p='pyproject.toml'; s=open(p).read()
old='structure-tests = [\n  "ViennaRNA>=2.6",\n]'
new='structure-tests = [\n  "ViennaRNA>=2.6",\n  "lineardesign",\n]'
assert old in s, "anchor 13"; open(p,'w').write(s.replace(old,new,1))
PY
control "dep/13: LinearDesign smuggled in via a test extra" \
  ".venv/bin/python -m pytest python/tests/test_repo_invariants.py -q -k forbidden"

# 14. A Rust crate takes a structure-engine dependency (the vendoring boundary).
python3 - <<'PY'
p='crates/cf-biophys-structure/Cargo.toml'; s=open(p).read()
open(p,'w').write(s + '\nviennarna-sys = "0.1"\n')
PY
control "dep/14: a Rust crate links a structure engine" \
  ".venv/bin/python -m pytest python/tests/test_repo_invariants.py -q -k forbidden"

# ---- DEP-SRCH-011 ----

# 15. Deferral becomes discarding in the one remaining implementation.
python3 - <<'PY'
p='crates/cf-landscape-api/src/batch.rs'; s=open(p).read()
old="        self.deferred = queue.into_iter().skip(budget).collect();"
new="        self.deferred = Vec::new();"
assert old in s, "anchor 15"; open(p,'w').write(s.replace(old,new,1))
PY
control "DEP-SRCH-011/15: the surviving deferral queue drops its tail" \
  "cargo test -p cf-landscape-api --test batch 2>&1"

# 16. The stranded set stops travelling in the snapshot.
python3 - <<'PY'
p='crates/cf-landscape-api/src/batch.rs'; s=open(p).read()
old='            "stranded": self.stranded.iter().copied().collect::<Vec<usize>>(),'
new='            "stranded": Vec::<usize>::new(),'
assert old in s, "anchor 16"; open(p,'w').write(s.replace(old,new,1))
PY
control "DEP-SRCH-011/16: CF-SRCH-016 falls silent across a resume" \
  "cargo test -p cf-landscape-api --test batch 2>&1"

# ---- AMD-013 grammar guard, against this round's new code ----
python3 - <<'PY'
p='python/codonforge/predictors/viennarna.py'; s=open(p).read()
old="class ViennaRNAEngine:"
new="# Intersected CFG over the codon DFA, with bifurcation rules.\nclass ViennaRNAEngine:"
assert old in s, "anchor 17"; open(p,'w').write(s.replace(old,new,1))
PY
control "AMD-013/17: the excluded pattern appears in this round's new code" \
  ".venv/bin/python -m pytest python/tests/test_repo_invariants.py -q -k amd_013"

echo "=================================================================="
echo "controls fired: $pass    silent: $fail"
git checkout -- . 2>/dev/null; rm -f docs/legal-clearance.yaml
exit 0
