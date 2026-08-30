#!/usr/bin/env bash
# Negative controls for round 15 (P2-07a, P2-06, P2-13).
#
# Each control injects a defect, runs the check that should catch it, and
# reverts. A control that stays silent means the check or the fixture is
# inadequate — per standing rule, assume the fixture first.
set -u
cd "$(dirname "$0")/.."

pass=0
fail=0

control() {
  local name="$1"; shift
  local cmd="$1"; shift
  echo "=================================================================="
  echo "CONTROL: $name"
  if eval "$cmd" >/tmp/control.out 2>&1; then
    echo "  SILENT — the check did not fire"
    fail=$((fail + 1))
  else
    echo "  CAUGHT:"
    grep -E 'FAILED|panicked at|assertion|error|FAIL ' /tmp/control.out | head -4 | sed 's/^/    /'
    pass=$((pass + 1))
  fi
  git checkout -- . 2>/dev/null
}

# ---------------------------------------------------------------- P2-07a ----

# 1. A parameter set with no licence is accepted.
python3 - <<'PY'
p='crates/cf-biophys-structure/src/parameters.rs'
s=open(p).read()
old = '''            licence: provenance_text(
                "licence",
                "AMD-009 exists because a licence nobody checked turned out to forbid the \\
                 distribution the project had assumed",
            )?,'''
new = '''            licence: provenance_map
                .get("licence")
                .and_then(Value::as_str)
                .unwrap_or("unknown")
                .to_string(),'''
assert old in s, "anchor 1 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-07a/1: an unattributed parameter set is accepted" \
  "cargo test -p cf-biophys-structure 2>&1"

# 2. The digest ignores the parameter content (digests the id only).
python3 - <<'PY'
p='crates/cf-biophys-structure/src/parameters.rs'
s=open(p).read()
old='        let computed = Digest::of_document_omitting(document, &["digest", "created_at"])?;'
new='        let computed = Digest::of_document_omitting(document, &["digest", "created_at", "stack", "temperature_k"])?;'
assert old in s, "anchor 2 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-07a/2: the digest stops covering the energies" \
  "cargo test -p cf-biophys-structure 2>&1"

# 3. A missing stacking energy defaults to zero instead of refusing.
python3 - <<'PY'
p='crates/cf-biophys-structure/src/parameters.rs'
s=open(p).read()
old='        self.stack.get(pair).copied().ok_or_else(|| {'
new='        if true { return Ok(self.stack.get(pair).copied().unwrap_or(0.0)); }\n        self.stack.get(pair).copied().ok_or_else(|| {'
assert old in s, "anchor 3 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-07a/3: an absent stack silently becomes zero" \
  "cargo test -p cf-biophys-structure 2>&1"

# 4. A parameter table gets compiled in.
python3 - <<'PY'
p='crates/cf-biophys-structure/src/parameters.rs'
s=open(p).read()
old='use std::collections::BTreeMap;'
new='use std::collections::BTreeMap;\n\nconst STACK_AU_UA: f64 = -1.10;'
assert old in s, "anchor 4 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-07a/4: a baked-in energy appears in the kernel" \
  "cargo test -p cf-biophys-structure 2>&1"

# ----------------------------------------------------------------- P2-06 ----

# 5. Truncation discards deferred individuals instead of carrying them.
python3 - <<'PY'
p='crates/cf-search/src/rationing.rs'
s=open(p).read()
old='        self.deferred = deferred.clone();'
new='        self.deferred = Vec::new();'
assert old in s, "anchor 5 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-06/5: deferral becomes discarding" \
  "cargo test -p cf-search --test rationing 2>&1"

# 6. The mandatory uniform sample is dropped when the budget is tight.
python3 - <<'PY'
p='crates/cf-search/src/rationing.rs'
s=open(p).read()
old='        let sample_size = ((population as f64) * self.sample_fraction).round() as usize;'
new='        let sample_size = 0usize;'
assert old in s, "anchor 6 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-06/6: the mandatory 10 % sample is dropped" \
  "cargo test -p cf-search --test rationing 2>&1"

# 7. The surrogate predicts before it is fit, returning zero.
python3 - <<'PY'
p='crates/cf-search/src/rationing.rs'
s=open(p).read()
old='        if !self.is_trained() {\n            return Err(CfError::new('
new='        if !self.is_trained() {\n            return Ok(0.0);\n        }\n        if false {\n            return Err(CfError::new('
assert old in s, "anchor 7 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-06/7: an untrained surrogate returns a number" \
  "cargo test -p cf-search --test rationing 2>&1"

# 8. Stranded individuals stop being reported.
python3 - <<'PY'
p='crates/cf-search/src/rationing.rs'
s=open(p).read()
old='        for &index in &deferred {\n            self.never_evaluated.insert(index);\n        }'
new='        for &index in &deferred {\n            let _ = index;\n        }'
assert old in s, "anchor 8 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-06/8: the budget-stranded set goes untracked" \
  "cargo test -p cf-search --test rationing 2>&1"

# 9. The deferral warning stops reaching job_status.
python3 - <<'PY'
p='python/codonforge/orchestrator/job.py'
s=open(p).read()
old='        if not stranded:\n            return []'
new='        if True:\n            return []'
assert old in s, "anchor 9 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-06/9: CF-SRCH-016 stops being emitted" \
  ".venv/bin/python -m pytest python/tests/test_orchestrator.py -q -k deferral 2>&1"

# ----------------------------------------------------------------- P2-13 ----

# 10. Tier 3 claims calibration.
python3 - <<'PY'
p='crates/cf-burden/src/capacity.rs'
s=open(p).read()
old='        calibrated: false,\n        rationale_tags: vec!["low_confidence_host"],'
new='        calibrated: true,\n        rationale_tags: vec!["low_confidence_host"],'
assert old in s, "anchor 10 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-13/10: tier 3 claims to be calibrated" \
  "cargo test -p cf-burden --test capacity 2>&1"

# 11. Tier 3 reports a growth-rate unit.
python3 - <<'PY'
p='crates/cf-burden/src/capacity.rs'
s=open(p).read()
old='            Self::FractionOfNativeDemand => "dimensionless",'
new='            Self::FractionOfNativeDemand => "relative_growth_rate",'
assert old in s, "anchor 11 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-13/11: tier 3 dresses a mismatch score as a growth prediction" \
  "cargo test -p cf-burden --test capacity 2>&1"

# 12. THE ONE THE AUDIT ASKED FOR: the estimator ignores the host's data and
#     returns a fixed profile. Every semantic check still passes; only
#     TST-BURD-006 can see it.
python3 - <<'PY'
p='crates/cf-burden/src/capacity.rs'
s=open(p).read()
old='''    let supply: BTreeMap<[u8; 3], f64> = native_demand
        .iter()
        .map(|(anticodon, demand)| (*anticodon, demand / total))
        .collect();'''
new='''    let count = native_demand.len() as f64;
    let supply: BTreeMap<[u8; 3], f64> = native_demand
        .keys()
        .map(|anticodon| (*anticodon, 1.0 / count))
        .collect();'''
assert old in s, "anchor 12 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-13/12: the estimator returns a plausible constant, ignoring the host" \
  "cargo test -p cf-burden --test capacity 2>&1"

# 13. Tier 1 accepts an inadequate panel.
python3 - <<'PY'
p='crates/cf-burden/src/capacity.rs'
s=open(p).read()
old='    if panel.len() < 20 {'
new='    if panel.len() < 2 {'
assert old in s, "anchor 13 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-13/13: tier 1 accepts a 12-construct panel" \
  "cargo test -p cf-burden --test capacity 2>&1"

# 14. Tier selection claims more than the data supports.
python3 - <<'PY'
p='crates/cf-burden/src/capacity.rs'
s=open(p).read()
old='        if self.has_burden_panel {\n            1\n        } else if self.has_ribosome_count {'
new='        if self.has_burden_panel || self.has_ribosome_count {\n            1\n        } else if self.has_ribosome_count {'
assert old in s, "anchor 14 missing"
open(p,'w').write(s.replace(old,new,1))
PY
control "P2-13/14: ribosome counts alone buy tier 1" \
  "cargo test -p cf-burden --test capacity 2>&1"

# ------------------------------------------------------------ validators ----

# 15. The planning validator catches an unregistered digest.
python3 - <<'PY'
import json
p='planning/02-schemas/examples/structure_parameters.example.json'
d=json.load(open(p))
d['digest']='sha256:'+'f'*64
json.dump(d, open(p,'w'), indent=2)
PY
control "validator/15: a wrong example digest" \
  ".venv/bin/python planning/tools/validate_planning_package.py 2>&1"

# 16. The schema-sync check catches planning/shipped drift.
python3 - <<'PY'
import json
p='planning/02-schemas/structure_parameters.schema.json'
d=json.load(open(p))
d['properties']['temperature_k']['description']='drift'
json.dump(d, open(p,'w'), indent=2)
PY
control "validator/16: planning and shipped schemas drift apart" \
  ".venv/bin/python tools/check_schema_sync.py --check 2>&1"

echo "=================================================================="
echo "controls fired: $pass    silent: $fail"
git checkout -- . 2>/dev/null
exit 0
