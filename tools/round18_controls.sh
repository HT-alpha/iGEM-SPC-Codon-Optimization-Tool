#!/usr/bin/env bash
# Negative controls for round 18 (P2-05a).
set -u
cd "$(dirname "$0")/.."
pass=0; fail=0
control() {
  local name="$1"; shift; local cmd="$1"; shift
  echo "=================================================================="
  echo "CONTROL: $name"
  # Rebuild *after* injecting. Python tests run against the compiled extension,
  # so without this the control measures the previous build and goes silent -
  # which is the exact failure `test_the_compiled_extension_is_not_older_than_
  # its_sources` exists to catch, reproduced here in the harness itself.
  make build >/dev/null 2>&1
  if eval "$cmd" >/tmp/c18.out 2>&1; then
    echo "  SILENT — the check did not fire"; fail=$((fail+1))
  else
    echo "  CAUGHT:"; grep -E 'AssertionError|FAILED|panicked|assertion|Error' /tmp/c18.out | head -3 | sed 's/^/    /'
    pass=$((pass+1))
  fi
  git checkout -- . 2>/dev/null; make build >/dev/null 2>&1
}
RS="cargo test -p cf-search --test batch_consumption"
PY=".venv/bin/python -m pytest python/tests/test_batch_crossing.py -q"

# (a) REQUIRED: a batch call returning stale scores for the wrong genotypes.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier3.rs'; s=open(p).read()
old="        for (index, value) in chosen.iter().zip(&values) {\n            population[*index].objectives[axis] = *value;"
new="        for (index, value) in chosen.iter().rev().zip(&values) {\n            population[*index].objectives[axis] = *value;"
assert old in s, "anchor a"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(a) results attached to the wrong genotypes" "$RS 2>&1"

# (b) REQUIRED: more than one crossing per generation.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier3.rs'; s=open(p).read()
old="        let values = scorer.score_batch(&genotypes)?;"
new="        let _warm = scorer.score_batch(&genotypes)?;\n        let values = scorer.score_batch(&genotypes)?;"
assert old in s, "anchor b"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(b) a second crossing sneaks into the same generation" "$PY"

# (b2) the per-individual loop this ticket replaced.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier3.rs'; s=open(p).read()
old="        let values = scorer.score_batch(&genotypes)?;"
new="        let mut values = Vec::new();\n        for g in &genotypes {\n            values.extend(scorer.score_batch(std::slice::from_ref(g))?);\n        }"
assert old in s, "anchor b2"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(b2) a per-individual loop instead of one batch" "$PY"

# (c) REQUIRED: CF-SRCH-016 fails to fire when the real path strands someone.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier3.rs'; s=open(p).read()
old="        stranded: request.l0.map(|_| protocol.stranded()).unwrap_or_default(),"
new="        stranded: Vec::new(),"
assert old in s, "anchor c"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(c) the stranded set stops reaching the report" "$PY"

# (c2) the stranded set stops travelling in the checkpoint.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier3.rs'; s=open(p).read()
old="            batch: request.l0.map(|_| protocol.snapshot()),"
new="            batch: None,"
assert old in s, "anchor c2"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(c2) protocol state dropped from the checkpoint" "$RS 2>&1"

# argument fidelity: the scorer stops being given the real population.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier3.rs'; s=open(p).read()
old="            .map(|index| &population[*index].genotype)"
new="            .map(|_| &population[0].genotype)"
assert old in s, "anchor d"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(d) every batch member becomes the same genotype" "$RS 2>&1"

# plausibility: an unscored individual gets a default instead of UNKNOWN.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier3.rs'; s=open(p).read()
old="            population[*index].objectives[axis] = *value;"
new="            population[*index].objectives[axis] = Some(value.unwrap_or(0.0));"
assert old in s, "anchor e"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(e) a declined item is imputed as zero" "$PY"

# a wrong-length reply is zipped instead of refused. Targets BatchProtocol,
# which is now the only mechanism enforcing this (the duplicate in PyL0Scorer
# was removed when this control went silent against it).
python3 - <<'PYEOF'
p='crates/cf-landscape-api/src/batch.rs'; s=open(p).read()
old="        if results.values.len() != expected {"
new="        if false {"
assert old in s, "anchor f"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(f) a short reply is aligned by position rather than refused" \
  "$PY -k malformed"

# the L0 axis silently overwrites an in-core objective.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier3.rs'; s=open(p).read()
old="        if axis < request.evaluator.objective_count() {"
new="        if false {"
assert old in s, "anchor g"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(g) an L0 scorer may overwrite an in-core axis" "$RS 2>&1"

# scoring after selection rather than before.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier3.rs'; s=open(p).read()
old="        if let Some(scorer) = request.l0 {\n            score_l0_batch(&mut combined, scorer, &mut protocol, &mut rationer, &mut rng)?;\n        }\n\n        members = select(combined, request.config.population, selection);"
new="        members = select(combined, request.config.population, selection);\n        if let Some(scorer) = request.l0 {\n            score_l0_batch(&mut members, scorer, &mut protocol, &mut rationer, &mut rng)?;\n        }"
assert old in s, "anchor h"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "(h) selection ranks before the L0 axis is filled" "$RS 2>&1"

echo "=================================================================="
echo "controls fired: $pass    silent: $fail"
git checkout -- . 2>/dev/null; make build >/dev/null 2>&1
exit 0
