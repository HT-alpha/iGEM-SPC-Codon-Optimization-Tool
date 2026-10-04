#!/usr/bin/env bash
# Negative controls for round 21 (CF-SRCH-008 sole gate, Crockford ids).
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

# 1. THE REMAINING GATE: Tier2Plan::check no longer refuses R > B.
#    The Python copy is gone (DEP-ORCH-002). If this is silent, the
#    production-path test is not hitting Rust.
python3 - <<'PYEOF'
p='crates/cf-search/src/tier2.rs'; s=open(p).read()
old='                    if codons > config.block {'
new='                    if false && codons > config.block {'
assert old in s, "anchor 1"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/1: Tier2Plan::check no longer refuses a term wider than the block" \
  "$PY python/tests/test_orchestrator.py -q -k wider_than_the_block"

# 2. Encoder alphabet includes I (not in the schema class). TST-CLI-013
#    must catch an out-of-alphabet character, not spot-check one lucky id.
python3 - <<'PYEOF'
p='python/codonforge/orchestrator/design.py'; s=open(p).read()
old='_CROCKFORD = "0123456789ABCDEFGHJKMNPQRSTVWXYZ"'
new='_CROCKFORD = "0123456789ABCDEFGHJKMNPQRSTVWXZI"'
assert old in s, "anchor 2"; open(p,'w').write(s.replace(old,new,1))
PYEOF
control "R/2: _CROCKFORD emits I, outside the schema class" \
  "$PY python/tests/test_derived_id.py -q"

control_summary
