#!/usr/bin/env bash
# Negative control for the Baldwinian repair path.
# Forcing the policy read to return true makes both arms Lamarckian, so the
# repair counts match and the test must fail.
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

python3 - <<'PYEOF' || exit 1
p = "crates/cf-search/src/tier3.rs"
s = open(p).read()
old = """fn repairs_are_lamarckian(request: &Tier3Request<'_>) -> bool {
    request.checker.set().repair_policy.lamarckian
}
"""
new = """fn repairs_are_lamarckian(request: &Tier3Request<'_>) -> bool {
    let _ = request;
    true
}
"""
assert old in s, "anchor 1"
assert s.count(old) == 1, "policy read is not unique"
open(p, "w").write(s.replace(old, new, 1))
PYEOF
control "R/1: Baldwinian repair counts match Lamarckian" \
  "$PY python/tests/test_exp05_arms.py -q"

control_summary
