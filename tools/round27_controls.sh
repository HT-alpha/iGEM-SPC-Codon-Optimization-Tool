#!/usr/bin/env bash
# Negative control for the shipped-corpus digest check.
# The accession guard does not see a replaced gzip whose header still says
# NC_000913.3. This check is the one that does.
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

python3 - <<'PYEOF' || exit 1
p = "python/codonforge/bench/corpus.py"
s = open(p).read()
old = """    if digest != expected:
"""
new = """    if False and digest != expected:
"""
assert old in s, "anchor 1"
assert s.count(old) == 1, "digest comparison is not unique"
open(p, "w").write(s.replace(old, new, 1))
PYEOF
control "R/1: a replaced corpus gzip is accepted" \
  "$PY python/tests/test_corpus.py -q -k replaced_corpus_file"

control_summary
