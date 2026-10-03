#!/usr/bin/env bash
# Negative controls for the corpus accession guard and the protein-body
# refusal. Each injection is the branch the check actually takes.
set -u
cd "$(dirname "$0")/.."
source tools/control_lib.sh

control_require_fresh_build

PY=".venv/bin/python -m pytest"

# 1. ACCESSION. A panel gene whose header is not NC_000913.3 is refused.
#    Disabling that predicate lets an invented FASTA through select_genes.
python3 - <<'PYEOF' || exit 1
p = "python/codonforge/bench/corpus.py"
s = open(p).read()
old = """    if _ACCESSION not in header or "[protein_id=" not in header or "[gene=" not in header:
"""
new = """    if False and (_ACCESSION not in header or "[protein_id=" not in header or "[gene=" not in header):
"""
assert old in s, "anchor 1"
assert s.count(old) == 1, "accession guard is not unique"
open(p, "w").write(s.replace(old, new, 1))
PYEOF
control "R/1: a CDS without NC_000913.3 enters the panel" \
  "$PY python/tests/test_corpus.py -q -k sequence_without_the_accession"

# 2. PROTEIN BODY. input.kind protein is CF-CLI-003, not a typo and not a
#    design. Emptying the kind set lets the body fall through to DesignRequest
#    parsing, which reports CF-CLI-001.
python3 - <<'PYEOF' || exit 1
p = "python/codonforge/api/server.py"
s = open(p).read()
old = """        if kind in {"protein", "construct", "single_cds"}:
"""
new = """        if kind in set():
"""
assert old in s, "anchor 2"
assert s.count(old) == 1, "protein-kind refusal is not unique"
open(p, "w").write(s.replace(old, new, 1))
PYEOF
control "R/2: a protein body is treated as a malformed DesignRequest" \
  "$PY python/tests/test_api.py -q -k protein_body"

control_summary
