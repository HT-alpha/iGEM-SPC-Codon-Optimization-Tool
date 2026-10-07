#!/bin/sh
# Negative control for the codon-usage table.
#
# Omitting a codon whose count is zero is the AGG gap: the table then has no
# entry, which is not the same record as relative adaptiveness 0. The
# completeness test must fail. Reverts only the file it edits.
set -eu
cd "$(dirname "$0")/.."
fail() { echo "CONTROL FAILED: $1" >&2; exit 1; }

python3 - <<'PY'
import sys
sys.path.insert(0, "tools")
import patch

patch.sub(
    "python/codonforge/ingest/host.py",
    """    for codon in _ALL_CODONS:
        count = counts.get(codon, 0)
""",
    """    for codon in _ALL_CODONS:
        count = counts.get(codon, 0)
        if count == 0:
            continue
""",
)
PY

set +e
.venv/bin/pytest -q -p no:cacheprovider \
  python/tests/test_ingest.py -k every_codon \
  >/tmp/usage-control.log 2>&1
code=$?
set -e
git checkout -- python/codonforge/ingest/host.py
if [ "$code" -eq 0 ]; then
  fail "omitted-codon injection was silent"
fi
echo "CAUGHT omitted codon (exit $code)"
grep -E 'AssertionError|FAILED' /tmp/usage-control.log | head -4
