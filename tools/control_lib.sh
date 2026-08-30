#!/usr/bin/env bash
# Shared harness for negative-control scripts.
#
# Part B. Last round this harness silently measured a *stale* build for seven
# controls: the Rust defect was injected, the Python test ran against the
# previously-compiled extension, the test passed, and the control reported
# "SILENT — the check did not fire". Read at face value that says the checks are
# redundant. It was the harness that was broken.
#
# That is the same failure `test_the_compiled_extension_is_not_older_than_its_
# sources` exists to catch for pytest, arriving in the tool whose job is
# catching failures. So the invariant is applied here too, and it fails loudly
# rather than producing a wrong conclusion quietly.

set -u

CONTROL_PASS=0
CONTROL_FAIL=0

#: The compiled extension a Python test would import.
_control_extension() {
  echo "python/codonforge/_core.abi3.so"
}

# Fail the whole run if the extension is older than any Rust source.
#
# Deliberately fatal rather than a warning. A warning in a script whose output
# is a list of SILENT/CAUGHT lines is exactly the thing a reader skims past, and
# the cost of skimming past this one is concluding that working checks are
# redundant and deleting them.
control_require_fresh_build() {
  local extension
  extension="$(_control_extension)"

  if [[ ! -f "$extension" ]]; then
    echo "HARNESS ERROR: $extension is missing; run 'make build'." >&2
    exit 3
  fi

  local newest
  newest="$(find crates -name '*.rs' -not -path '*/target/*' -newer "$extension" -print -quit 2>/dev/null)"
  if [[ -n "$newest" ]]; then
    echo "==================================================================" >&2
    echo "HARNESS ERROR: the compiled extension is older than its sources." >&2
    echo "  newer than $extension: $newest" >&2
    echo >&2
    echo "  Every Python control would run against the previous build, so an" >&2
    echo "  injected Rust defect would not be present in the code under test." >&2
    echo "  Controls would report SILENT and the honest reading of that is" >&2
    echo "  'these checks are redundant', which would be wrong and expensive." >&2
    echo "  This happened for seven controls in round 18." >&2
    echo >&2
    echo "  Run 'make build' first." >&2
    echo "==================================================================" >&2
    exit 3
  fi
}

# Run one control: inject (already done by the caller), rebuild, measure, revert.
#
# `rebuild` is not optional and not the caller's responsibility. Making it a
# step each control script has to remember is how round 18 happened.
control() {
  local name="$1"; shift
  local cmd="$1"; shift

  echo "=================================================================="
  echo "CONTROL: $name"

  if ! make build >/tmp/control-build.out 2>&1; then
    # A build failure after injection is not a caught defect. It means the
    # injection was malformed, and reporting it as CAUGHT would credit the
    # check for catching the harness's own mistake.
    echo "  HARNESS ERROR: the build failed after injection, so nothing was"
    echo "  measured. The injected code probably does not compile; fix the"
    echo "  injection so the *test* does the catching."
    tail -5 /tmp/control-build.out | sed 's/^/    /'
    CONTROL_FAIL=$((CONTROL_FAIL + 1))
    git checkout -- . 2>/dev/null
    return
  fi

  if eval "$cmd" >/tmp/control.out 2>&1; then
    echo "  SILENT — the check did not fire"
    CONTROL_FAIL=$((CONTROL_FAIL + 1))
  else
    echo "  CAUGHT:"
    grep -E 'AssertionError|FAILED|panicked|assertion|BLOCKED|Error' /tmp/control.out \
      | head -3 | sed 's/^/    /'
    CONTROL_PASS=$((CONTROL_PASS + 1))
  fi

  git checkout -- . 2>/dev/null
}

control_summary() {
  echo "=================================================================="
  echo "controls fired: $CONTROL_PASS    silent: $CONTROL_FAIL"
  git checkout -- . 2>/dev/null
  make build >/dev/null 2>&1
}
