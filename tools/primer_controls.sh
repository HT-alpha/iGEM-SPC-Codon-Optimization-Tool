#!/bin/sh
# Negative controls for the primer guards. Each injection must make a test fail.
# Reverts only the files it edits. Run from a clean tree.
set -eu
cd "$(dirname "$0")/.."
fail() { echo "CONTROL FAILED: $1" >&2; exit 1; }

echo "control 1: digest guard forced open"
python3 - <<'PY'
from pathlib import Path
path = Path("crates/cf-primer-thermo/src/parameters.rs")
text = path.read_text()
old = '''    if computed.as_str() != digest {
        return Err(ParameterError(format!(
            "CF-PRIM-005 primer parameter digest {digest} does not match {computed}"
        )));
    }'''
new = "    let _ = (computed, digest);"
if old not in text:
    raise SystemExit("digest guard not found")
path.write_text(text.replace(old, new, 1))
PY
set +e
cargo test -p cf-primer-thermo --lib a_replaced_parameter_document_is_refused --offline >/tmp/primer-control-digest.log 2>&1
code=$?
set -e
git checkout -- crates/cf-primer-thermo/src/parameters.rs
if [ "$code" -eq 0 ]; then
  fail "digest guard injection was silent"
fi
echo "CAUGHT digest guard (exit $code)"

echo "control 2: homopolymer filter forced open"
python3 - <<'PY'
from pathlib import Path
path = Path("crates/cf-primer-search/src/lib.rs")
text = path.read_text()
old = "longest_homopolymer(seq) < 5"
new = "longest_homopolymer(seq) < 100"
if old not in text:
    raise SystemExit("homopolymer filter not found")
path.write_text(text.replace(old, new, 1))
PY
set +e
cargo test -p cf-primer-search --lib each_prune_filter_removes_what_it_should --offline >/tmp/primer-control-homo.log 2>&1
code=$?
set -e
git checkout -- crates/cf-primer-search/src/lib.rs
if [ "$code" -eq 0 ]; then
  fail "homopolymer filter injection was silent"
fi
echo "CAUGHT homopolymer filter (exit $code)"

echo "control 2b: GC filter forced open"
python3 - <<'PY'
from pathlib import Path
path = Path("crates/cf-primer-search/src/lib.rs")
text = path.read_text()
old = "(30.0..=70.0).contains(&gc)"
new = "(0.0..=100.0).contains(&gc)"
if old not in text:
    raise SystemExit("gc filter not found")
path.write_text(text.replace(old, new, 1))
PY
set +e
cargo test -p cf-primer-search --lib each_prune_filter_removes_what_it_should --offline >/tmp/primer-control-gc.log 2>&1
code=$?
set -e
git checkout -- crates/cf-primer-search/src/lib.rs
if [ "$code" -eq 0 ]; then
  fail "GC filter injection was silent"
fi
echo "CAUGHT GC filter (exit $code)"

echo "control 2c: Tm prefilter forced open"
python3 - <<'PY'
from pathlib import Path
path = Path("crates/cf-primer-search/src/lib.rs")
text = path.read_text()
old = "(tm - target).abs() <= 8.0"
new = "(tm - target).abs() <= 1000.0"
if old not in text:
    raise SystemExit("tm filter not found")
path.write_text(text.replace(old, new, 1))
PY
set +e
cargo test -p cf-primer-search --lib each_prune_filter_removes_what_it_should --offline >/tmp/primer-control-tm.log 2>&1
code=$?
set -e
git checkout -- crates/cf-primer-search/src/lib.rs
if [ "$code" -eq 0 ]; then
  fail "Tm filter injection was silent"
fi
echo "CAUGHT Tm filter (exit $code)"

echo "control 2d: product-size filter forced open"
python3 - <<'PY'
from pathlib import Path
path = Path("crates/cf-primer-search/src/lib.rs")
text = path.read_text()
old = "(min..=max).contains(&size)"
new = "(0..=usize::MAX).contains(&size)"
if old not in text:
    raise SystemExit("product filter not found")
path.write_text(text.replace(old, new, 1))
PY
set +e
cargo test -p cf-primer-search --lib each_prune_filter_removes_what_it_should --offline >/tmp/primer-control-product.log 2>&1
code=$?
set -e
git checkout -- crates/cf-primer-search/src/lib.rs
if [ "$code" -eq 0 ]; then
  fail "product-size filter injection was silent"
fi
echo "CAUGHT product-size filter (exit $code)"

echo "control 3: refuse zone forced closed"
python3 - <<'PY'
from pathlib import Path
path = Path("crates/cf-primer-thermo/src/duplex.rs")
text = path.read_text()
old = "if delta_g < -9.0 {"
new = "if false && delta_g < -9.0 {"
if old not in text:
    raise SystemExit("refuse gate not found")
path.write_text(text.replace(old, new, 1))
PY
set +e
cargo test -p cf-primer-thermo --lib dg_zones_match_the_idt_cutoffs --offline >/tmp/primer-control-dg.log 2>&1
code=$?
set -e
git checkout -- crates/cf-primer-thermo/src/duplex.rs
if [ "$code" -eq 0 ]; then
  fail "delta-G refuse injection was silent"
fi
echo "CAUGHT delta-G refuse gate (exit $code)"
echo "ALL PRIMER CONTROLS CAUGHT"
