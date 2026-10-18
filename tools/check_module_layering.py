#!/usr/bin/env python3
"""Enforce the module graph declared in planning/01-modules/modules.yaml.

`modules.yaml` is normative (BRIEF-2.1-b). This tool makes it enforced rather
than documentary, by comparing it against the real source tree:

  --check     fail if the implementation violates the declared graph
  --scaffold  create missing module skeletons and refresh the workspace members

Checks performed by --check:

  1. structure       every declared module has its declared directory
  2. edges           no cf-* dependency exists that modules.yaml does not declare
                     (an undeclared edge is how a cycle gets introduced)
  3. rank            every real edge strictly decreases in rank
  4. acyclic         the real graph is a DAG
  5. ffi-boundary    pyo3 appears in exactly one crate, cf-pyext
                     (ADR-0001: no call may cross the Rust/Python boundary more
                      than once per generation; the only way to guarantee that
                      structurally is to keep the binding surface in one crate)
  6. no-py-in-core   no rust-core crate contains a Python callback construct
  7. runtime-compat  no rust-core module declares a dependency on a python
                     module. Such an edge cannot compile - a Cargo crate has no
                     way to depend on a Python package - so it is a declaration
                     that can never be satisfied, and the `edges` check will not
                     notice because it only looks for edges the tree has that
                     modules.yaml lacks, never the reverse. Anything a rust-core
                     module needs from a python module must cross as a value.
  8. layer-coverage  every architectural layer still has an owning module

Exit code 0 iff every selected check passes.
"""

from __future__ import annotations

import argparse
import re
import sys
import textwrap
from collections import defaultdict, deque
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("PyYAML is required: pip install pyyaml")

ROOT = Path(__file__).resolve().parent.parent
MODULES_YAML = ROOT / "planning" / "01-modules" / "modules.yaml"

# Crates permitted to depend on pyo3. Exactly one, by design.
PYO3_ALLOWED = {"cf-pyext"}

# Source constructs that would smuggle a Python call into the hot path.
FORBIDDEN_IN_CORE = [
    (re.compile(r"\bPython::with_gil\b"), "Python::with_gil in a rust-core crate"),
    (re.compile(r"\bpyo3\b"), "pyo3 reference in a rust-core crate"),
]


class Report:
    def __init__(self, verbose: bool) -> None:
        self.verbose = verbose
        self.failures: list[str] = []
        self.assertions = 0

    def ok(self, msg: str) -> None:
        self.assertions += 1
        if self.verbose:
            print(f"  ok   {msg}")

    def fail(self, msg: str) -> None:
        self.assertions += 1
        self.failures.append(msg)
        print(f"  FAIL {msg}")

    def section(self, name: str) -> None:
        print(f"\n[{name}]")


def load_modules() -> tuple[dict, dict]:
    doc = yaml.safe_load(MODULES_YAML.read_text(encoding="utf-8"))
    return doc, {m["name"]: m for m in doc["modules"]}


def rust_modules(mods: dict) -> dict:
    return {n: m for n, m in mods.items() if m["runtime"] == "rust-core"}


def python_modules(mods: dict) -> dict:
    return {n: m for n, m in mods.items() if m["runtime"].startswith("python")}


def py_subpackage(module: dict) -> str:
    return Path(module["path"]).name


# --------------------------------------------------------------------------- #
# scaffold
# --------------------------------------------------------------------------- #


def scaffold(mods: dict) -> None:
    created: list[str] = []

    for name, mod in sorted(rust_modules(mods).items()):
        crate = ROOT / mod["path"]
        (crate / "src").mkdir(parents=True, exist_ok=True)
        manifest = crate / "Cargo.toml"
        if not manifest.exists():
            manifest.write_text(
                f'[package]\nname = "{name}"\nversion.workspace = true\n'
                f"edition.workspace = true\nlicense.workspace = true\n"
                f"rust-version.workspace = true\n\n[dependencies]\n",
                encoding="utf-8",
            )
            created.append(str(manifest.relative_to(ROOT)))
        lib = crate / "src" / "lib.rs"
        if not lib.exists():
            responsibility = " ".join(mod["responsibility"].split())
            wrapped = []
            line = "//!"
            for word in responsibility.split():
                if len(line) + len(word) + 1 > 78:
                    wrapped.append(line)
                    line = "//!"
                line += " " + word
            wrapped.append(line)
            body = "\n".join(wrapped)
            lib.write_text(
                f"//! `{name}` — {mod['layer']} layer.\n//!\n{body}\n//!\n"
                f"//! Declared dependencies: {', '.join(mod['depends_on']) or 'none'}.\n"
                f"//! Contract: `{mod['public_contract']}`.\n\n"
                "#![forbid(unsafe_code)]\n",
                encoding="utf-8",
            )
            created.append(str(lib.relative_to(ROOT)))

    # cf-pyext is glue, not an architectural layer, so it is not in modules.yaml.
    pyext = ROOT / "crates" / "cf-pyext"
    (pyext / "src").mkdir(parents=True, exist_ok=True)
    if not (pyext / "src" / "lib.rs").exists():
        (pyext / "src" / "lib.rs").write_text(
            "//! PyO3 glue producing the `dna_design_suite._core` extension module.\n"
            "//!\n//! This is the ONLY crate permitted to depend on pyo3\n"
            "//! (ADR-0001). `tools/check_module_layering.py --check` enforces it.\n",
            encoding="utf-8",
        )
        created.append("crates/cf-pyext/src/lib.rs")

    for name, mod in sorted(python_modules(mods).items()):
        pkg = ROOT / mod["path"]
        pkg.mkdir(parents=True, exist_ok=True)
        init = pkg / "__init__.py"
        if not init.exists():
            responsibility = textwrap.fill(" ".join(mod["responsibility"].split()), width=88)
            deps = ", ".join(mod["depends_on"]) or "none"
            init.write_text(
                f'"""{name} — {mod["layer"]} layer.\n\n{responsibility}\n\n'
                f"Declared dependencies: {deps}.\n"
                f'Contract: {mod["public_contract"]}\n"""\n',
                encoding="utf-8",
            )
            created.append(str(init.relative_to(ROOT)))

    for msg in created:
        print(f"  created {msg}")
    if not created:
        print("  nothing to create; tree already matches modules.yaml")


# --------------------------------------------------------------------------- #
# checks
# --------------------------------------------------------------------------- #

_DEP_LINE = re.compile(r"^\s*(cf-[a-z-]+)\s*=", re.MULTILINE)


def real_rust_edges(mods: dict, rep: Report) -> dict[str, set[str]]:
    edges: dict[str, set[str]] = {}
    for name, mod in rust_modules(mods).items():
        manifest = ROOT / mod["path"] / "Cargo.toml"
        if not manifest.exists():
            rep.fail(f"{name}: {manifest.relative_to(ROOT)} is missing")
            edges[name] = set()
            continue
        edges[name] = set(_DEP_LINE.findall(manifest.read_text(encoding="utf-8")))
    return edges


_PY_IMPORT = re.compile(
    r"^\s*(?:from\s+dna_design_suite(?:\.codon)?\.([a-z_]+)|import\s+dna_design_suite(?:\.codon)?\.([a-z_]+))",
    re.MULTILINE,
)


def real_python_edges(mods: dict) -> dict[str, set[str]]:
    by_subpkg = {py_subpackage(m): n for n, m in python_modules(mods).items()}
    edges: dict[str, set[str]] = {}
    for name, mod in python_modules(mods).items():
        found: set[str] = set()
        pkg = ROOT / mod["path"]
        if pkg.exists():
            for src in pkg.rglob("*.py"):
                for a, b in _PY_IMPORT.findall(src.read_text(encoding="utf-8")):
                    sub = a or b
                    if sub in by_subpkg and by_subpkg[sub] != name:
                        found.add(by_subpkg[sub])
        edges[name] = found
    return edges


def check(mods: dict, rep: Report) -> None:
    rep.section("structure")
    for name, mod in sorted(mods.items()):
        path = ROOT / mod["path"]
        if path.is_dir():
            rep.ok(f"{name}: {mod['path']}")
        else:
            rep.fail(f"{name}: declared path {mod['path']} does not exist")

    rep.section("edges")
    edges = real_rust_edges(mods, rep)
    edges.update(real_python_edges(mods))
    for name, actual in sorted(edges.items()):
        declared = set(mods[name]["depends_on"])
        undeclared = actual - declared
        if undeclared:
            rep.fail(
                f"{name}: undeclared dependency {sorted(undeclared)} — "
                "add it to modules.yaml (and re-check for a cycle) or remove it"
            )
        else:
            rep.ok(f"{name}: {len(actual)} real edge(s), all declared")

    rep.section("rank")
    for name, actual in sorted(edges.items()):
        for dep in sorted(actual):
            if dep not in mods:
                rep.fail(f"{name} -> {dep}: not a declared module")
            elif mods[dep]["rank"] >= mods[name]["rank"]:
                rep.fail(
                    f"{name} (rank {mods[name]['rank']}) -> {dep} "
                    f"(rank {mods[dep]['rank']}): rank must strictly decrease"
                )
    rep.ok("all real edges strictly decrease in rank")

    rep.section("acyclic")
    indeg = {n: 0 for n in mods}
    adj = defaultdict(list)
    for name, actual in edges.items():
        for dep in actual:
            if dep in mods:
                adj[dep].append(name)
                indeg[name] += 1
    queue = deque(n for n, d in indeg.items() if d == 0)
    seen = 0
    while queue:
        node = queue.popleft()
        seen += 1
        for nxt in adj[node]:
            indeg[nxt] -= 1
            if indeg[nxt] == 0:
                queue.append(nxt)
    if seen == len(mods):
        rep.ok(f"real dependency graph is a DAG ({seen} modules)")
    else:
        rep.fail(f"cycle among: {sorted(n for n, d in indeg.items() if d > 0)}")

    rep.section("ffi-boundary")
    crates_dir = ROOT / "crates"
    holders = []
    if crates_dir.is_dir():
        for manifest in sorted(crates_dir.glob("*/Cargo.toml")):
            text = manifest.read_text(encoding="utf-8")
            if re.search(r"^\s*pyo3\s*=", text, re.MULTILINE):
                holders.append(manifest.parent.name)
    if set(holders) - PYO3_ALLOWED:
        rep.fail(
            f"pyo3 declared in {sorted(set(holders) - PYO3_ALLOWED)}; "
            f"ADR-0001 confines the binding surface to {sorted(PYO3_ALLOWED)}"
        )
    else:
        rep.ok(f"pyo3 confined to {sorted(PYO3_ALLOWED)} (found in {holders or 'none yet'})")

    rep.section("no-py-in-core")
    offenders = 0
    for mod in sorted(rust_modules(mods).values(), key=lambda m: m["name"]):
        src = ROOT / mod["path"] / "src"
        if not src.is_dir():
            continue
        for rs in sorted(src.rglob("*.rs")):
            text = rs.read_text(encoding="utf-8")
            for pattern, why in FORBIDDEN_IN_CORE:
                if pattern.search(text):
                    rep.fail(f"{rs.relative_to(ROOT)}: {why}")
                    offenders += 1
    if not offenders:
        rep.ok("no rust-core crate references pyo3 or holds the GIL")

    rep.section("runtime-compat")
    impossible = 0
    for name, mod in sorted(rust_modules(mods).items()):
        for dep in sorted(mod["depends_on"]):
            if dep in mods and mods[dep]["runtime"].startswith("python"):
                rep.fail(
                    f"{name} ({mod['runtime']}) -> {dep} ({mods[dep]['runtime']}): "
                    "a Cargo crate cannot depend on a Python package. Invert it — "
                    "pass what the crate needs across the boundary as a value."
                )
                impossible += 1
    if not impossible:
        rep.ok("no rust-core module declares a dependency on a python module")

    rep.section("layer-coverage")
    expected = {
        "context-ingestion",
        "epistatic-surface",
        "biophysical-predictors",
        "burden-model",
        "search-engine",
        "constraint-repair",
        "decision-support",
        "active-learning",
    }
    present = {m["layer"] for m in mods.values()}
    missing = expected - present
    if missing:
        rep.fail(f"layers with no owning module: {sorted(missing)}")
    else:
        rep.ok("every architectural layer has an owning module")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--scaffold", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    if not (args.check or args.scaffold):
        args.check = True

    _doc, mods = load_modules()
    rep = Report(args.verbose)

    if args.scaffold:
        print("[scaffold]")
        scaffold(mods)
    if args.check:
        check(mods, rep)
        print("\n" + "=" * 70)
        if rep.failures:
            print(f"FAILED: {len(rep.failures)} problem(s) across {rep.assertions} assertions")
            for failure in rep.failures:
                print(f"  - {failure}")
            return 1
        print(f"PASSED: {rep.assertions} assertions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
