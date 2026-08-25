#!/usr/bin/env python3
"""Consistency checker for the CodonForge planning package.

This is planning-package tooling, not CodonForge backend code. It exists so that
a package made of schemas, graphs and cross-referenced documents can be reviewed
mechanically instead of by eye. It answers, for every check, a question a
reviewer would otherwise have to answer manually:

  schemas       Are the JSON Schemas actually valid draft 2020-12, and do the
                worked examples validate against them?
  deps          Is the module dependency graph a DAG, with no dangling edges and
                no rank inversion?  (BRIEF-2.1-b)
  refs          Does every ADR / AMD / EXP / RISK / GATE / TST identifier cited
                anywhere in the package actually exist?
  links         Does every relative Markdown link resolve to a real file?
  trace         Does every requirement in the assignment map to an artifact that
                exists?  (BRIEF-5 acceptance criteria)
  tags          Are all rationale tags used in examples part of the closed
                vocabulary?

Usage:
    python3 planning/tools/validate_planning_package.py [--check NAME ...] [-v]

Exit code 0 iff every selected check passes.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict, deque
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("PyYAML is required: pip install pyyaml")

try:
    from jsonschema import Draft202012Validator
    from referencing import Registry, Resource
except ImportError:  # pragma: no cover
    sys.exit("jsonschema>=4.18 and referencing are required: pip install jsonschema")

PLANNING = Path(__file__).resolve().parent.parent
ROOT = PLANNING.parent

# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #


class Report:
    def __init__(self, verbose: bool) -> None:
        self.verbose = verbose
        self.failures: list[str] = []
        self.checks_run = 0
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
        self.checks_run += 1
        print(f"\n[{name}]")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def schema_files() -> list[Path]:
    return sorted(PLANNING.glob("02-schemas/*.schema.json"))


def build_registry() -> Registry:
    """Register every local schema under its $id so relative $refs resolve."""
    resources = []
    for path in schema_files():
        doc = load_json(path)
        if "$id" not in doc:
            continue
        resources.append((doc["$id"], Resource.from_contents(doc)))
    return Registry().with_resources(resources)


def markdown_files() -> list[Path]:
    return sorted(PLANNING.rglob("*.md"))


def all_text_files() -> list[Path]:
    out = []
    for pattern in ("*.md", "*.yaml", "*.yml", "*.json", "*.pyi"):
        out.extend(PLANNING.rglob(pattern))
    return sorted(p for p in out if "tools" not in p.parts)


# --------------------------------------------------------------------------- #
# checks
# --------------------------------------------------------------------------- #


def check_schemas(rep: Report) -> None:
    rep.section("schemas")
    registry = build_registry()

    for path in schema_files():
        doc = load_json(path)
        try:
            Draft202012Validator.check_schema(doc)
            rep.ok(f"{path.relative_to(ROOT)} is valid draft 2020-12")
        except Exception as exc:  # noqa: BLE001
            rep.fail(f"{path.relative_to(ROOT)} is not a valid schema: {exc}")

    example_dir = PLANNING / "02-schemas" / "examples"
    examples = sorted(example_dir.glob("*.example.json"))
    if not examples:
        rep.fail("no example documents found under 02-schemas/examples/")
        return

    for example in examples:
        stem = example.name[: -len(".example.json")]
        schema_path = PLANNING / "02-schemas" / f"{stem}.schema.json"
        if not schema_path.exists():
            rep.fail(f"{example.name} has no matching {stem}.schema.json")
            continue
        schema = load_json(schema_path)
        instance = load_json(example)
        validator = Draft202012Validator(schema, registry=registry)
        errors = sorted(validator.iter_errors(instance), key=lambda e: list(e.path))
        if errors:
            for err in errors[:6]:
                loc = "/".join(str(p) for p in err.absolute_path) or "<root>"
                rep.fail(f"{example.name}: at {loc}: {err.message}")
            if len(errors) > 6:
                rep.fail(f"{example.name}: ... and {len(errors) - 6} more errors")
        else:
            rep.ok(f"{example.name} validates against {schema_path.name}")

    # Every schema should have a worked example: an unexercised contract is an
    # unreviewed contract.
    have = {e.name[: -len(".example.json")] for e in examples}
    for path in schema_files():
        stem = path.name[: -len(".schema.json")]
        if stem in ("common",):
            continue
        if stem not in have:
            rep.fail(f"{path.name} has no worked example in 02-schemas/examples/")
        else:
            rep.ok(f"{path.name} has a worked example")


def check_sequence_invariants(rep: Report) -> None:
    """Invariants I1-I4 from sequence.schema.json, checked on the example."""
    rep.section("sequence-invariants")
    codon_table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
        "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
        "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
        "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
        "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
        "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }

    def translate(cds: str) -> str:
        return "".join(codon_table[cds[i:i + 3]] for i in range(0, len(cds), 3))

    example = PLANNING / "02-schemas" / "examples" / "sequence.example.json"
    doc = load_json(example)
    cds = doc["cds"]

    if len(cds) == 3 * doc["codon_count"]:
        rep.ok("I1 len(cds) == 3 * codon_count")
    else:
        rep.fail(f"I1 violated: len(cds)={len(cds)}, codon_count={doc['codon_count']}")

    if translate(cds) == doc["protein"]:
        rep.ok(f"I2 translation matches declared protein ({doc['protein']})")
    else:
        rep.fail(f"I2 violated: translate(cds)={translate(cds)!r} != {doc['protein']!r}")

    for lock in doc.get("locked_codons", []):
        span = lock["span"]
        if 0 <= span["start"] < span["end"] <= doc["codon_count"]:
            rep.ok(f"I4 locked span {span['start']}..{span['end']} in range")
        else:
            rep.fail(f"I4 violated: locked span {span} outside [0,{doc['codon_count']})")

    # The Pareto front example must also carry in-frame, translatable CDSs.
    front = load_json(PLANNING / "02-schemas" / "examples" / "pareto_front.example.json")
    for cand in front["candidates"]:
        for seq in cand["sequences"]:
            s = seq["cds"]
            if len(s) % 3:
                rep.fail(f"{cand['candidate_id']}: cds not a whole number of codons")
                continue
            try:
                prot = translate(s)
            except KeyError as exc:
                rep.fail(f"{cand['candidate_id']}: bad codon {exc}")
                continue
            if prot.count("*") != 1 or not prot.endswith("*"):
                rep.fail(f"{cand['candidate_id']}: internal or missing stop in {prot!r}")
            else:
                rep.ok(f"{cand['candidate_id']}/{seq['cistron_id']} translates to {prot}")


def check_deps(rep: Report) -> None:
    rep.section("deps")
    path = PLANNING / "01-modules" / "modules.yaml"
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    modules = {m["name"]: m for m in doc["modules"]}

    if len(modules) == len(doc["modules"]):
        rep.ok(f"{len(modules)} modules declared, names unique")
    else:
        rep.fail("duplicate module names in modules.yaml")

    runtimes = set(doc["runtimes"])
    for name, mod in modules.items():
        if mod["runtime"] not in runtimes:
            rep.fail(f"{name}: undeclared runtime {mod['runtime']!r}")

    # dangling edges + rank monotonicity
    for name, mod in modules.items():
        for dep in mod["depends_on"]:
            if dep not in modules:
                rep.fail(f"{name} -> {dep}: dependency is not a declared module")
                continue
            if modules[dep]["rank"] >= mod["rank"]:
                rep.fail(
                    f"{name} (rank {mod['rank']}) -> {dep} (rank {modules[dep]['rank']}): "
                    "rank must strictly decrease along an edge"
                )
    rep.ok("all dependency edges reference declared modules with strictly lower rank")

    # acyclicity via Kahn
    indeg = {n: 0 for n in modules}
    adj = defaultdict(list)
    for name, mod in modules.items():
        for dep in mod["depends_on"]:
            if dep in modules:
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
    if seen == len(modules):
        rep.ok(f"dependency graph is a DAG ({seen} modules topologically ordered)")
    else:
        remaining = [n for n, d in indeg.items() if d > 0]
        rep.fail(f"dependency graph contains a cycle among: {sorted(remaining)}")

    # every layer of the architecture has exactly one owning module
    expected_layers = {
        "context-ingestion", "epistatic-surface", "biophysical-predictors",
        "burden-model", "search-engine", "constraint-repair",
        "decision-support", "active-learning",
    }
    present = {m["layer"] for m in modules.values()}
    missing = expected_layers - present
    if missing:
        rep.fail(f"architecture layers with no owning module: {sorted(missing)}")
    else:
        rep.ok("every architectural layer has at least one owning module")

    # public contracts exist
    for name, mod in modules.items():
        contract = ROOT / mod["public_contract"]
        if contract.exists():
            rep.ok(f"{name}: public contract {mod['public_contract']} exists")
        else:
            rep.fail(f"{name}: public contract {mod['public_contract']} does not exist")

    # declared test ids are defined in the test plan
    plan = PLANNING / "05-testing" / "01-test-plan.md"
    plan_text = plan.read_text(encoding="utf-8") if plan.exists() else ""
    for name, mod in modules.items():
        for tid in mod.get("tests", []):
            if tid in plan_text:
                rep.ok(f"{name}: test {tid} defined in the test plan")
            else:
                rep.fail(f"{name}: test {tid} is not defined in 05-testing/01-test-plan.md")


def _collect_ids(pattern: str, text: str) -> set[str]:
    return set(re.findall(pattern, text))


def check_refs(rep: Report) -> None:
    rep.section("refs")
    corpus = {p: p.read_text(encoding="utf-8") for p in all_text_files()}
    joined = "\n".join(corpus.values())

    # ADRs: an ADR-XXXX cited anywhere must exist as a file.
    adr_files = {p.name.split("-")[0] + "-" + p.name.split("-")[1]
                 for p in (PLANNING / "07-decisions").glob("ADR-*.md")}
    for cited in sorted(_collect_ids(r"\bADR-(\d{4})\b", joined)):
        if f"ADR-{cited}" in adr_files:
            rep.ok(f"ADR-{cited} exists")
        else:
            rep.fail(f"ADR-{cited} is cited but no 07-decisions/ADR-{cited}-*.md exists")

    def defined_in(relpath: str, pattern: str) -> set[str]:
        path = PLANNING / relpath
        if not path.exists():
            rep.fail(f"missing definition file: planning/{relpath}")
            return set()
        return _collect_ids(pattern, path.read_text(encoding="utf-8"))

    checks = [
        ("AMD", r"\bAMD-(\d{3})\b", "00-scope/03-architecture-amendments.md"),
        ("EXP", r"\bEXP-(\d{2})\b", "08-experiments/experiment-register.md"),
        ("RISK", r"\bRISK-([A-Z]{3,4}-\d{3})\b", "06-build/06-risk-register.md"),
        ("GATE", r"\bGATE-([A-Z]{3,6}-\d{3})\b", "05-testing/02-acceptance-thresholds.md"),
        ("TST", r"\bTST-([A-Z]{3,5}-\d{3})\b", "05-testing/01-test-plan.md"),
        # BRIEF-4.b: every departure must appear in the departure index.
        ("DEP", r"\bDEP-([A-Z]{2,4}-\d{3})\b", "00-scope/04-departure-index.md"),
        ("SYN", r"\bSYN-(\d{2})\b", "05-testing/03-synthetic-data-suite.md"),
    ]
    for label, pattern, relpath in checks:
        defined = defined_in(relpath, pattern)
        cited = _collect_ids(pattern, joined)
        for ident in sorted(cited - defined):
            rep.fail(f"{label}-{ident} is cited but not defined in planning/{relpath}")
        if cited and not (cited - defined):
            rep.ok(f"all {len(cited)} {label}-* identifiers resolve to planning/{relpath}")

    # Every phase ticket id must live in its own phase file.
    for phase in range(1, 6):
        path = PLANNING / "06-build" / f"0{phase}-phase-{phase}-tickets.md"
        if not path.exists():
            rep.fail(f"missing planning/06-build/0{phase}-phase-{phase}-tickets.md")
            continue
        ids = _collect_ids(rf"\bP{phase}-(\d{{2}})\b", path.read_text(encoding="utf-8"))
        if ids:
            rep.ok(f"phase {phase} defines {len(ids)} tickets")
        else:
            rep.fail(f"phase {phase} ticket file defines no P{phase}-* tickets")


def check_links(rep: Report) -> None:
    rep.section("links")
    link_re = re.compile(r"\[[^\]]*\]\(([^)#]+)(?:#[^)]*)?\)")
    bad = 0
    total = 0
    for md in markdown_files():
        for target in link_re.findall(md.read_text(encoding="utf-8")):
            if target.startswith(("http://", "https://", "mailto:")):
                continue
            total += 1
            resolved = (md.parent / target).resolve()
            if not resolved.exists():
                rep.fail(f"{md.relative_to(ROOT)}: broken link -> {target}")
                bad += 1
    if not bad:
        rep.ok(f"all {total} relative Markdown links resolve")


def check_openapi(rep: Report) -> None:
    rep.section("openapi")
    path = PLANNING / "04-api" / "openapi.yaml"
    if not path.exists():
        rep.fail("planning/04-api/openapi.yaml is missing")
        return
    try:
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:  # noqa: PERF203
        rep.fail(f"openapi.yaml is not valid YAML: {exc}")
        return
    rep.ok(f"openapi.yaml parses (openapi {doc.get('openapi')})")

    if not doc.get("paths"):
        rep.fail("openapi.yaml declares no paths")
    else:
        rep.ok(f"{len(doc['paths'])} paths declared")

    # every external $ref must resolve to a file that exists, and to a JSON
    # pointer inside it where one is given
    def walk(node, refs):
        if isinstance(node, dict):
            for key, value in node.items():
                if key == "$ref" and isinstance(value, str):
                    refs.append(value)
                else:
                    walk(value, refs)
        elif isinstance(node, list):
            for item in node:
                walk(item, refs)

    refs: list[str] = []
    walk(doc, refs)
    external = [r for r in refs if not r.startswith("#")]
    for ref in sorted(set(external)):
        file_part, _, pointer = ref.partition("#")
        target = (path.parent / file_part).resolve()
        if not target.exists():
            rep.fail(f"openapi.yaml: $ref target does not exist: {ref}")
            continue
        if pointer:
            doc_json = load_json(target)
            node = doc_json
            ok = True
            for part in pointer.strip("/").split("/"):
                part = part.replace("~1", "/").replace("~0", "~")
                if isinstance(node, dict) and part in node:
                    node = node[part]
                else:
                    ok = False
                    break
            if not ok:
                rep.fail(f"openapi.yaml: JSON pointer does not resolve: {ref}")
                continue
        rep.ok(f"$ref resolves: {ref}")

    internal = [r for r in refs if r.startswith("#")]
    for ref in sorted(set(internal)):
        node = doc
        ok = True
        for part in ref.lstrip("#").strip("/").split("/"):
            if isinstance(node, dict) and part in node:
                node = node[part]
            else:
                ok = False
                break
        if ok:
            rep.ok(f"$ref resolves: {ref}")
        else:
            rep.fail(f"openapi.yaml: internal $ref does not resolve: {ref}")


def check_trace(rep: Report) -> None:
    rep.section("trace")
    path = PLANNING / "00-scope" / "traceability.yaml"
    if not path.exists():
        rep.fail("planning/00-scope/traceability.yaml is missing")
        return
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))

    seen_ids: set[str] = set()
    for req in doc["requirements"]:
        rid = req["id"]
        if rid in seen_ids:
            rep.fail(f"duplicate requirement id {rid}")
        seen_ids.add(rid)
        artifacts = req.get("artifacts") or []
        if not artifacts:
            rep.fail(f"{rid}: no artifact satisfies this requirement")
            continue
        missing = [rel for rel in artifacts if not (ROOT / rel).exists()]
        for rel in missing:
            rep.fail(f"{rid}: artifact {rel} does not exist")
        if not missing:
            rep.ok(f"{rid}: satisfied by {len(artifacts)} artifact(s)")

    # Every section of the brief must be represented.
    required_prefixes = [
        "BRIEF-2.1", "BRIEF-2.2", "BRIEF-2.3", "BRIEF-2.4", "BRIEF-2.5",
        "BRIEF-2.6", "BRIEF-3.1", "BRIEF-3.2", "BRIEF-3.3", "BRIEF-3.4",
        "BRIEF-3.5", "BRIEF-3.6", "BRIEF-3.7", "BRIEF-4", "BRIEF-5",
    ]
    for prefix in required_prefixes:
        if any(rid.startswith(prefix) for rid in seen_ids):
            rep.ok(f"{prefix}* is covered")
        else:
            rep.fail(f"no requirement covering {prefix}*")


def check_tags(rep: Report) -> None:
    rep.section("tags")
    vocab = load_json(PLANNING / "02-schemas" / "rationale_tags.json")
    known = {t["tag"] for t in vocab["tags"]}
    scopes = {t["tag"]: t["scope"] for t in vocab["tags"]}
    rep.ok(f"{len(known)} rationale tags in the closed vocabulary")

    front = load_json(PLANNING / "02-schemas" / "examples" / "pareto_front.example.json")
    for cand in front["candidates"]:
        for tag in cand["rationale"]["candidate_tags"]:
            if tag not in known:
                rep.fail(f"{cand['candidate_id']}: unknown candidate tag {tag!r}")
            elif scopes[tag] != "candidate":
                rep.fail(f"{cand['candidate_id']}: tag {tag!r} has scope {scopes[tag]}, not candidate")
            else:
                rep.ok(f"{cand['candidate_id']}: candidate tag {tag!r} valid")
        for ann in cand["rationale"]["annotations"]:
            tag = ann["tag"]
            if tag not in known:
                rep.fail(f"{cand['candidate_id']}: unknown annotation tag {tag!r}")
            elif scopes[tag] == "candidate":
                rep.fail(f"{cand['candidate_id']}: candidate-scope tag {tag!r} used as an annotation")
            else:
                rep.ok(f"{cand['candidate_id']}: annotation tag {tag!r} valid")


CHECKS = {
    "schemas": check_schemas,
    "sequence": check_sequence_invariants,
    "deps": check_deps,
    "refs": check_refs,
    "links": check_links,
    "openapi": check_openapi,
    "trace": check_trace,
    "tags": check_tags,
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="append", choices=sorted(CHECKS), default=None)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    rep = Report(args.verbose)
    selected = args.check or list(CHECKS)
    for name in selected:
        CHECKS[name](rep)

    print("\n" + "=" * 70)
    if rep.failures:
        print(f"FAILED: {len(rep.failures)} problem(s) across {rep.assertions} assertions")
        for failure in rep.failures:
            print(f"  - {failure}")
        return 1
    print(f"PASSED: {rep.assertions} assertions across {len(selected)} check group(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
