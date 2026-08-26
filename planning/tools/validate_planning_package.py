#!/usr/bin/env python3
"""Consistency checker for the CodonForge planning package.

This is planning-package tooling, not CodonForge backend code. It exists so that
a package made of schemas, graphs and cross-referenced documents can be reviewed
mechanically instead of by eye. It answers, for every check, a question a
reviewer would otherwise have to answer manually:

  schemas       Are the JSON Schemas actually valid draft 2020-12, and do the
                worked examples validate against them?
  digests       Is every digest in the worked examples genuine — a real
                canonicalize-then-hash of its own payload, or a correctly
                threaded reference to one — and is every digest accounted for?
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
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }

    def translate(cds: str) -> str:
        return "".join(codon_table[cds[i : i + 3]] for i in range(0, len(cds), 3))

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
        "context-ingestion",
        "epistatic-surface",
        "biophysical-predictors",
        "burden-model",
        "search-engine",
        "constraint-repair",
        "decision-support",
        "active-learning",
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
    adr_files = {
        p.name.split("-")[0] + "-" + p.name.split("-")[1]
        for p in (PLANNING / "07-decisions").glob("ADR-*.md")
    }
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
        "BRIEF-2.1",
        "BRIEF-2.2",
        "BRIEF-2.3",
        "BRIEF-2.4",
        "BRIEF-2.5",
        "BRIEF-2.6",
        "BRIEF-3.1",
        "BRIEF-3.2",
        "BRIEF-3.3",
        "BRIEF-3.4",
        "BRIEF-3.5",
        "BRIEF-3.6",
        "BRIEF-3.7",
        "BRIEF-4",
        "BRIEF-5",
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
                rep.fail(
                    f"{cand['candidate_id']}: tag {tag!r} has scope {scopes[tag]}, not candidate"
                )
            else:
                rep.ok(f"{cand['candidate_id']}: candidate tag {tag!r} valid")
        for ann in cand["rationale"]["annotations"]:
            tag = ann["tag"]
            if tag not in known:
                rep.fail(f"{cand['candidate_id']}: unknown annotation tag {tag!r}")
            elif scopes[tag] == "candidate":
                rep.fail(
                    f"{cand['candidate_id']}: candidate-scope tag {tag!r} used as an annotation"
                )
            else:
                rep.ok(f"{cand['candidate_id']}: annotation tag {tag!r} valid")


def check_digests(rep: Report) -> None:
    """Every digest in the examples is genuine, threaded, or declared external.

    The worked examples originally carried hand-written placeholders, five of
    which meant two or three unrelated things at once. That made them usable as
    schema-shape fixtures but not as digest fixtures, and it is exactly the
    class of problem that looks fine until something depends on it.

    The load-bearing part of this check is not that the listed digests verify —
    it is that a digest the map does NOT account for is a failure. Without that,
    a newly-added placeholder would simply not be checked.
    """
    rep.section("digests")

    sys.path.insert(0, str(PLANNING / "tools"))
    import jcs  # noqa: PLC0415
    from refresh_example_digests import (  # noqa: PLC0415
        EXAMPLES,
        MAP_PATH,
        pointer_parts,
        resolve,
        resolve_parent,
    )

    spec = yaml.safe_load(MAP_PATH.read_text(encoding="utf-8"))
    try:
        tables = jcs.load_genetic_codes(ROOT)
    except FileNotFoundError as exc:
        rep.fail(str(exc))
        return

    docs = {
        path.name: json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(EXAMPLES.glob("*.example.json"))
    }

    # ---- 1. every occurrence is accounted for -----------------------------
    def walk(node, pointer=""):
        if isinstance(node, dict):
            for key, value in node.items():
                escaped = key.replace("~", "~0").replace("/", "~1")
                yield from walk(value, f"{pointer}/{escaped}")
        elif isinstance(node, list):
            for index, value in enumerate(node):
                yield from walk(value, f"{pointer}/{index}")
        elif isinstance(node, str) and node.startswith("sha256:"):
            yield pointer, node

    declared: dict[tuple[str, str], str] = {}
    for group in ("computed", "reference", "external"):
        for entry in spec.get(group) or []:
            key = (entry["file"], entry["at"])
            if key in declared:
                rep.fail(f"{entry['file']}{entry['at']} is classified twice")
            declared[key] = group

    found = {(name, ptr) for name, doc in docs.items() for ptr, _ in walk(doc)}
    unaccounted = found - set(declared)
    for name, ptr in sorted(unaccounted):
        rep.fail(
            f"{name}{ptr} holds a digest that digest_map.yaml does not classify; "
            "add it as computed, reference or external"
        )
    stale = set(declared) - found
    for name, ptr in sorted(stale):
        rep.fail(f"digest_map.yaml classifies {name}{ptr}, which no longer exists")
    if not unaccounted and not stale:
        rep.ok(f"all {len(found)} digest occurrences across {len(docs)} examples are classified")

    # ---- 2. well-formedness ----------------------------------------------
    malformed = 0
    for name, doc in docs.items():
        for ptr, value in walk(doc):
            hexpart = value.removeprefix("sha256:")
            if len(hexpart) != 64 or not all(c in "0123456789abcdef" for c in hexpart):
                rep.fail(f"{name}{ptr}: {value!r} is not sha256 + 64 lower-case hex digits")
                malformed += 1
    if not malformed:
        rep.ok("every digest is well-formed")

    # ---- 3. computed digests actually verify ------------------------------
    for entry in spec["computed"]:
        doc = docs.get(entry["file"])
        if doc is None:
            rep.fail(f"digest_map.yaml references missing example {entry['file']}")
            continue
        rule = entry["rule"]
        try:
            if rule == "self_document":
                parent, _ = resolve_parent(doc, entry["at"])
                expected = jcs.digest_of_document(parent, "digest")
                what = "self"
            elif rule == "codon_sequence_document":
                cds = resolve(doc, entry["cds"])
                table = (
                    entry["table_const"] if "table_const" in entry else resolve(doc, entry["table"])
                )
                expected = jcs.codon_sequence_document(cds, int(table), tables)["digest"]
                what = "codon sequence"
            elif rule == "raw_bytes":
                expected = jcs.digest_of_bytes(resolve(doc, entry["source"]).encode("utf-8"))
                what = "raw bytes"
            else:
                rep.fail(f"{entry['file']}{entry['at']}: unknown rule {rule!r}")
                continue
            actual = resolve(doc, entry["at"])
        except (KeyError, IndexError, ValueError) as exc:
            rep.fail(f"{entry['file']}{entry['at']}: cannot evaluate ({exc})")
            continue
        if actual == expected:
            rep.ok(f"{entry['file']}{entry['at']} ({what}) matches a fresh hash of its payload")
        else:
            rep.fail(
                f"{entry['file']}{entry['at']} ({what}) is stale: "
                f"stored {actual[:22]}…, recomputed {expected[:22]}…"
            )

    # ---- 4. references point at the right value ---------------------------
    for entry in spec["reference"]:
        target_file, _, target_pointer = entry["same_as"].partition("#")
        try:
            expected = resolve(docs[target_file], target_pointer)
            actual = resolve(docs[entry["file"]], entry["at"])
        except (KeyError, IndexError, ValueError) as exc:
            rep.fail(f"{entry['file']}{entry['at']}: cannot resolve reference ({exc})")
            continue
        if actual == expected:
            rep.ok(f"{entry['file']}{entry['at']} -> {entry['same_as']}")
        else:
            rep.fail(
                f"{entry['file']}{entry['at']} should equal {entry['same_as']} "
                f"({expected[:22]}…) but holds {actual[:22]}…"
            )

    # ---- 5. externals are declared with a reason --------------------------
    for entry in spec.get("external") or []:
        if not entry.get("reason"):
            rep.fail(f"{entry['file']}{entry['at']}: external digests must state a reason")
        else:
            rep.ok(f"{entry['file']}{entry['at']} declared external: {entry['reason'][:60]}")

    # ---- 6. same dataset_id, same digest ----------------------------------
    if (spec.get("consistency") or {}).get("same_dataset_id_same_digest"):
        by_dataset: dict[str, set[str]] = defaultdict(set)
        where: dict[str, list[str]] = defaultdict(list)

        def collect(node, name, pointer=""):
            if isinstance(node, dict):
                if isinstance(node.get("dataset_id"), str) and isinstance(node.get("digest"), str):
                    by_dataset[node["dataset_id"]].add(node["digest"])
                    where[node["dataset_id"]].append(f"{name}{pointer}")
                for key, value in node.items():
                    collect(value, name, f"{pointer}/{key}")
            elif isinstance(node, list):
                for index, value in enumerate(node):
                    collect(value, name, f"{pointer}/{index}")

        for name, doc in docs.items():
            collect(doc, name)
        for dataset_id, digests in sorted(by_dataset.items()):
            if len(digests) > 1:
                rep.fail(
                    f"dataset {dataset_id} carries {len(digests)} different digests at "
                    f"{', '.join(where[dataset_id])}"
                )
            else:
                rep.ok(
                    f"dataset {dataset_id}: one digest across {len(where[dataset_id])} reference(s)"
                )

    # ---- 7. no value means two different things ---------------------------
    # The original failure: one placeholder serving as the self digest of two
    # unrelated documents. Two self digests colliding now would mean two
    # documents with identical payloads, which is worth knowing about.
    self_entries = [e for e in spec["computed"] if e["rule"] == "self_document"]
    seen: dict[str, str] = {}
    collision = False
    for entry in self_entries:
        try:
            value = resolve(docs[entry["file"]], entry["at"])
        except (KeyError, IndexError):
            continue
        origin = f"{entry['file']}{entry['at']}"
        if value in seen and pointer_parts(entry["at"]) == ["digest"]:
            rep.fail(f"self digest collision: {origin} and {seen[value]} share a value")
            collision = True
        seen.setdefault(value, origin)
    if not collision:
        rep.ok(f"{len(self_entries)} self digests are distinct where they must be")


CHECKS = {
    "schemas": check_schemas,
    "digests": check_digests,
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
