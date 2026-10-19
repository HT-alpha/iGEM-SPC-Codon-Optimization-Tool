"""The `dna-design-suite` command line.

Ticket P1-30. Specification: `planning/04-api/02-cli-and-library.md` §3.

The CLI builds a `DesignRequest` and calls `dna_design_suite.api.design`. It holds no
validation, no defaults and no objective presets of its own — every one of
those would be a second code path, and a second code path is how three surfaces
that are supposed to agree stop agreeing.

## Exit codes

From §3.1's table, which exists so a script can tell an infeasible constraint
set from a crash:

    0  success            3  infeasible constraints
    1  design failed      4  quota/budget
    2  invalid input      5  internal

## Cancellation

`DEP-CLI-002`. `SIGINT` sets `P1-29`'s cooperative cancellation flag, so the job
stops at the next segment boundary with a checkpoint written and exits `0` in
the `cancelled` state. A second `SIGINT` aborts immediately: a user pressing it
twice has stopped asking politely, and the checkpoint from the last completed
segment is already on disk, so nothing is lost that was not already unwritten.
"""

from __future__ import annotations

import argparse
import json
import signal
import sys
from pathlib import Path
from typing import Any

from dna_design_suite import _core
from dna_design_suite.cli import writers
from dna_design_suite.codon.orchestrator import BUILTIN_TERMS, DesignError, DesignRequest, design
from dna_design_suite.modelreg import ArtifactStore

EXIT_OK = 0
EXIT_DESIGN_FAILED = 1
EXIT_INVALID_INPUT = 2
EXIT_INFEASIBLE = 3
EXIT_BUDGET = 4
EXIT_INTERNAL = 5

#: Codes that mean "your constraints cannot be satisfied", not "we broke".
_INFEASIBLE_CODES = {"CF-FEAS-001", "CF-FEAS-002", "CF-FEAS-003", "CF-SRCH-006"}
_BUDGET_CODES = {"CF-ORCH-002", "CF-BUDGET-001", "CF-BUDGET-002"}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="dna-design-suite", description=__doc__)
    parser.add_argument("--version", action="store_true", help="print the engine version")
    sub = parser.add_subparsers(dest="verb")

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--output",
        choices=["table", "json", "fasta", "genbank", "csv", "explain", "sbol"],
        default="table",
    )
    common.add_argument("--store", default=".dna-design-suite/artifacts")

    d = sub.add_parser("design", parents=[common], help="run a design job")
    d.add_argument("--json", dest="json_file", help="a request document; equals the REST body")
    d.add_argument("--host")
    d.add_argument("--cds")
    d.add_argument("--host-profile", help="a host_profile.schema.json document")
    d.add_argument(
        "--source-profile",
        help="codon usage of the organism the CDS comes from, for usage.harmonize_rca",
    )
    d.add_argument("--objectives", default="usage.cai")
    d.add_argument("--constraints")
    d.add_argument("--budget", choices=["fast", "standard", "exact"], default="standard")
    d.add_argument("--seed", type=int, default=0)
    d.add_argument("--generations", type=int, default=12)
    d.add_argument("--job-id")
    d.add_argument("--checkpoint-dir", help="DEP-CLI-001: where checkpoints are written")
    d.add_argument("--resume", action="store_true", help="continue from the checkpoint")
    d.add_argument("--starting-cds", help="synonymous CDS to search from")
    d.add_argument(
        "--structure-engine",
        choices=["viennarna", "nussinov"],
        help="viennarna, or the opt-in approximate nussinov matching",
    )

    f = sub.add_parser("front", parents=[common], help="show a stored front")
    f.add_argument("--front", required=True)

    s = sub.add_parser("status", parents=[common], help="show a stored job status")
    s.add_argument("--status", required=True)

    feas = sub.add_parser("feasibility", parents=[common], help="analyse a constraint set")
    feas.add_argument("--cds", required=True)
    feas.add_argument("--constraints", required=True)
    feas.add_argument("--translation-table", type=int, default=11)

    sub.add_parser("terms", parents=[common], help="list built-in terms and directions")

    seeds = sub.add_parser("seeds", help="import a starting CDS")
    seeds_sub = seeds.add_subparsers(dest="seeds_verb")
    seed_import = seeds_sub.add_parser(
        "import", parents=[common], help="check a CDS against a protein"
    )
    seed_import.add_argument("--tool", default="manual")
    seed_import.add_argument("--input", required=True, help="tool output or a raw CDS")
    seed_import.add_argument("--target-protein")
    seed_import.add_argument("--target-cds")
    seed_import.add_argument("--translation-table", type=int, default=11)
    seed_import.add_argument("--host")

    hosts = sub.add_parser("hosts", help="import a codon-usage table")
    hosts_sub = hosts.add_subparsers(dest="hosts_verb")
    host_import = hosts_sub.add_parser(
        "import",
        parents=[common],
        help="read Kazusa text, a delimited or JSON table, or CDS",
    )
    host_import.add_argument("--input", required=True)
    host_import.add_argument(
        "--format",
        default="auto",
        choices=["auto", "kazusa", "table", "cds", "json"],
    )
    host_import.add_argument("--translation-table", type=int, default=11)
    host_import.add_argument("--organism", default="")
    host_import.add_argument(
        "--label",
        default="user_supplied",
        choices=["user_supplied", "highly_expressed"],
    )

    screen = sub.add_parser("vendor-screen", help="dry-run synthesis-rule screen")
    screen.add_argument("--cds", required=True)
    screen.add_argument(
        "--set",
        dest="sets",
        action="append",
        choices=["dds_synthesis", "dds_assembly"],
        help="repeat to choose a set; default is both",
    )

    fold = sub.add_parser("fold", help="labelled global fold, separate from windowed accessibility")
    fold.add_argument("--cds", required=True)
    fold.add_argument("--temperature", type=float, default=37.0)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(_core.version())
        return EXIT_OK
    if args.verb is None:
        parser.print_help()
        return EXIT_INVALID_INPUT

    try:
        return _dispatch(args)
    except DesignError as error:
        return _report(error.code, error.detail, _exit_for(error))
    except FileNotFoundError as error:
        return _report("CF-CLI-001", str(error), EXIT_INVALID_INPUT)
    except ValueError as error:
        return _report(*_from_value_error(error))


def _dispatch(args: argparse.Namespace) -> int:
    if getattr(args, "output", None) == "sbol" and args.verb != "design":
        raise DesignError(
            "CF-SBOL-002",
            "sbol export is a design result: the recommended candidate, the sequence "
            "document, the host profile, and the constraint set. "
            "Run design --output sbol. This command does not carry those documents.",
            exit_code=EXIT_INVALID_INPUT,
        )
    if args.verb == "terms":
        return _terms(args)
    if args.verb == "front":
        return _show(Path(args.front), args)
    if args.verb == "status":
        return _show(Path(args.status), args)
    if args.verb == "feasibility":
        return _feasibility(args)
    if args.verb == "seeds":
        return _seeds(args)
    if args.verb == "hosts":
        return _hosts(args)
    if args.verb == "vendor-screen":
        return _vendor_screen(args)
    if args.verb == "fold":
        return _fold(args)
    return _design(args)


def _seeds(args: argparse.Namespace) -> int:
    if getattr(args, "seeds_verb", None) != "import":
        raise DesignError(
            "CF-CLI-001",
            "use: seeds import --input <file>",
            exit_code=EXIT_INVALID_INPUT,
        )
    from dna_design_suite.ingest.seeds import KNOWN_TOOLS, import_seeds

    tool = args.tool or "manual"
    if tool not in KNOWN_TOOLS:
        raise DesignError(
            "CF-CLI-001",
            f"tool {tool!r} is not one of {', '.join(sorted(KNOWN_TOOLS))}",
            exit_code=EXIT_INVALID_INPUT,
        )
    payload = Path(args.input).read_text(encoding="utf-8")
    protein = args.target_protein
    if not protein:
        if not args.target_cds:
            raise DesignError(
                "CF-CLI-001",
                "seeds import needs --target-protein or --target-cds",
                exit_code=EXIT_INVALID_INPUT,
            )
        protein = json.loads(
            _core.sequence_document(args.target_cds, args.translation_table, None, None, None)
        )["protein"]
    document = import_seeds(
        protein,
        args.translation_table,
        [{"tool": tool, "payload": payload}],
        host_id=args.host,
    )
    print(json.dumps(document, indent=2, sort_keys=True))
    accepted = [item for item in document["seeds"] if item["acceptance"]["accepted"]]
    return EXIT_OK if accepted else EXIT_INVALID_INPUT


def _hosts(args: argparse.Namespace) -> int:
    if getattr(args, "hosts_verb", None) != "import":
        raise DesignError(
            "CF-CLI-001",
            "use: hosts import --input <file>",
            exit_code=EXIT_INVALID_INPUT,
        )
    from dna_design_suite.ingest.errors import IngestError
    from dna_design_suite.ingest.usage_table import parse_usage_table

    text = Path(args.input).read_text(encoding="utf-8")
    try:
        profile = parse_usage_table(
            text,
            translation_table=args.translation_table,
            source_format=args.format,
            organism=args.organism,
            label=args.label,
        )
    except IngestError as exc:
        raise DesignError(exc.code, exc.detail, exit_code=EXIT_INVALID_INPUT) from exc
    print(json.dumps(profile, indent=2, sort_keys=True))
    return EXIT_OK


def _design(args: argparse.Namespace) -> int:
    request = _request(args)
    store = ArtifactStore(Path(args.store))
    interrupts = {"count": 0}

    def on_interrupt(_signum: int, _frame: Any) -> None:
        interrupts["count"] += 1
        if interrupts["count"] > 1:
            raise KeyboardInterrupt
        print(
            "cancellation requested; stopping at the next segment boundary "
            "(interrupt again to abort now)",
            file=sys.stderr,
        )

    previous = signal.getsignal(signal.SIGINT)
    signal.signal(signal.SIGINT, on_interrupt)
    try:
        result = design(
            request,
            store,
            on_segment=lambda job: setattr(job, "cancel_requested", bool(interrupts["count"])),
        )
    finally:
        signal.signal(signal.SIGINT, previous)

    if result.job.state == "failed":
        return _failed(result, args.output)
    _emit(result, args.output)
    return EXIT_OK


def _failed(result: Any, output: str) -> int:
    """A failed design exits non-zero and does not print the input as a CDS.

    The status document still says which tiers ran. The problem detail is the
    thing a script branches on, so it goes to stderr as `code: detail`.
    """
    status = result.status
    problem = status.get("problem") or {}
    code = str(problem.get("code") or "CF-CLI-002")
    detail = str(problem.get("detail") or "design failed")
    if output == "json":
        print(writers.to_json({"status": status, "front": None}))
    else:
        _render_status(result)
        print(f"problem  {code}")
        print(f"    {detail}")
    return _report(code, detail, _exit_for_code(code))


def _render_status(result: Any) -> None:
    """Tiers and warnings, without a candidate sequence."""
    status = result.status
    print(f"job      {status['job_id']}  {status['state']}")
    print(f"lockfile {status['landscape_lockfile_digest']}")
    print("tiers")
    for tier in status["progress"]["tiers"]:
        detail = ""
        if tier.get("generation") is not None:
            detail = f"  {tier['generation']}/{tier.get('max_generations')}"
        print(f"  {tier['tier']:<20} {tier['state']}{detail}")
    if status["warnings"]:
        print("warnings")
        for warning in status["warnings"]:
            print(f"  {warning['code']}  {warning['stage']}/{warning['subject']}")
            print(f"      {warning['title']}")


def _exit_for_code(code: str) -> int:
    if code in _INFEASIBLE_CODES:
        return EXIT_INFEASIBLE
    if code in _BUDGET_CODES:
        return EXIT_BUDGET
    return EXIT_DESIGN_FAILED


def _emit(result: Any, output: str) -> None:
    """Render a finished design in the requested format.

    Every branch below hands the *documents* to a writer. None of them reaches
    into the job for a direction, a score or a tier state, because a value
    computed twice is a value that can disagree with itself.
    """
    if output == "json":
        print(writers.to_json({"status": result.status, "front": result.front}))
    elif output == "fasta":
        print(writers.to_fasta(result.front), end="")
    elif output == "genbank":
        print(writers.to_genbank(result.front, result.sequence_document()), end="")
    elif output == "csv":
        print(writers.to_csv(result.front), end="")
    elif output == "explain":
        print(writers.explain(result.front, result.status), end="")
    elif output == "sbol":
        print(result.to_sbol(), end="")
    else:
        _render(result)


def _request(args: argparse.Namespace) -> DesignRequest:
    if args.json_file:
        return DesignRequest.from_json(Path(args.json_file).read_text(encoding="utf-8"))
    missing = [flag for flag, value in (("--host", args.host), ("--cds", args.cds)) if not value]
    if missing:
        raise DesignError(
            "CF-CLI-001",
            f"missing required {' and '.join(missing)} (or pass --json)",
            exit_code=EXIT_INVALID_INPUT,
        )
    profile = _load(args.host_profile)
    if args.source_profile:
        source = _load(args.source_profile)
        block = source.get("codon_usage") if isinstance(source.get("codon_usage"), dict) else source
        profile = dict(profile)
        profile["source_codon_usage"] = block
    return DesignRequest(
        host=args.host,
        cds=args.cds,
        host_profile=profile,
        objectives=[t for t in args.objectives.split(",") if t],
        constraints=_load(args.constraints, _EMPTY_CONSTRAINTS),
        budget=args.budget,
        seed=args.seed,
        generations=args.generations,
        job_id=args.job_id,
        checkpoint_dir=args.checkpoint_dir,
        resume=args.resume,
        starting_cds=args.starting_cds or "",
        structure_engine=args.structure_engine or "",
    )


def _render(result: Any) -> None:
    """The human view. Per-tier, because progress is per-tier.

    Not flattened to a percentage: a job that ran tier 1 and skipped tier 2 and
    one that ran both are different outcomes, and a single number cannot say
    which happened.
    """
    status = result.status
    print(f"job      {status['job_id']}  {status['state']}")
    print(f"lockfile {status['landscape_lockfile_digest']}")
    print("objectives")
    for objective in status["resolved_objectives"]:
        print(f"  {objective['term_id']:<32} {objective['direction']}")
    print("tiers")
    for tier in status["progress"]["tiers"]:
        detail = ""
        if tier.get("generation") is not None:
            detail = f"  {tier['generation']}/{tier.get('max_generations')}"
        print(f"  {tier['tier']:<20} {tier['state']}{detail}")
    if status["warnings"]:
        print("warnings")
        for warning in status["warnings"]:
            print(f"  {warning['code']}  {warning['stage']}/{warning['subject']}")
            print(f"      {warning['title']}")
    front = result.front
    reduction = front["reduction"]["method"]
    print(f"front    {len(front['candidates'])} candidate(s), reduction={reduction}")
    for cand in front["candidates"]:
        values = ", ".join(
            # UNKNOWN is a value, not a formatting problem. An out-of-core term
            # the budget never reached has no number, and printing a placeholder
            # would make an unscored axis look measured.
            f"{term}={'unknown' if value is None else format(value, '.4f')}"
            for term, value in zip(front["objective_order"], cand["objective_vector"], strict=True)
        )
        print(f"  {cand['candidate_id']}  {values}")
        print(f"    {cand['sequences'][0]['cds']}")


def _terms(args: argparse.Namespace) -> int:
    rows = [(term, _core.term_directions(_probe(term))) for term in BUILTIN_TERMS]
    parsed = [(term, json.loads(text)[term]) for term, text in rows]
    if args.output == "json":
        print(json.dumps(dict(parsed), indent=2, sort_keys=True))
    else:
        for term, direction in parsed:
            print(f"{term:<32} {direction}")
    return EXIT_OK


def _probe(term: str) -> str:
    """A minimal request that constructs one term, to ask its direction.

    `usage.harmonize_rca` compares the input CDS with a source table. The
    probe uses the same uniform table as both, and pins `ATGTAA` as the
    input. That is enough to read the direction. It is not a design.
    """
    profile = dict(_UNIFORM_PROFILE)
    cds = "ATGTAA"
    if term == "usage.harmonize_rca":
        profile["source_codon_usage"] = profile["codon_usage"]
        profile["harmonize_native_cds"] = cds
    return json.dumps(
        {
            "cds": cds,
            "translation_table": 1,
            "five_prime_flank": "",
            "three_prime_flank": "",
            "constraint_set": _EMPTY_CONSTRAINTS,
            "host_profile": profile,
            "host_id": "escherichia_coli:511145/BL21",
            "host_digest": "sha256:" + "0" * 64,
            "terms": [term],
            "weights": [1.0],
            "root_seed": 0,
        }
    )


def _feasibility(args: argparse.Namespace) -> int:
    from dna_design_suite.codon.orchestrator import local_structure

    constraints = _load(args.constraints, _EMPTY_CONSTRAINTS)
    try:
        local_structure.validate(constraints)
        local_structure.require_engine(constraints, False)
    except DesignError as error:
        print(f"{error.code}: {error.detail}", file=sys.stderr)
        return error.exit_code
    report = json.loads(
        _core.constraint_violations(
            json.dumps(constraints),
            args.cds,
            [(0, len(args.cds))],
            "assembled",
        )
    )
    evaluated = local_structure.evaluated_ids(constraints)
    report["undetectable"] = [item for item in report["undetectable"] if item not in evaluated]
    cache: dict[tuple[Any, ...], float] = {}
    for item in local_structure.rules(constraints):
        if not item["blocking"]:
            continue
        value = local_structure.measure(args.cds, item["rule"], cache)
        if local_structure.holds(value, item["rule"]):
            continue
        span = item["rule"]["span"]
        report["violations"].append(
            {
                "constraint_id": item["id"],
                "class": "blocking",
                "start_nt": int(span["start_nt"]),
                "end_nt": int(span["end_nt"]),
                "detail": f"{item['rule']['metric']} is {value}",
                "actionable": True,
            }
        )
    violations = report["violations"]
    if args.output == "json":
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        for violation in violations:
            mark = "actionable" if violation["actionable"] else "locked"
            print(f"{violation['constraint_id']}  [{mark}]  {violation['detail']}")
        for note in report["undetectable"]:
            print(f"undetectable: {note}")
    # Exit 3 only for violations the design could actually have avoided. A
    # violation confined to immutable sequence is a property of the input, not
    # a failure of the constraint set, and DEP-CONS-002 draws that line
    # everywhere else too.
    return EXIT_INFEASIBLE if any(v["actionable"] for v in violations) else EXIT_OK


def _show(path: Path, args: argparse.Namespace) -> int:
    document = json.loads(path.read_text(encoding="utf-8"))
    if args.output in ("fasta", "csv") and "candidates" in document:
        print(
            writers.to_fasta(document) if args.output == "fasta" else writers.to_csv(document),
            end="",
        )
    else:
        print(writers.to_json(document))
    return EXIT_OK


def _vendor_screen(args: argparse.Namespace) -> int:
    from dna_design_suite.vendor.screen import screen_sequence

    try:
        report = screen_sequence(args.cds, sets=args.sets)
    except ValueError as exc:
        raise DesignError("CF-CLI-001", str(exc), exit_code=EXIT_INVALID_INPUT) from exc
    print(json.dumps(report, indent=2, sort_keys=True))
    return EXIT_OK


def _fold(args: argparse.Namespace) -> int:
    from dna_design_suite.codon.predictors.global_fold import global_fold_report

    try:
        report = global_fold_report(args.cds, temperature_c=args.temperature)
    except ValueError as exc:
        raise DesignError("CF-CLI-001", str(exc), exit_code=EXIT_INVALID_INPUT) from exc
    print(json.dumps(report, indent=2, sort_keys=True))
    return EXIT_OK


def _load(path: str | None, empty: dict[str, Any] | None = None) -> dict[str, Any]:
    if not path:
        return dict(empty) if empty else {}
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _exit_for(error: DesignError) -> int:
    if error.code in _INFEASIBLE_CODES:
        return EXIT_INFEASIBLE
    if error.code in _BUDGET_CODES:
        return EXIT_BUDGET
    return error.exit_code


def _from_value_error(error: ValueError) -> tuple[str, str, int]:
    """Map a Rust error, which arrives as a JSON problem document."""
    try:
        problem = json.loads(str(error))
        code = problem["code"]
    except (ValueError, KeyError):
        return ("CF-CLI-001", str(error), EXIT_INVALID_INPUT)
    if code in _INFEASIBLE_CODES:
        return (code, problem.get("detail", ""), EXIT_INFEASIBLE)
    if code in _BUDGET_CODES:
        return (code, problem.get("detail", ""), EXIT_BUDGET)
    return (code, problem.get("detail", ""), EXIT_INVALID_INPUT)


def _report(code: str, detail: str, exit_code: int) -> int:
    print(f"{code}: {detail}", file=sys.stderr)
    return exit_code


#: A valid, empty constraint set. `ConstraintSet` requires an id even when it
#: constrains nothing, so "no constraints" still has to be spelled out.
_EMPTY_CONSTRAINTS = {
    "schema_version": "1.0.0",
    "constraint_set_id": "cst_none",
    "constraints": [],
}

#: A flat table, used only to construct a term so it can be asked its direction.
#: Never used to score anything — a uniform CAI table is meaningless as a
#: reference set, and `CodonUsageTable` refuses one for scoring.
_UNIFORM_PROFILE = {
    "codon_usage": {
        "subset": "highly_expressed",
        "table": {
            a + b + c: {"relative_adaptiveness": 0.5}
            for a in "TCAG"
            for b in "TCAG"
            for c in "TCAG"
        },
    }
}


if __name__ == "__main__":
    raise SystemExit(main())
