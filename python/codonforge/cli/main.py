"""The `codonforge` command line.

Ticket P1-30. Specification: `planning/04-api/02-cli-and-library.md` §3.

The CLI builds a `DesignRequest` and calls `codonforge.api.design`. It holds no
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

from codonforge import _core
from codonforge.cli import writers
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator import BUILTIN_TERMS, DesignError, DesignRequest, design

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
    parser = argparse.ArgumentParser(prog="codonforge", description=__doc__)
    parser.add_argument("--version", action="store_true", help="print the engine version")
    sub = parser.add_subparsers(dest="verb")

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--output",
        choices=["table", "json", "fasta", "genbank", "csv", "explain"],
        default="table",
    )
    common.add_argument("--store", default=".codonforge/artifacts")

    d = sub.add_parser("design", parents=[common], help="run a design job")
    d.add_argument("--json", dest="json_file", help="a request document; equals the REST body")
    d.add_argument("--host")
    d.add_argument("--cds")
    d.add_argument("--host-profile", help="a host_profile.schema.json document")
    d.add_argument("--objectives", default="usage.cai")
    d.add_argument("--constraints")
    d.add_argument("--budget", choices=["fast", "standard", "exact"], default="standard")
    d.add_argument("--seed", type=int, default=0)
    d.add_argument("--generations", type=int, default=12)
    d.add_argument("--job-id")
    d.add_argument("--checkpoint-dir", help="DEP-CLI-001: where checkpoints are written")
    d.add_argument("--resume", action="store_true", help="continue from the checkpoint")

    f = sub.add_parser("front", parents=[common], help="show a stored front")
    f.add_argument("--front", required=True)

    s = sub.add_parser("status", parents=[common], help="show a stored job status")
    s.add_argument("--status", required=True)

    feas = sub.add_parser("feasibility", parents=[common], help="analyse a constraint set")
    feas.add_argument("--cds", required=True)
    feas.add_argument("--constraints", required=True)
    feas.add_argument("--translation-table", type=int, default=11)

    sub.add_parser("terms", parents=[common], help="list built-in terms and directions")
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
    if args.verb == "terms":
        return _terms(args)
    if args.verb == "front":
        return _show(Path(args.front), args)
    if args.verb == "status":
        return _show(Path(args.status), args)
    if args.verb == "feasibility":
        return _feasibility(args)
    return _design(args)


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

    _emit(result, args.output)
    return EXIT_OK


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
    return DesignRequest(
        host=args.host,
        cds=args.cds,
        host_profile=_load(args.host_profile),
        objectives=[t for t in args.objectives.split(",") if t],
        constraints=_load(args.constraints, _EMPTY_CONSTRAINTS),
        budget=args.budget,
        seed=args.seed,
        generations=args.generations,
        job_id=args.job_id,
        checkpoint_dir=args.checkpoint_dir,
        resume=args.resume,
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
            f"{term}={value:.4f}"
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
    """A minimal request that constructs one term, to ask its direction."""
    return json.dumps(
        {
            "cds": "ATGTAA",
            "translation_table": 1,
            "five_prime_flank": "",
            "three_prime_flank": "",
            "constraint_set": _EMPTY_CONSTRAINTS,
            "host_profile": _UNIFORM_PROFILE,
            "host_id": "escherichia_coli:511145/BL21",
            "host_digest": "sha256:" + "0" * 64,
            "terms": [term],
            "weights": [1.0],
            "root_seed": 0,
        }
    )


def _feasibility(args: argparse.Namespace) -> int:
    report = json.loads(
        _core.constraint_violations(
            json.dumps(_load(args.constraints, _EMPTY_CONSTRAINTS)),
            args.cds,
            [(0, len(args.cds))],
            "assembled",
        )
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
