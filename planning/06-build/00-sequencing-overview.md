# Build Sequencing — Overview

Satisfies the framing half of `BRIEF-2.6`; the ticket-level breakdown lives in the per-phase files.

| Phase | Architecture roadmap item | Tickets | Exit gates |
| --- | --- | --- | --- |
| 1 | Hard-constraint filter + baseline optimizer | [`01-phase-1-tickets.md`](01-phase-1-tickets.md) (35) | `GATE-DET-001`, Phase-1 test set |
| 2 | Structure + burden Pareto objectives | [`02-phase-2-tickets.md`](02-phase-2-tickets.md) (23) | `GATE-BEN-001/002`, `GATE-BURD-001`, `EXP-04` decided |
| 3 | Epistatic surface | [`03-phase-3-tickets.md`](03-phase-3-tickets.md) (19) | `GATE-EPI-001/002`, `GATE-UNC-001`, `EXP-01`/`EXP-08` decided |
| 4 | Systems-level joint search | [`04-phase-4-tickets.md`](04-phase-4-tickets.md) (12) | `GATE-BEN-003/004`, `GATE-PERF-001` |
| 5 | Active-learning loop | [`05-phase-5-tickets.md`](05-phase-5-tickets.md) (13) | `GATE-ACT-*`, `GATE-REPRO-001`, `EXP-10` running |

Risk register: [`06-risk-register.md`](06-risk-register.md).

---

## Ticket format

Every ticket carries: an id, a title, the owning module, a concrete deliverable, an acceptance criterion expressed as named tests, its dependencies, and a size class.

**Sizes are scope classes, not schedules.** `S` = one module, no new external dependency, no new contract. `M` = one module plus a contract change, or a new dependency. `L` = several modules, or a new algorithm with a correctness argument to establish. `XL` = an algorithm with no reference implementation plus its validation.

---

## Why the phases are ordered this way

The architecture proposal's roadmap order is kept, and detailed planning confirms it for a reason worth recording: **each phase produces the ground truth the next one needs to be validated.**

- Phase 1's constraint layer and exhaustive toy space give Phase 2's multi-objective search a correctness oracle.
- Phase 2's structure and burden objectives, plus the benchmark harness, give Phase 3's epistatic surface something to be measured *against* — without them, "the epistatic surface improves designs" has no denominator.
- Phase 3's per-codon dwell prediction is a direct input to Phase 4's ledger feedback and to the kinetics proxy.
- Phase 4's full objective set is what Phase 5's wet-lab loop recalibrates.

Reordering would mean validating each layer against a moving target.

## Two departures from a naive reading of the roadmap

`DEP-SEQ-001` — **Seed import lands in Phase 1**, not later as a migration convenience. The benchmark harness needs it to place competitor output on CodonForge's axes ([ADR-0007](../07-decisions/ADR-0007-external-tool-seed-import.md) §"Why this earns its place"), so deferring it would force a parallel, unvalidated ingestion path for benchmarking.

`DEP-SEQ-002` — **The batch-yield evaluation protocol and the worker process land in Phase 2**, before any model needs them, exercised by a stub L0 term. The protocol is the load-bearing structural decision of [ADR-0001](../07-decisions/ADR-0001-language-and-runtime-split.md); retrofitting it in Phase 3 alongside the LM would mean changing the search engine's control flow at the same time as introducing a new model, which is the combination most likely to produce an undiagnosable performance regression.

## Cross-phase standing work

Not tickets in any one phase; carried throughout:

- **Schema changes** go through the planning package first: update `planning/02-schemas`, run `validate_planning_package.py`, then propagate to `schemas/v1`. The schemas stay normative rather than becoming documentation of whatever the code does.
- **Every `DEP-*` marker** in this package is re-checked at each phase exit — a departure that stops being justified should be reverted.
- **Every experiment's decision rule is pre-registered** before its data is collected.
