# Phase 2 — Structure + Burden Pareto Objectives

**What Phase 2 delivers:** the multi-objective machinery, the two biophysical objectives that do not need learned models, the burden model, the decision-support layer, and the benchmark harness. At the end of Phase 2 CodonForge is a credible multi-objective codon optimizer that can be compared head-to-head against LinearDesign and CDSfold — which is what makes every later claim measurable.

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P2-01` | Indicator math | `cf-metrics` | Fast non-dominated sort, crowding distance, Das–Dennis (single and two-layer), WFG hypervolume, IGD+, ε-indicator | `TST-METR-001`, `TST-METR-002`, `TST-METR-003` | Phase 1 | M |
| `P2-02` | NSGA-II loop | `cf-search` | The `M ≤ 3` path, with constraint-domination and the `UNKNOWN` rule | `TST-SRCH-007` | `P2-01` | M |
| `P2-03` | NSGA-III loop | `cf-search` | Normalization (ideal point, extreme points, intercepts, degeneracy guards), reference-direction association, niching selection; automatic `M ≥ 4` switch | `TST-SRCH-007`, `TST-METR-002` | `P2-02` | L |
| `P2-04` | Termination and progress | `cf-search` | Hypervolume-stagnation termination; `hypervolume_history`; per-tier progress reporting | `contract_*` against `job_status.schema.json` | `P2-03` | S |
| `P2-05` | Batch-yield protocol | `cf-search`, `cf-landscape` | The resumable `step`/`absorb`/`snapshot`/`resume` state machine; Arrow IPC batch transport; `cf-worker` process; exercised end-to-end by a **stub** L0 term | `TST-LM-001` applied to the stub; one boundary crossing per generation asserted | `P2-03` | L |
| `P2-06` | Surrogate rationing | `cf-search` | Ridge surrogate on k-mer + L1-term features; the four-way `chosen` set including the mandatory random 10 %; periodic refit | `TST-SRCH-005` | `P2-05` | M |
| `P2-07` | Front certification | `cf-search` | Exact re-evaluation of the reduced front before emission; warning when certification changes the front | `TST-SRCH-005` | `P2-06` | S |
| `P2-08` | Native structure kernel | `cf-biophys-structure` | Permissively-licensed windowed accessibility (`O(n·w²)`), temperature-parameterized | `TST-STRU-002` | Phase 1 | XL |
| `P2-09` | Structure engine seam | `cf-biophys-structure` | `StructureEngine` trait; runtime detection of user-installed ViennaRNA / LinearFold; **no vendoring**; engine recorded in the lockfile | `TST-STRU-001`, `TST-STRU-003` | `P2-08` | M |
| `P2-10` | Incremental structure cache | `cf-biophys-structure` | Window-keyed cache; `O(w³)` recomputation per 3-nt edit; exactness against full recomputation | `TST-STRU-002` | `P2-08` | L |
| `P2-11` | Structure terms | `cf-biophys-structure` | `structure.local_accessibility` (L1) and `structure.mfe_global` (L0, `exact` tier only); `local_structure` constraint support | `contract_*`; delta-correctness test | `P2-10` | M |
| `P2-12` | Wobble table and ledger | `cf-burden` | `crick_wobble_v1` decoding table with inosine; demand accumulation; native-proteome demand; `O(1)` incremental update | `TST-BURD-001`, `TST-BURD-002` | Phase 1 | M |
| `P2-13` | Capacity estimators | `cf-burden` | Tiers 1, 2 and 3, with tier selection from `data_completeness`; tier-3 semantics enforcement (unit, `calibrated: false`, tag) | `TST-BURD-005` | `P2-12` | L |
| `P2-14` | Burden penalty | `cf-burden` | Bounded saturating `g` for reporting; unbounded soft hinge `H` for search; soft-max aggregation; frozen-availability delta | `TST-BURD-003`, `TST-BURD-004` | `P2-13` | M |
| `P2-15` | Front reduction | `cf-decision` | Reference-direction representatives, knee selection, clustering, ε-grid; extreme-point retention | `TST-DEC-003` | `P2-01` | M |
| `P2-16` | Rationale tagger | `cf-decision` | The closed vocabulary, trigger evaluation from `per_position` tracks, template rendering, `diff_vs_input` | `TST-DEC-002` | `P2-15` | M |
| `P2-17` | Kinetics v1 proxy | `cf-biophys-kinetics` | Domain-map ingestion (annotation → predicted → refuse); desired-window geometry; concordance metric; wired as a **`diagnostic`** term only | `TST-KIN-001`, `TST-KIN-002` | `P2-11` | L |
| `P2-18` | REST API | `cf-api` | Every endpoint in the OpenAPI document; `problem+json`; idempotency; SSE progress; tenancy middleware (single-tenant mode until Phase 5) | `TST-API-001`, `TST-API-002`, `TST-API-004` | `P2-04` | L |
| `P2-19` | Benchmark harness | `cf-bench` | Dataset manifests; baseline adapters; `metrics.parquet`; `gates.json`; `report.html`; comparison refusal rules | `TST-BEN-001`, `TST-BEN-002` | `P2-18` | L |
| `P2-20` | `bench_mfe_cai_frontier` | `cf-bench` | Two-objective head-to-head with a competitor `λ` sweep, shared folding engine, feasibility reported alongside | `GATE-BEN-001`, `GATE-BEN-002` | `P2-19` | M |
| `P2-21` | `bench_burden_backtest` + `EXP-04` | `cf-bench` | Burden back-test at every capacity tier, with the tier-3 path run with abundance data withheld | `GATE-BURD-001`; `EXP-04` decided | `P2-19` | M |
| `P2-22` | `SYN-03`, `SYN-06` | `cf-testkit` | Two-cistron contention instance with closed-form margin; three-tier synthetic host | `TST-SYN-003` | `P2-12` | M |
| `P2-23` | Licence and patent review | — | Written review of every dependency's licence and of the LinearDesign patent position; confirmation that no non-redistributable code is vendored; `NOTICE` file | Signed-off review document in `docs/` | `P2-09` | S |

### `P2-01` / `P2-02` — the audit's corrections

**Crowding distance had no acceptance test.** `P2-01` lists it as a deliverable and no `TST-*` covered it, which matters more than an ordinary omission because crowding distance is where the tie behaviour lives: boundary members take infinity, and a front whose members share an objective value has a zero range that a naive normalization divides by. `TST-METR-004` now covers both, plus a fully degenerate front, a single-member front and an unknown term. Writing it found that a one-member front took distance zero and could be discarded by selection.

**`P2-02`'s only acceptance criterion was a comparator truth table.** A loop that ranks correctly, loses its population on resume, and reports one member of a tied front would have passed `TST-SRCH-007`. `TST-SRCH-016` (checkpoint completeness across five splits of the same total) and `TST-SRCH-017` (a tied front reported whole, and independent of member order) close that.

**The reported front is not the scalar argmax of a front.** `CF-ORCH-007`'s guarantee is extended rather than reduced: `CandidateArchive::front` returns the non-dominated set and `CandidateArchive::best` returns the scalarization's winner, and neither is derived from the other. Collapsing a trade-off curve to one point hands a caller who asked for three objectives a single design, with nothing in the document saying the rest were found and discarded. Feasibility partitions the front the same way it outranks score in `best`, and an `UNKNOWN` term makes a candidate incomparable rather than dominated away.

**Objective vectors travel in the checkpoint, not just the genotypes.** Recomputing them on resume would be correct only if evaluation were pure and cheap, and it is neither: an `UNKNOWN` is a fact about a budget already spent, and re-running could turn it into a value and change the ranking.

**Absent, not stubbed.** §7's loop calls a surrogate (`P2-06`), reference-direction niching (`P2-03`) and certification (`P2-07`). `P2-02` evaluates exactly, selects on crowding distance, emits an uncertified front, and refuses `M > 3` with `CF-SRCH-014` naming `P2-03` — rather than approximating any of the three, since a stub returning a plausible number is indistinguishable from a working component until something depends on it.

`P2-23` is a ticket rather than a background task because `AMD-009` identifies it as a blocking risk for the open-source goal, and because the cost of discovering it at release rather than at Phase 2 is a redesign of the structure layer.

## Phase 2 exit checklist

- [ ] Phase 1 checklist still green
- [ ] `EXIT-003` decided — `EXP-05` (Lamarckian vs Baldwinian repair), moved here from Phase 1 because its rule turns on hypervolume and population diversity and Phase 1 had no population. `DEP-EXP-001`.
- [ ] `TST-METR-*`, `TST-STRU-*`, `TST-BURD-*`, `TST-DEC-*`, `TST-KIN-*`, `TST-API-*`, `TST-BEN-*` green
- [ ] `GATE-BEN-001`, `GATE-BEN-002`, `GATE-BURD-001` pass
- [ ] `EXP-04` decided; the tier-3 burden term's role is settled
- [ ] `P2-23` licence review complete with no unresolved item
- [ ] Kinetics term ships as `diagnostic`, never as a Pareto axis (`AMD-002`)
