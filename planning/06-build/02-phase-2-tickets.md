# Phase 2 — Structure + Burden Pareto Objectives

**What Phase 2 delivers:** the multi-objective machinery, the two biophysical objectives that do not need learned models, the burden model, the decision-support layer, and the benchmark harness. At the end of Phase 2 CodonForge is a credible multi-objective codon optimizer that can be compared head-to-head against LinearDesign and CDSfold — which is what makes every later claim measurable.

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P2-01` | Indicator math | `cf-metrics` | Fast non-dominated sort, crowding distance, Das–Dennis (single and two-layer), WFG hypervolume, IGD+, ε-indicator | `TST-METR-001`, `TST-METR-002`, `TST-METR-003` | Phase 1 | M |
| `P2-02` | NSGA-II loop | `cf-search` | The `M ≤ 3` path, with constraint-domination and the `UNKNOWN` rule | `TST-SRCH-007` | `P2-01` | M |
| `P2-03` | NSGA-III loop | `cf-search` | Normalization (ideal point, extreme points, intercepts, degeneracy guards), reference-direction association, niching selection; automatic `M ≥ 4` switch | `TST-METR-005`, `TST-METR-006`, `TST-SRCH-018` | `P2-02` | L |
| `P2-04` | Termination and progress | `cf-search` | Hypervolume-stagnation termination; `hypervolume_history`; per-tier progress reporting; **tier 3 wired into the job's execution path** | `TST-SRCH-019`, `TST-CLI-008`, `contract_*` against `job_status.schema.json` | `P2-03` | S |
| `P2-05` | Batch-yield protocol | `cf-landscape-api` | The resumable `step`/`absorb`/`snapshot`/`resume` state machine and the `score_batch` seam, exercised end-to-end by a **stub** L0 term. Arrow IPC and the `cf-worker` process are deferred (`DEP-LAND-004`) | `TST-LM-001` applied to the stub; one boundary crossing per generation asserted | `P2-03` | L |

### `P2-05` — synchronous by specification, and what that settles

§6's pseudocode is explicit that `results <- WORKER.SCORE_BATCH(chosen)` blocks: the batching amortizes *model* cost over many sequences rather than overlapping work with the search. So the protocol needs no worker pool, which is the capability `P1-29` declined and `CF-CLI-003` still reports absent — the audit question about an undeclared concurrency dependency has a clean answer.

**"One crossing per generation" is a property of the type, not a discipline.** `step` refuses while a batch is outstanding and `advance` refuses too, so a second crossing in one generation is a `CF-LAND-006` rather than something a caller has to remember not to do.

**Truncation defers rather than discards.** §6's component (d): an individual dropped for budget would be `UNKNOWN` forever with nothing in the output saying it had ever been a candidate. Deferred items go to the head of the next generation's queue, and an item deferred and re-proposed is not scored twice.

**A protocol snapshotted mid-batch resumes still waiting**, which is the state a checkpoint is most likely to be taken in — the search is idle precisely while the scorer runs. A resume that forgot the outstanding batch would step again and issue a second crossing for the same generation.
| `P2-05a` | Tier-3 batch consumption | `cf-search` | **Done.** Joins `P2-06`'s ranker, `P2-05`'s protocol and the generation loop, which had each been tested against themselves and never together. Adds the `L0Scorer` seam, `PyL0Scorer` for a real PyO3 crossing, and the protocol's state to the tier-3 checkpoint so `CF-SRCH-016` survives a resume. `DEP-SRCH-012`: the L0 axis is appended by tier 3 rather than declared by the evaluator, and the crossing is the **PyO3** boundary `ADR-0001` actually caps — `DEP-LAND-004`'s out-of-process worker stays deferred | `TST-SRCH-023`, `TST-SRCH-024` | `P2-05`, `P2-06` | L |
| `P2-06` | Surrogate rationing | `cf-search` | Ridge surrogate on k-mer + L1-term features; the four-way `chosen` set including the mandatory random 10 %; periodic refit | `TST-SRCH-005` | `P2-05` | M |
| `P2-07` | Front certification | `cf-search` | Exact re-evaluation of the reduced front before emission; warning when certification changes the front | `TST-SRCH-005` | `P2-06` | S |
| `P2-08` | Structure kernel | `cf-biophys-structure` | **Re-scoped by `AMD-013`.** Bind the windowed-accessibility term to a user-installed **ViennaRNA**, rather than implementing folding from scratch; record engine and parameter-set identity in the lockfile. A from-scratch kernel is no longer required, since the licence that ruled ViennaRNA out as a default did not in fact do so for a freely-distributed release | `TST-STRU-004`, `TST-STRU-005` | Phase 1, `P2-07a`, `P2-09` | L |
| `P2-07a` | Structure parameter provenance | `cf-biophys-structure` | `structure_parameters` schema; provenance and citation fields; digest recorded in the landscape lockfile; refusal with `CF-STRU-001` when absent | `TST-STRU-005` | Phase 1 | M |
| `P2-09` | Structure engine seam | `cf-biophys-structure` | `StructureEngine` trait; runtime detection of user-installed ViennaRNA / LinearFold; **no vendoring**; engine recorded in the lockfile. **Now precedes `P2-08`** per `AMD-013`: if the default engine is a binding rather than an implementation, the seam is what `P2-08` binds through, so building the kernel first would mean building it against nothing | `TST-STRU-001`, `TST-STRU-003` | `P2-07a` | M |
| `P2-10` | Incremental structure cache | `codonforge.predictors` | Window-keyed cache; bounded recomputation per 3-nt edit; exactness against full recomputation. **Lives with the engine (Python), not in the Rust crate** — `AMD-013` moved the engine out of process, and a cache of engine results cannot sit on the other side of that boundary | `TST-STRU-002`, `TST-STRU-006` | `P2-08` | M |
| `P2-11` | Structure terms | `cf-biophys-structure` | **Capability corrected by the P2-10 audit (`DEP-STRU-003`).** `structure.local_accessibility` was declared **L1**, which assumed the in-process Rust kernel `AMD-013` removed; an out-of-process Python engine cannot serve an exact per-move delta inside the crossing budget, so it is **L0**. Promoting it from `diagnostic` to a Pareto axis needs `P2-05a`. `structure.mfe_global` stays L0, `exact` tier only. Also: `local_structure` constraint support. **Did not activate `P1-26` (`DEP-SRCH-013`):** the term is windowed (`R = 14`) and out of core, so it never enters the tier-2 objective whose partition the overlapping-block fallback would repair. | `contract_*`; the term is an axis only once tier 3 can batch | `P2-10`, `P2-05a` | M |
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
| `P2-23` | Licence and patent review | — | Written review of every dependency's licence and of the LinearDesign patent position per `AMD-013`; confirmation that nothing is vendored; `NOTICE` file carrying ViennaRNA's condition-2 credit (authors + Institute for Theoretical Chemistry, University of Vienna) and LinearFold's required paper citation. **Must be signed off by a faculty advisor, PI, or institutional tech-transfer/legal contact** — `AMD-013`'s patent reasoning is an engineer's reading of claim text, not a legal opinion, and grant status could not be verified | Signed-off review document in `docs/`, `NOTICE` present | `P2-09` | S |

### `P2-04` — what wiring tier 3 into the job exposed

**The pipeline is now monotone, and `CF-ORCH-007`'s warning correctly stops firing.** `P1-30` found the objective-free mutation loop moving *away* from tier 1's optimum, with the front reporting the loss. Tier 3 is elitist — NSGA-II keeps the best front — so the regression it was compensating for is gone. The archive guarantee is unchanged and still tested directly against constructed regressions, because a future stage may well regress; what changed is that the pipeline no longer produces one.

**`serde_json` needs its `float_roundtrip` feature, and every checkpoint depended on it.** Without it the parser takes a fast path that can land one ulp from the value written, so `0.40476190476190477` came back as `0.4047619047619048`. Every checkpoint here carries floats — objective vectors, hypervolume history, violations — so a resumed run diverged from an uninterrupted one by a rounding error. Found by a tier-3 resume test, not by a control, and it is exactly what `GATE-REPRO-001` forbids.

**The generation cap has to be checked before a generation runs, not after.** Checking afterwards let a resumed segment run one past it: an uninterrupted request for nine against a cap of eight stopped at eight, while the same run split 8 + 1 ran a ninth. Both were internally consistent and disagreed with each other.

**The tier-3 checkpoint carries the tabu ledger.** The fourth thing found missing from a checkpoint in this repository, and the first found by a test rather than by a negative control.

**A single-objective fixture cannot exercise tier 3 at all.** Tier 1's DP is exact for an additive objective, so the population starts at the optimum, converges to one sequence, and hypervolume is constant — every assertion about the population loop passes by doing nothing. The orchestrator and CLI fixtures are two-objective for this reason.

**`current_front_summary.best_per_objective` is a per-axis extremum**, filled from tier 3's own report rather than recomputed, so the number a UI plots is the one termination was decided on. It is deliberately not a third notion of "best" beside `front()` and `best()`: those answer "which designs", this answers "how far did each axis get".

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
- [ ] **Release gate live**: `tools/check_release_clearance.py` runs as the `legal-clearance` job in `.github/workflows/wheels.yml`, blocking every tagged release until `docs/legal-clearance.yaml` records sign-off from a faculty advisor, PI, or institutional tech-transfer/legal contact. Expected to *fail* until someone signs; that is the gate working, not a defect to route around. `AMD-013`.
- [ ] `NOTICE` credits ViennaRNA (authors and the Institute for Theoretical Chemistry, University of Vienna) and LinearFold (its required paper citation), quoting both licences rather than paraphrasing them
- [ ] Kinetics term ships as `diagnostic`, never as a Pareto axis (`AMD-002`)
