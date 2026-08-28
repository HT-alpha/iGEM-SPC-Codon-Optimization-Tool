# Test Plan

Satisfies: `BRIEF-2.5-a` (mapping each validation unit in the architecture's §5 to specific automated tests or benchmark scripts) and `BRIEF-5.8` (every module has at least one specified test).
Thresholds: [`02-acceptance-thresholds.md`](02-acceptance-thresholds.md). Synthetic fixtures: [`03-synthetic-data-suite.md`](03-synthetic-data-suite.md). Benchmarks: [`04-benchmark-harness.md`](04-benchmark-harness.md).

---

## 1. Tiers and where they run

| Tier | Naming | Budget | Runs |
| --- | --- | --- | --- |
| **unit** | `unit_*` | < 100 ms each, no I/O | pre-commit, every PR |
| **property** | `property_*` | < 30 s per property | every PR |
| **contract** | `contract_*` | < 1 s each | every PR — validates a module's public contract against its JSON Schema using golden examples |
| **integration** | `it_*` | < 5 min total | every PR |
| **benchmark** | `bench_*` | minutes to hours | nightly; weekly for the full sweep |
| **performance regression** | `perf_*` | ~30 min | nightly |

`contract_*` is the tier that makes "independently swappable" real: a module can be replaced if and only if its replacement passes the same contract tests against the same golden documents.

---

## 2. Architecture §5 validation units → tests

This is the mapping `BRIEF-2.5-a` asks for. Each of the four named validation units becomes a concrete, automated, thresholded check.

| Architecture §5 validation unit | Test / benchmark | Gate |
| --- | --- | --- |
| Held-out ribosome-profiling prediction | `TST-EPI-003` (unit-scale) + `bench_epistasis_heldout` (full) | `GATE-EPI-001`, `GATE-EPI-002` |
| MFE/CAI trade-off benchmarking against LinearDesign / CDSfold | `bench_mfe_cai_frontier` | `GATE-BEN-001`, `GATE-BEN-002` |
| Burden-model back-testing against a Codon-Harmony-Index-style dataset | `bench_burden_backtest` | `GATE-BURD-001` |
| Pareto-front hypervolume comparisons | `bench_hypervolume_vs_baselines` | `GATE-BEN-003` |

Two supporting units the architecture implies but does not name, added here because without them the four above are not interpretable:

| Added unit | Test | Gate | Why it is needed |
| --- | --- | --- | --- |
| Uncertainty calibration | `TST-UNC-002` | `GATE-UNC-001` | Every hypervolume and dominance claim in the plan depends on ranked candidates, and [ADR-0005](../07-decisions/ADR-0005-uncertainty-and-pareto-ranking.md) makes ranking depend on calibrated intervals. |
| End-to-end performance | `TST-PERF-001` | `GATE-PERF-001` | `BRIEF-2.3-g` sets performance budgets; without a regression suite they decay silently. |

---

## 3. Test catalogue

### `cf-types`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-CORE-001` | property | Sequence invariants I1–I5: length/codon-count agreement, translation matches the declared protein, protein invariance under every exposed operation, locked spans in range, digest stability. |
| `TST-CORE-002` | unit | Codon-index bijection (`16·b₀+4·b₁+b₂`) round-trips for all 64 codons, and the 27 defined NCBI translation tables match a vendored reference (ids 7, 8, 17-20 are retired and must be rejected, not defaulted). |
| `TST-CORE-003` | unit | Canonical JSON (RFC 8785) serialization is stable across key insertion order, float formatting and platform. A digest that changes with map iteration order silently breaks every cache and every lockfile. |

### `cf-metrics`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-METR-001` | property | Non-dominated sorting and hypervolume agree with brute force on 10⁴ random small fronts (`M ≤ 5`, `N ≤ 30`). |
| `TST-METR-002` | unit | Das–Dennis directions: count equals `C(M+p−1, p)`, all sum to 1, all unique; the two-layer construction covers the interior. |
| `TST-METR-003` | property | Hypervolume is monotone under front improvement and invariant to candidate ordering; IGD+ is zero iff the front equals the reference. |

### `cf-landscape-api`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-LAPI-001` | contract | Every golden query/response document round-trips Rust → JSON → Python → JSON byte-identically. Catches the type-duplication drift that [ADR-0001](../07-decisions/ADR-0001-language-and-runtime-split.md) accepts as a cost. |
| `TST-LAPI-002` | contract | A predictor's declared `capability` matches its actual method surface: an L1 claim without a working `delta_table` fails. |
| `TST-LAPI-003` | property | The three contract invariants composition depends on. (a) *Ordering*: `objective_vector.order` equals the caller's requested axis order and is independent of registry insertion order, so the response never leaves the order to be inferred from `terms` key iteration. (b) *Seed determinism*: two evaluations of the same query with the same `seed` agree bit-for-bit, which is what makes `TST-LAND-001`'s "a cache hit returns exactly what a cold evaluation returns" a meaningful claim rather than a vacuous one. (c) *Term independence*: evaluating a set of terms in any permutation yields identical per-term values, so composition may not be order-sensitive. Added during the `P1-20` sequencing audit — `TST-LAPI-001`/`002` check the shape of the contract, and none of these three is a shape property. |

### `cf-constraints`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-CONS-001` | property | Incremental detectors agree exactly with a naive `O(n)` rescan over 10⁵ random edits, for every constraint class. |
| `TST-CONS-002` | unit + property | Repair minimality on planted violations with a known unique fix; and, over 1 000 seeds on a violation with several minimal fixes, more than one distinct fix is returned (the diversity property). |
| `TST-CONS-003` | property | Repair always terminates and never exceeds its declared work bound, including on adversarial overlapping and mutually-conflicting violations. |
| `TST-CONS-004` | property | Repair never changes the translated protein and never touches a locked codon (10⁶ cases). |
| `TST-CONS-005` | integration | Anti-thrashing: a scripted mutation/repair cycle designed to thrash is broken within `tabu_tenure` generations. **Scheduled with `P1-18b`/`P1-23`** — it needs a mutation operator and a generation clock, which `cf-constraints` does not have. The ledger's own behaviour (record, tick, leak) is unit-tested with `P1-18a`. |
| `TST-CONS-006` | unit | Each of the five feasibility checks fires on a purpose-built infeasible input and stays silent on a feasible one. One test per check (`locked_span`, `counting_bound`, `per_position`, `automaton_reachability`, `gc_attainability`), plus one asserting a satisfiable set produces no finding — so a check that fired on everything would fail. |
| `TST-CONS-007` | integration | **Legacy differential (`P1-27`).** The checker agrees with `codon_6.0.py` on feasibility for 10⁴ synonymous variants, with the comparison configured to the legacy semantics. The oracle is the prototype itself, loaded and called — not a transcription, which could share a bug with the checker and agree for the wrong reason. Requires both outcomes to be well represented, so a checker that called everything feasible would fail. |
| `TST-CONS-008` | unit | A blocking violation confined to immutable sequence is reported but marked non-actionable, and does not make a candidate infeasible. See `docs/legacy-parity.md` §3. |
| `TST-CONS-010` | unit | A violation confined to a locked span is reported but not actionable, does not make a candidate infeasible, and *is* reported by feasibility check 1. `Violation::actionable` was computed against the whole coding region, so a forbidden motif spelled entirely by locked codons satisfied it — and check 1, which exists to report that case once at submission, uses the same predicate and stayed silent too. The case fell through both layers and out into tier 1, whose beam emptied and reported `CF-SRCH-006`: a claim about the design space rather than about the input. Also asserts that a construct with no locks is unaffected, so the `P1-27a` legacy differential is untouched. Added by `P1-25`'s audit. |
| `TST-CONS-009` | property | `check_region` equals `check` restricted to the same span, over randomized constructs and every codon-aligned region. Exactness is the method's whole value — a windowed check that is merely usually right yields deltas that are usually right, and a search steered by those converges confidently on the wrong answer. Includes the case of a homopolymer run wider than the fixed radius, and asserts that a k-mer repeat rule reports *no* bounded neighbourhood rather than a large one. Fixture density is asserted, so the test cannot pass by comparing empty lists. Added by `P1-21`, which needs it for `constraint.violation_count`'s exact delta. |

### `cf-modelreg`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-MREG-001` | unit | The resolution ladder resolves in documented order and records `resolved_by` correctly at every rung. |
| `TST-MREG-002` | unit | Replaying a lockfile with a missing digest fails with `CF-RESOLVE-006` and never substitutes a newer artifact. |
| `TST-MREG-003` | unit | Fallback admissibility is enforced; admissible fallbacks set `fallback: true`, inflate `sigma`, and set `calibrated: false`. |
| `TST-MREG-004` | unit | Mixing `contract_major` values across terms in one job is refused. |
| `TST-MREG-005` | unit | A `yanked` artifact never resolves through a channel and always resolves by digest. |

### `cf-testkit`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-SYN-001` | unit | The toy design space (8 codons, ≤3 synonyms, 6 561 sequences) enumerates exhaustively and its exact Pareto front is computable, giving every search test a ground truth. |
| `TST-SYN-002` | unit | The synthetic Potts generator produces sequences whose empirical conditionals match the generating parameters. |
| `TST-SYN-003` | unit | The two-cistron contention instance has an analytically known joint optimum that is strictly better than the best independent-optimization result. |
| `TST-SYN-004` | unit | Planted-violation generators produce sequences with exactly the intended violations and a known minimal repair distance. |

### `cf-ingest`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-INGE-001` | property | FASTA/GenBank → `CodonSequence` → FASTA/GenBank round-trips; frame is preserved; RNA input is transcoded; ambiguity codes and lower case are handled or rejected explicitly. |
| `TST-INGE-002` | unit | Host loaders populate `data_completeness` correctly, and Ribo-seq QC gates flag low periodicity, low depth and replicate discordance. |
| `TST-INGE-003` | unit | Seed importers accept valid output from each supported tool and reject a protein-mismatched seed **with a correct positional diff**. |

### `cf-epistasis`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-EPI-001` | integration | Parameter recovery from a known synthetic Potts surface: Pearson `r ≥ 0.95` on `J` in zero-sum gauge. Catches gauge, indexing and gradient bugs with no biology involved. |
| `TST-EPI-002` | property | **Delta correctness**: `delta == score(after) − score(before)` within `1e-9` over 10⁴ random substitutions. The highest-value test in the package — a wrong delta corrupts every tier-2 result while leaving tier-3 results looking plausible. |
| `TST-EPI-003` | benchmark | Held-out dwell prediction beats the codon-usage baseline (`GATE-EPI-001`). |

### `cf-burden`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-BURD-001` | property | Ledger conservation: total amino-acid demand is bit-identical before and after 10⁵ random synonymous substitutions. |
| `TST-BURD-002` | property | Incremental and full ledger agree to `1e-12` over 10⁵ edits. |
| `TST-BURD-003` | integration | Frozen-availability delta error stays below 1 % of the term's reference range. |
| `TST-BURD-004` | property | `g` is monotone in each `u_α` and bounded in `[0, g_max]` over 10⁶ random ledgers. |
| `TST-BURD-005` | unit | Tier-3 semantics: unit is `dimensionless`, `calibrated` is `false`, and `low_confidence_host` is attached. Keeps the honesty requirement of [ADR-0004](../07-decisions/ADR-0004-resource-ledger-granularity.md) from eroding. |

### `cf-biophys-structure`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-STRU-001` | integration | Engine parity: on a reference sequence set, the native kernel's accessibility correlates with an installed ViennaRNA reference above a fixed threshold. Skipped, not failed, when ViennaRNA is absent. |
| `TST-STRU-002` | property | Incremental accessibility equals full recomputation to `1e-9` over 10⁴ edits. |
| `TST-STRU-003` | unit | With no optional engine installed, the native engine is selected, the lockfile records it, and no import of an absent package is attempted. |

### `cf-biophys-kinetics`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-KIN-001` | unit | Desired pause windows land at `d+25 … d+45` and never overlap a boundary. |
| `TST-KIN-002` | integration | Reconciliation with the epistasis term: a pause inside a designated window does not increase `dwell_deviation`; one 100 codons away does. |

### `cf-biophys-lm`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-LM-001` | property | `score_batch` equals a loop over `score` to within float tolerance, for every batch size including 1 and `E_max`. |
| `TST-LM-002` | unit | Given a fixed seed, MC-dropout uncertainty is reproducible across runs and across batch composition. |

### `cf-landscape`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-LAND-001` | property | Cache correctness: a cache hit returns exactly what a cold evaluation returns; the cache key includes the lockfile digest so cross-version contamination is impossible. |
| `TST-LAND-002` | unit | `status: partial` propagates as `UNKNOWN` and is never imputed; the deferred term appears in `diagnostics.terms_deferred`. |
| `TST-LAND-003` | unit | `cargo tree -p cf-landscape` contains no predictor crate — `cf-epistasis`, `cf-biophys-structure`, `cf-biophys-kinetics`, `cf-biophys-lm`, `cf-burden` — named rather than left implicit. The mechanical enforcement of the swappable-layer property. `cf-constraints` is present and permitted: it is a first-party core layer, not a swappable predictor, and `constraint.violation_count` is a built-in term. |
| `TST-LAND-004` | unit | Mixed-capability composition degrades per term, not per query. In a registry where some terms are L1 and at least one is L0-only, a delta query returns exact deltas for every L1 term and falls back to full rescoring for the L0 term *alone*; the number of full rescores equals the number of L0-only terms requested. Guards the no-gradients binding decision from the other side: the cheap path must not be abandoned wholesale because one term cannot take it. Added during the `P1-22` sequencing audit — `TST-LAND-001`/`002`/`003` all pass under a composition that degrades globally. |
| `TST-LAND-005` | unit | A term whose role is `diagnostic` is evaluated and reported but never enters `objective_vector`, and toggling it cannot change a ranking. Enforces `AMD-002` structurally for folding kinetics ahead of `EXP-06`. |
| `TST-LAND-006` | unit | The delta path costs what it claims. An L1 provider's `score` is never called during a sweep of delta queries against a pinned base, and a delta table pins once rather than once per alternative. Asserted on **call counts**, not on values: the defect `DEP-LAND-003` fixed returned correct numbers, reported `full_rescores: 0`, and differed from a full rescore only in cost — so a value assertion cannot see it. Added during the `P1-23` audit. |

### `cf-search`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-SRCH-001` | property | Tier-2 monotonicity: the scalarized surrogate never decreases across a half-sweep, over 500 randomized runs. |
| `TST-SRCH-002` | integration | On the toy exhaustive space, tier 1 + tier 2 find the exact decomposable optimum in ≥99 % of seeds. |
| `TST-SRCH-003` | property | Seam reconciliation never worsens the incumbent (10⁴ configurations). |
| `TST-SRCH-004` | integration | Bit-identical output across 1, 4 and 16 threads. |
| `TST-SRCH-005` | integration | Every emitted candidate carries an exact, non-surrogate evaluation. |
| `TST-SRCH-006` | property | Protein and lock invariance after 10⁶ randomized operator applications, checked after *every* application rather than at the end — a violation introduced and then reverted is invisible to a final check. Run against a genotype with a locked span, since an operator that ignored locks would still preserve the protein and would pass a protein-only assertion. Crossover is covered separately, because it is the operator most likely to write a codon directly. Both halves assert that a useful fraction of applications actually changed something, so neither can pass by doing nothing. |
| `TST-SRCH-012` | unit | Tier 2's red-black independence, asserted structurally rather than through answers. Two observables: the landscape is queried once per scored term per **half-sweep** rather than once per block (six blocks in the fixture, so a per-block implementation differs sixfold), and reversing the order in which same-colour blocks are visited leaves the incumbent hash and the whole surrogate trace unchanged. The second is the decisive one — if any same-colour block could observe another, order would matter — and it is only assertable because each block's RNG stream is keyed on its index rather than on iteration order. A coupled implementation returns a good answer and is simply not the coordinate-ascent step it reports being. Added by `P1-25`. |
| `TST-SRCH-011` | unit | Tier 1's cost, not its answer. The DP takes **one delta table per scored term against the base** and reads increments out of it; asserted on landscape call counts, because a DP that queried per `(state, synonym)` would return the same seeds while issuing `O(L·b·S)` queries — tens of thousands on the test fixture, roughly 120 000 on a 666-codon protein at `b = 64` — and no assertion on values can tell the two apart. Also pins the tier-planning refusal: a term whose provider is L0 is rejected before the tier runs rather than discovered mid-run. Added by `P1-24`, in the shape `TST-LAND-006` established. |
| `TST-SRCH-010` | property | An operator move *is* a landscape delta query: `Move` carries `cf_landscape_api::Edit` values, so `to_delta_target` hands them over rather than translating them, and every proposed move applies cleanly to the base it was proposed against. Guards the interop the exact-delta decision depends on — a genotype recording only its final codons would force an `O(L)` diff to recover the edit list, putting the query back on the expensive path by construction. Asserts that single-substitution moves, the shape tier-2 accept/reject uses, are well represented. Added by `P1-23`. |
| `TST-SRCH-007` | unit | Constraint-domination ordering matches Deb's definition on an enumerated truth table, including `UNKNOWN` cases. |
| `TST-SRCH-008` | integration | Joint-search oscillation is detected and the fully-joint fallback produces a result at least as good as the decomposed path. |
| `TST-SRCH-009` | integration | On the synthetic contention instance, joint search beats independent per-cistron optimization by the analytically known margin. The test that justifies the systems-level layer's existence. |

### `cf-decision`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-DEC-001` | property | Probabilistic dominance agrees with exhaustive enumeration on small discrete predictive distributions, including the correlated case. |
| `TST-DEC-002` | unit | Every rationale tag fires on a fixture that meets its trigger and stays silent on one that does not; every emitted tag is in the closed vocabulary with the right scope. |
| `TST-DEC-003` | property | Front reduction always retains the per-objective extreme points, so a reduced front never hides the best value on any axis. |

### `cf-orchestrator`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-ORCH-001` | integration | The lockfile is sealed before any search work begins, and no query in the job resolves outside it. |
| `TST-ORCH-002` | integration | Checkpoint/resume: a job killed mid-run resumes to a bit-identical result. Compared against an *uninterrupted run of the same seed plan*, state for state — incumbent, tabu ledger and per-generation history — not against itself. A resume that restarted its RNG stream counter at zero would be perfectly reproducible against itself and quietly different from an uninterrupted run, and a final-state comparison would also miss a run that diverged and reconverged. |
| `TST-ORCH-004` | unit | The segment runner: every split of the same generation total produces an identical checkpoint - incumbent, tabu ledger and per-generation history. Seven splits, not one, because a run that resumed correctly at one boundary and not another would pass a single-split test. |
| `TST-ORCH-005` | unit | The tabu ledger survives the checkpoint, and the fixture asserts it was live at the boundary. An empty ledger makes "resume carries it" and "resume drops it" the same test. |
| `TST-ORCH-006` | unit | The ledger's checkpoint snapshot is ordered. `HashMap` iteration order is randomly seeded per process, so an unsorted snapshot surfaces as an *intermittent* resume-comparison failure; asserted directly by inserting the same entries in opposite orders. |
| `TST-ORCH-007` | unit | A malformed checkpoint is refused with `CF-ORCH-001` naming the missing field, rather than partly believed. |
| `TST-ORCH-008` | integration | Every `job_status.warnings` entry carries a `stage` and a `subject`, and the three obligation sources remain separable in one document. |
| `TST-ORCH-009` | integration | A job checkpointed across a model promotion resumes on its sealed digests. The promotion is published *while the job is stopped*, which is the case the non-blocking-recalibration guarantee has to survive and a continuous run never exercises. |
| `TST-ORCH-003` | integration | `codonforge replay --assert-identical` reproduces a front exactly. |

### `cf-active`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-ACT-001` | unit | Each ingestion validation check fires and stays silent appropriately. |
| `TST-ACT-002` | integration | Isolation invariant over a simulated multi-tenant history: no tenant artifact's manifest lists a foreign tenant's dataset. |
| `TST-ACT-003` | integration | **Non-blocking**: a job running while a recalibration promotes a new artifact produces a bit-identical front to the same job with no concurrent promotion. The direct executable test of `BRIEF-2.4-b`. |
| `TST-ACT-004` | integration | An artifact failing any required gate is unreachable through a channel and reachable only by explicit digest. |

### `cf-api` / `cf-cli` / `cf-bench`

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-API-001` | contract | Every response validates against its declared schema, driven from the OpenAPI document. |
| `TST-API-002` | integration | Idempotency: same key + same body replays; same key + different body returns `409`. |
| `TST-API-003` | integration | Tenancy: tenant A cannot read, pin or replay any tenant-B resource. Enumerates every endpoint. |
| `TST-API-004` | unit | Every error path returns `application/problem+json` with a documented `code`. |
| `TST-CLI-001` | integration | REST, CLI and library produce identical lockfiles and fronts for a corpus of requests under a fixed seed. |
| `TST-CLI-002` | unit | Exit codes match the documented table, in particular that infeasible constraints exit `3` and not `1`. |
| `TST-BEN-001` | integration | The benchmark harness reproduces a stored reference result exactly from its lockfile. |
| `TST-BEN-002` | unit | The harness refuses to compare fronts computed at different budget tiers or with different lockfiles. |

### Cross-cutting

| Id | Tier | Assertion |
| --- | --- | --- |
| `TST-UNC-001` | contract | Every term's `Uncertainty` declares `method` and `calibrated`; `method: none` appears only for deterministic terms. |
| `TST-UNC-002` | benchmark | Empirical coverage of the nominal 90 % interval on held-out data lands in `[0.87, 0.93]` (`GATE-UNC-001`). |
| `TST-PERF-001` | perf | End-to-end wall time for the 2 kb and 15 kb reference jobs at all three tiers; fails on a >25 % regression versus the trailing 7-day median. |

---

## 4. What is deliberately not tested by unit tests

Stating this prevents a false sense of coverage:

- **Whether the models are biologically right.** That is what the experiments in [`../08-experiments/experiment-register.md`](../08-experiments/experiment-register.md) are for. A unit test can prove the Potts fitter recovers its own synthetic parameters; only `EXP-01` can say whether the fitted surface means anything.
- **Whether designs actually express better.** `EXP-10`.
- **Absolute predictive accuracy.** Gates are relative, against named baselines — see [`02-acceptance-thresholds.md`](02-acceptance-thresholds.md) §1 for why.
