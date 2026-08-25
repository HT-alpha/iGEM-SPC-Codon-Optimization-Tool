# Assignment Traceability and Acceptance

Machine-readable companion: [`traceability.yaml`](traceability.yaml), checked by `planning/tools/validate_planning_package.py --check trace`.

---

## 1. Deliverables (brief §2) → artifacts

| Requirement | Artifact(s) |
| --- | --- |
| **2.1-a** Repository layout, one module per layer, independently testable and swappable | [`../01-modules/01-repository-layout.md`](../01-modules/01-repository-layout.md), [`../01-modules/modules.yaml`](../01-modules/modules.yaml) |
| **2.1-b** Explicit dependency graph, cycles caught at design time | [`../01-modules/02-module-dependency-graph.md`](../01-modules/02-module-dependency-graph.md) + a DAG assertion in the validator |
| **2.1-c** Language and runtime split, justified | [ADR-0001](../07-decisions/ADR-0001-language-and-runtime-split.md) |
| **2.2-a** Canonical codon-indexed, frame-safe sequence representation | [`sequence.schema.json`](../02-schemas/sequence.schema.json) |
| **2.2-b** Landscape query interface schema | [`landscape_query.schema.json`](../02-schemas/landscape_query.schema.json), [`landscape_response.schema.json`](../02-schemas/landscape_response.schema.json) |
| **2.2-c** Host data, construct context, Pareto front + rationale tags | [`host_profile.schema.json`](../02-schemas/host_profile.schema.json), [`riboseq_track.schema.json`](../02-schemas/riboseq_track.schema.json), [`construct_context.schema.json`](../02-schemas/construct_context.schema.json), [`pareto_front.schema.json`](../02-schemas/pareto_front.schema.json), [`rationale_tags.json`](../02-schemas/rationale_tags.json) |
| **2.2-d** Artifact versioning for reproducibility | [`model_manifest.schema.json`](../02-schemas/model_manifest.schema.json), [`landscape_lockfile.schema.json`](../02-schemas/landscape_lockfile.schema.json), [ADR-0009](../07-decisions/ADR-0009-artifact-versioning-and-resolution.md) |
| **2.3-a…f** Pseudocode for the six named algorithms | [`../03-algorithms/01`](../03-algorithms/01-epistatic-surface-fit.md)…[`06`](../03-algorithms/06-constraint-repair.md) |
| **2.3-g** Complexity analysis and performance budgets | [`../03-algorithms/08-complexity-and-performance-budgets.md`](../03-algorithms/08-complexity-and-performance-budgets.md) |
| **2.4-a** REST/CLI/library surface | [`openapi.yaml`](../04-api/openapi.yaml), [`../04-api/01-rest-api.md`](../04-api/01-rest-api.md), [`../04-api/02-cli-and-library.md`](../04-api/02-cli-and-library.md) |
| **2.4-b** Active-learning ingestion, non-blocking recalibration | [`../04-api/03-active-learning-ingestion.md`](../04-api/03-active-learning-ingestion.md), [`observation.schema.json`](../02-schemas/observation.schema.json) |
| **2.4-c** Query-time model version selection | [`../04-api/04-model-resolution.md`](../04-api/04-model-resolution.md) |
| **2.5-a** Test plan mapping §5 validation units, with thresholds | [`../05-testing/01-test-plan.md`](../05-testing/01-test-plan.md) §2, [`../05-testing/02-acceptance-thresholds.md`](../05-testing/02-acceptance-thresholds.md) |
| **2.5-b** Synthetic-data suite with known ground truth | [`../05-testing/03-synthetic-data-suite.md`](../05-testing/03-synthetic-data-suite.md) |
| **2.6-a** Ticket-level task breakdown per phase | [`../06-build/01`](../06-build/01-phase-1-tickets.md)…[`05`](../06-build/05-phase-5-tickets.md) |
| **2.6-b** Risk register | [`../06-build/06-risk-register.md`](../06-build/06-risk-register.md) |

---

## 2. Open questions (brief §3) → answers

Each has an explicit, justified answer. Where a *parameter* within the answer genuinely requires data that does not yet exist, the answer fixes the **selection procedure, the grid, the default and the fallback**, and names the experiment that fills in the number — which is what `BRIEF-4` asks for and what `BRIEF-5.9` accepts.

| Q | Question | Answer | Deferred parameter, if any |
| --- | --- | --- | --- |
| 3.1 | Epistasis window and order | [ADR-0003](../07-decisions/ADR-0003-epistasis-window-and-order.md): host-dependent window selected by held-out likelihood on a fixed grid with a one-standard-error rule and a data-volume tiering table; pairwise order with triplets behind a hard per-host gate; data-poor hosts shrink to a pooled cross-host prior and lose the dwell head rather than fitting it badly. Default `W = 10`, anchored on the ribosome footprint. | The fitted `W*` per host and whether the triplet gate ever passes → `EXP-08`; whether pooling helps → `EXP-01` |
| 3.2 | Oracle vs. differentiable surrogate | [ADR-0002](../07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md): black-box oracle is the contract; exact deltas (L1) — not gradients — are what make local search fast, at ~640× cheaper per move at 2 kb and ~4 800× at 15 kb. Gradients stay optional and unused. | none |
| 3.3 | Resource-ledger granularity | [ADR-0004](../07-decisions/ADR-0004-resource-ledger-granularity.md): per-isoacceptor pool with wobble shares, plus aggregate ribosome occupancy, at a single steady state. Amino-acid and energy pools are *not* tracked because they are invariant under synonymous redesign. Capacity by three tiers down to a genome-only fallback with explicitly different output semantics. | Whether per-timepoint resolution changes designs → `EXP-03`; whether tier 3 predicts burden → `EXP-04` |
| 3.4 | Uncertainty representation and propagation | [ADR-0005](../07-decisions/ADR-0005-uncertainty-and-pareto-ranking.md): per-term base estimators (analytic / deep ensemble / bootstrap / MC dropout) wrapped in split conformal, Mondrian for per-position. Propagated by probabilistic dominance with a modelled shared/idiosyncratic error split, a reported robustness score, and uncertainty-aware front reduction. | Empirical coverage → `TST-UNC-002`/`GATE-UNC-001` |
| 3.5 | Active-learning data isolation | [ADR-0006](../07-decisions/ADR-0006-active-learning-data-isolation.md): tenant data never enters a shared weight set; adapters over a public base; four enforcement layers ending in a promotion-time manifest assertion; opt-in per-observation community pool with review and revocation semantics. | none |
| 3.6 | Seed import from existing tools | [ADR-0007](../07-decisions/ADR-0007-external-tool-seed-import.md): yes, in Phase 1, because the benchmark harness needs the same interface. Strict acceptance with a positional protein diff on mismatch; 20 % population cap; provenance surfaced in the front. | none |
| 3.7 | Compute budget and scaling | [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md): single workstation (8–16 cores, 32 GB, optional GPU); a three-tier degradation ladder in which `fast` drops L0 terms entirely, restricts structure to the initiation window and skips joint refinement; auto-downgrade with a warning rather than failure. | Per-tier `E_max` → `EXP-11` |

---

## 3. Format constraints (brief §4)

| Constraint | How it is met |
| --- | --- |
| **4.a** Schemas in a concrete serializable form, reviewable line by line | 15 JSON Schema (draft 2020-12) files plus a closed tag vocabulary, each with a worked example, all machine-validated; plus typed interface stubs. See [`../02-schemas/README.md`](../02-schemas/README.md). |
| **4.b** Every departure flagged with a one-line rationale | 28 `DEP-*` markers, flagged in place and collected in [`04-departure-index.md`](04-departure-index.md). |
| **4.c** No production implementation code | The package contains schemas, pseudocode listings, interface stubs and one validator script. The validator is planning-package tooling — it checks this document's own consistency and contains no CodonForge domain logic. Pseudocode is written in a deliberately non-executable style. |
| **4.d** Name the experiment where data does not yet exist | 11 entries in [`../08-experiments/experiment-register.md`](../08-experiments/experiment-register.md), each with a **pre-registered decision rule** and a fallback. |

---

## 4. Acceptance criteria (brief §5)

| Criterion | Status | Evidence |
| --- | --- | --- |
| **5.8** Every module has an interface, an I/O schema and ≥1 specified test | Met | `modules.yaml` declares a `public_contract` and ≥1 test id per module; the validator asserts both the file exists and the test id is defined in the test plan. |
| **5.9** Every open question has an explicit answer or a named resolving experiment | Met | §2 above. All seven have an answer; five additionally name an experiment for a parameter inside the answer, each gated to run before the phase that needs it. |
| **5.10** Phase 1 startable immediately with no further design discussion | Met | [`../06-build/01-phase-1-tickets.md`](../06-build/01-phase-1-tickets.md): 35 tickets, each with module, deliverable, acceptance tests, dependencies and size. Every algorithm Phase 1 needs has pseudocode; every schema it needs is written and validated. |
| **5.11** The risk register names the technical bet and its falsifier for each of the four unaddressed capabilities | Met | [`../06-build/06-risk-register.md`](../06-build/06-risk-register.md) §1: `RISK-EPI-001`, `RISK-BURD-001`, `RISK-SYS-001`, `RISK-KIN-001`, each with bet, falsifier, mitigation and fallback. |

---

## 5. What this package deliberately does not contain

- **Backend implementation code.** Out of scope by `BRIEF-4.c`.
- **The architecture proposal itself.** The companion document is the fixed target, not part of this deliverable. Where it was silent, this package decided and flagged the decision; where it was unworkable, [`03-architecture-amendments.md`](03-architecture-amendments.md) says so and proposes a specific change.
- **Fitted models or real datasets.** Dataset *manifests* and their acceptance thresholds are specified; the data itself is Phase 2–3 work.
- **Calendar estimates.** Ticket sizes are scope classes (`S`/`M`/`L`/`XL`), not schedules.
