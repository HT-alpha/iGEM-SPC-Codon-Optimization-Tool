# CodonForge — Backend Planning Package

This is the deliverable for the assignment in *Detailed Backend Planning Assignment: CodonForge* (August 2026): an implementation-ready engineering plan produced **before any backend code is written**. It converts the CodonForge architecture — the epistasis-aware, systems-level, multi-objective codon optimization framework — into schemas, interfaces, pseudocode, complexity budgets, tests and a ticket-level build sequence that an implementation team could pick up without making further architectural judgment calls.

It contains no backend implementation code, per the brief's phase boundary.

---

## Start here

| If you want to… | Read |
| --- | --- |
| Check the assignment is fully covered | [`00-scope/01-assignment-traceability.md`](00-scope/01-assignment-traceability.md) |
| **See where the architecture was found impractical** | [`00-scope/03-architecture-amendments.md`](00-scope/03-architecture-amendments.md) |
| See the answers to the seven open questions | [`07-decisions/`](07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md), summarized in [`00-scope/01-assignment-traceability.md`](00-scope/01-assignment-traceability.md) §2 |
| Start building | [`06-build/01-phase-1-tickets.md`](06-build/01-phase-1-tickets.md) |
| Review the contracts | [`02-schemas/README.md`](02-schemas/README.md) |
| Understand the algorithms | [`03-algorithms/`](03-algorithms/01-epistatic-surface-fit.md) |
| Know what could go wrong | [`06-build/06-risk-register.md`](06-build/06-risk-register.md) |
| Look up a symbol or an id prefix | [`00-scope/02-glossary-and-notation.md`](00-scope/02-glossary-and-notation.md) |

---

## Contents

```
planning/
├── 00-scope/
│   ├── 01-assignment-traceability.md    brief -> artifact map; acceptance status
│   ├── 02-glossary-and-notation.md      terms, symbols, id prefixes
│   ├── 03-architecture-amendments.md    12 places the architecture needs changing
│   ├── 04-departure-index.md            28 local departures, flagged in place
│   └── traceability.yaml                machine-readable, validated
├── 01-modules/
│   ├── 01-repository-layout.md          directory tree, layer -> module map
│   ├── 02-module-dependency-graph.md    the graph, and the 4 cycles it breaks
│   └── modules.yaml                     normative; asserted acyclic
├── 02-schemas/                          15 JSON Schemas + examples + typed stubs
├── 03-algorithms/                       pseudocode + complexity for (a)-(f) + kinetics + budgets
├── 04-api/                              OpenAPI, REST/CLI/library, ingestion, model resolution
├── 05-testing/                          test plan, thresholds, synthetic suite, benchmarks
├── 06-build/                            5 phases of tickets + the risk register
├── 07-decisions/                        ADR-0001..0009
├── 08-experiments/                      11 experiments with pre-registered decision rules
└── tools/validate_planning_package.py   consistency checker for this package
```

---

## The eight decisions that shape everything else

1. **A Rust core for the per-move hot path, Python for training, orchestration and I/O**, with a hard rule that no call crosses the boundary more than once per generation. The move-count arithmetic forces this. ([ADR-0001](07-decisions/ADR-0001-language-and-runtime-split.md))
2. **Exact deltas, not gradients.** The landscape contract is a black-box oracle; what makes local search fast is an optional *exact single-substitution delta* capability, ~640× cheaper per move than rescoring at 2 kb. ([ADR-0002](07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md))
3. **A synonymous-restricted Potts model.** Normalizing the partition function over synonymous sets only makes the fit concave with a unique optimum, keeps it at ~41 k parameters, gives `O(W·S)` exact deltas, and makes the parameters transfer across proteins. ([`03-algorithms/01`](03-algorithms/01-epistatic-surface-fit.md))
4. **A red-black block partition** with `B ≥ W`, which makes same-colour blocks conditionally independent and eliminates the overlap-reconciliation step entirely. ([`03-algorithms/02`](03-algorithms/02-block-local-coordinate-ascent.md))
5. **Shadow-price decomposition for the joint search.** Cross-cistron coupling runs entirely through a ~40-number vector, so cistrons optimize independently against prices, and a final genuinely-joint refinement drops the mean-field approximation. ([`03-algorithms/04`](03-algorithms/04-systems-level-joint-search.md))
6. **Windowed local accessibility instead of global MFE.** Global folding inside the loop costs ~53× the *entire* job's compute budget by itself; the windowed objective is ~300× cheaper per evaluation, incrementally updatable, and better matched to the mechanistic claims. ([`AMD-001`](00-scope/03-architecture-amendments.md))
7. **The lockfile.** One construct delivers reproducibility, query-time model selection, and non-blocking recalibration — because a job that pins its artifacts at submission simply cannot see a mid-flight model promotion. ([ADR-0009](07-decisions/ADR-0009-artifact-versioning-and-resolution.md))
8. **Uncertainty is a ranking input, not an annotation.** Probabilistic dominance over a joint predictive distribution with a modelled shared/idiosyncratic error split, and a reported per-candidate robustness score. ([ADR-0005](07-decisions/ADR-0005-uncertainty-and-pareto-ranking.md))

---

## Validating this package

```bash
pip install jsonschema pyyaml
python3 planning/tools/validate_planning_package.py -v
```

The checker verifies that the schemas are valid and their examples conform, that the module dependency graph is acyclic with no rank inversion, that every cross-referenced identifier resolves, that every Markdown link resolves, that the OpenAPI document parses and its external `$ref`s exist, that every requirement in the brief maps to an artifact, and that every rationale tag used is in the closed vocabulary.

It is a review aid, not backend code — it contains no CodonForge domain logic beyond the genetic code needed to check the worked examples translate correctly.

---

## Relationship to this repository

The repository currently contains a working prototype (`codon_6.0.py`) of the hard-constraint layer: a stack-based backtracking solver with auto-detected forbidden motifs and overlapping k-mer repeat limits. The plan treats it as a real asset rather than something to discard:

- Its **overlapping** k-mer counting semantics are preserved explicitly in `constraint_set.schema.json` (`count_overlapping`), because non-overlapping counting undercounts exactly the tandem repeats a synthesis vendor rejects.
- Ticket `P1-27` ports it as `baseline_backtracking` and keeps it as a **differential-test oracle** for the new constraint layer — an independently-written implementation of the same semantics is a genuine correctness check.
- Its failure mode — exhausting a restart budget and reporting only `Exceeded max restarts` — is what motivates the five-check feasibility analyser (`AMD-010`), which always names a constraint to relax.
