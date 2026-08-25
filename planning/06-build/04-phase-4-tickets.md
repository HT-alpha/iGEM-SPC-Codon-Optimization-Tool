# Phase 4 — Systems-Level Joint Search

**What Phase 4 delivers:** multi-cistron design against a shared resource ledger — the capability that distinguishes CodonForge from every per-gene codon optimizer — plus the performance work needed to hit the 15 kb budget.

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P4-01` | Multi-cistron end-to-end | `cf-orchestrator` | Construct-level job execution: several editable CDSs, per-cistron locks, construct-scope constraints, per-cistron outputs in one candidate | A 5-cistron construct designs end-to-end with independent per-cistron optimization | Phase 3 | L |
| `P4-02` | Shared ledger | `cf-burden` | Cross-cistron demand accumulation with `f_k` weighting and plasmid copy number; native-proteome baseline demand | `TST-BURD-001` extended to multi-cistron | `P4-01` | M |
| `P4-03` | Availability fixed point | `cf-burden` | Damped iteration to self-consistent availability; convergence tolerance; iteration cap | Fixed point converges on `SYN-03` and on `bm_constructs_20` | `P4-02` | M |
| `P4-04` | Shadow prices | `cf-burden` | Analytic `π_α` for both the burden and elongation channels | Analytic prices match finite differences to `1e-6` | `P4-03` | M |
| `P4-05` | Price-decomposed outer loop | `cf-search` | Alternating optimization with cistrons in parallel against a frozen price vector; convergence test on `‖ΔD‖` | Decomposed path beats independent optimization on `SYN-03` | `P4-04` | L |
| `P4-06` | Oscillation detection and fallback | `cf-search` | Non-decreasing-`‖ΔD‖` detector; damping reduction; fully-joint fallback; `CF-JOINT-001` warning | `TST-SRCH-008` | `P4-05` | M |
| `P4-07` | Joint refinement | `cf-search` | Tier-3 NSGA-III over the concatenated genotype with the ledger evaluated exactly per individual | `TST-SRCH-009` | `P4-05` | L |
| `P4-08` | Structure performance work | `cf-biophys-structure` | Strided interior sampling; per-generation structure rationing under `E_max`; profiling-driven kernel optimization | `GATE-PERF-001` for the 15 kb case | `P4-07` | L |
| `P4-09` | `EXP-11` runner | `cf-bench` | `E_max` sweep against final hypervolume and certification-change rate | `EXP-11` decided; per-tier `E_max` set | `P4-08` | M |
| `P4-10` | `EXP-03` runner | `cf-bench` | Two-phase vs. steady-state ledger comparison on an inducible-system benchmark set | `EXP-03` decided; the two-phase extension is built or dropped | `P4-03` | M |
| `P4-11` | `bench_joint_search` + `bench_ablation` | `cf-bench` | Joint vs. independent on real constructs; the full ablation arm set | `GATE-BEN-003`, `GATE-BEN-004` | `P4-07` | M |
| `P4-12` | Distributed evaluation transport | `cf-orchestrator` | Optional remote worker pool for the per-generation evaluation batch; the batch-yield protocol already serializes it, so this is transport only | A cluster run reproduces a single-node run bit-identically | `P4-08` | M |

`P4-12` is optional for the workstation target of [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md) and is included because the batch-yield protocol makes it cheap. The island-model parallelization of the generation loop is deliberately **not** in Phase 4: it changes the search's convergence characteristics and would need its own validation.

## Phase 4 exit checklist

- [ ] Phase 3 checklist still green
- [ ] `TST-SRCH-008`, `TST-SRCH-009` green
- [ ] `GATE-BEN-003`, `GATE-BEN-004`, `GATE-PERF-001`, `GATE-PERF-002` pass
- [ ] `EXP-03` and `EXP-11` decided
- [ ] The 15 kb reference job meets its budget on the `workstation` profile, or the shortfall is documented with the mitigation actually applied
- [ ] Ablation numbers published, including any layer that contributes little
