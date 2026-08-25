# ADR-0008 — Target Hardware Profile and the Degradation Ladder

- **Status:** Accepted
- **Resolves:** open question 3.7 (`BRIEF-3.7`)
- **Related:** [`../03-algorithms/08-complexity-and-performance-budgets.md`](../03-algorithms/08-complexity-and-performance-budgets.md), `AMD-001`

## Question as posed

> What is the target hardware profile (single workstation vs. cluster) for a typical job, and where does the design need to make different algorithmic choices (e.g. reduced population size, approximate structure prediction) to stay within that budget?

## Decision

**The target is a single workstation: 8–16 physical cores, 32 GB RAM, optional single consumer GPU. A cluster is supported but never required.** Everything in the plan — the Rust core, windowed folding, surrogate-rationed L0 evaluation, price-decomposed joint search — exists to make that target achievable.

### Why workstation-first

The intended users are academic and iGEM labs and small synthetic-biology groups. A tool that needs a cluster is a tool most of them cannot run, and one that is hard to reproduce results with. The architecture's open-source goal implies a locally-runnable tool; a cluster requirement would quietly contradict it.

Model *training* is a different matter and is explicitly not workstation-scoped: fitting an epistatic surface is a per-host, per-month operation that may use a bigger machine, and its output is a 164 KB artifact anyone can then use.

### Reference profiles

| Profile | Cores | RAM | GPU | Intended use |
| --- | --- | --- | --- | --- |
| `laptop` | 4 | 8 GB | none | `fast` tier, single gene |
| **`workstation` (target)** | 8–16 | 32 GB | optional | `standard` tier, single gene and multi-cistron |
| `server` | 32–64 | 128 GB | 1 | `exact` tier, batch jobs, benchmark suite |
| `cluster` | many nodes | — | — | Model training, benchmark sweeps, `exact` tier at 15 kb. Never required for a user job. |

### Budgets (from [`../03-algorithms/08-complexity-and-performance-budgets.md`](../03-algorithms/08-complexity-and-performance-budgets.md) §4)

| Case | Tier | Target | Estimated |
| --- | --- | --- | --- |
| Single gene ≤2 kb | `fast` | ≤ 30 s | ~8 s |
| Single gene ≤2 kb | `standard` | ≤ 600 s (GPU) / ≤ 900 s (CPU) | ~142 s / ~275 s |
| Multi-cistron ≤15 kb | `fast` | ≤ 120 s | ~50 s |
| Multi-cistron ≤15 kb | `standard` | ≤ 3 600 s on 8 cores | ~885 s |
| Either | `exact` | best effort, no target | 35 min / ~4 h |

`GATE-PERF-001` makes these acceptance criteria for Phase 4, not aspirations, and `TST-PERF-001` fails the nightly build on a >25 % regression.

## The degradation ladder

The tier is the single user-facing knob, and it changes *algorithmic choices*, not just parameters.

| | `fast` | `standard` | `exact` |
| --- | --- | --- | --- |
| **L0 terms** | none — L1 terms only | LM at `E_max = 64`/generation via surrogate rationing | LM on the whole population; global MFE on the reduced front |
| **Structure engine** | initiation window only (~100 nt) | full-CDS windowed accessibility (`w = 80`) | windowed + global MFE re-rank of the front |
| **Population `N`** | 64 | 212 | 320 |
| **Generations `G`** | 80 | 300 | 600 |
| **Tier-2 instances** | 32 | 200 | 320 |
| **Beam width `b`** | 16 | 64 | 128 |
| **Uncertainty** | point estimates, `calibrated: false` | 5-member ensemble + conformal | 5-member ensemble + Mondrian conformal + 2 000-sample probabilistic dominance |
| **Epistasis window** | host default, capped at 5 | host-fitted | host-fitted |
| **Joint search** | price decomposition only, no joint refinement | decomposition + joint refinement | fully joint |
| **Repair depth** | 1 | 3 | 3 with escalation |

Three of these are genuine algorithmic changes rather than parameter tweaks, and they are where the budget is actually won:

1. **Dropping L0 terms entirely at `fast`.** This removes the Python/GPU round trip from the job, so `fast` is pure Rust and roughly 20× faster than `standard`. The cost is that LM naturalness does not steer the search at all.
2. **Restricting structure to the initiation window at `fast`.** ~5× cheaper on the term that dominates the multi-cistron budget. The cost is that interior structure is unmodelled.
3. **Skipping joint refinement at `fast`.** The mean-field price decomposition alone is used, so the cross-cistron interaction is approximated rather than optimized.

Every one of these is recorded in the lockfile and surfaced in `job_status`, because a `fast` result and a `standard` result are not comparable and users must not be able to confuse them.

### Auto-selection

`budget.tier` defaults to `standard`. At submission the orchestrator checks the detected hardware against the tier's requirements and, if it cannot meet them, **downgrades with a `CF-BUDGET-002` warning rather than failing**. A user on a 4-core laptop asking for `standard` gets `fast` and is told so; they do not get a job that runs for six hours.

## Where the design would need to change for a cluster

Recorded so that scaling out later is not a redesign:

- **Tier-2 instances** are embarrassingly parallel and already independent; distributing them needs only a work queue.
- **The tier-3 evaluation batch** is the natural distribution unit — `E_max` candidates per generation, each independent. The batch-yield protocol already serializes them, so a remote worker pool is a transport change, not an architecture change.
- **The generation loop itself is sequential** and would remain on one node. An island-model variant (several populations with periodic migration) is the standard way to parallelize further and is registered as a v2 option, not v1 — it changes the search's convergence characteristics and would need its own validation.
- **Model training** already assumes a bigger machine.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Cluster-first, workstation as a degraded mode | Excludes most of the intended user base and makes reproduction hard. |
| GPU-required | The GPU only helps the LM term; everything else is CPU-bound integer/float work. Requiring one for a 3 % share of the objective vector is a poor trade. |
| One tier, auto-tuned | Attractive, but a design produced under an auto-chosen configuration is not comparable to another one, and users would have no way to say "spend more". Named tiers are comparable and requestable. |
| Anytime algorithm with a wall-clock budget only | Partly adopted — `max_wall_seconds` and cooperative cancellation already give anytime behaviour. Rejected as the *only* control because the algorithmic choices (which terms exist at all) cannot be made continuously. |

## Consequences

- Every performance claim in the plan is stated against the `workstation` profile; other profiles are scaled from it.
- `fast` tier is a genuinely different objective set (no LM, initiation-window structure), so `ParetoFront.summary.budget_tier` must be checked before comparing two fronts. `TST-BEN-002` asserts the benchmark harness refuses to compare across tiers.
- The structure term dominates at 15 kb, so if the budget is ever missed the first two mitigations are already designed: strided interior sampling, and raising the surrogate's share of structure evaluation.

## Departures

`DEP-BUDG-001` — The architecture proposal does not state a hardware target. This plan fixes it at a single workstation and derives the algorithmic degradation ladder from that constraint, rather than treating performance as a tuning exercise after the fact.
