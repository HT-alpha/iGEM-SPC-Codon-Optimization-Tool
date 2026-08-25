# Acceptance Thresholds and Gates

Satisfies: `BRIEF-2.5-a` — "with pass/fail or acceptance thresholds defined up front".

---

## 1. Why almost every gate is relative, not absolute

An absolute threshold set before the data exists is a guess wearing the authority of a number. "Spearman ≥ 0.5 for codon-level dwell prediction" sounds rigorous, but nobody planning this system knows whether 0.5 is trivially easy or physically unattainable for a given host and dataset — codon-level ribosome density is noisy, and a large share of its variance is protocol artefact rather than biology.

A **relative** gate does not have that problem. "Beats a codon-usage-only baseline by ≥0.05 Spearman with a bootstrap CI excluding zero" is meaningful the moment it is evaluated, on any host, at any noise level, and it tests the claim that actually matters: *does the extra machinery earn its place?*

So the rule for this package is:

> **Gates compare against a named baseline and require a confidence interval excluding zero. Absolute thresholds are used only where the quantity has a definition-driven correct value** — coverage of a 90 % interval should be 90 %; a delta must equal the difference it claims to be; a wall-clock budget is a budget.

Every gate below names its baseline, its statistic, its threshold, its uncertainty treatment and the consequence of failing.

---

## 2. Model-quality gates

### `GATE-EPI-001` — Epistatic surface beats codon usage on held-out dwell

| | |
| --- | --- |
| Statistic | Spearman ρ between predicted and observed per-codon A-site density, held-out genes, split **by gene** |
| Baseline | `codon_usage_only` (relative adaptiveness + position spline, no context terms) |
| Threshold | Δρ ≥ **0.05**, bootstrap 95 % CI (2 000 resamples, resampled by gene) excluding 0 |
| Applies to | every `epistasis_dwell` artifact at promotion |
| On failure | The artifact is not promoted. If it fails for a host repeatedly, that host's `supported_terms` drops `epistasis.dwell_deviation`. |

Resampling **by gene**, not by codon, because codons within a gene are strongly dependent and a codon-level bootstrap gives absurdly narrow intervals.

### `GATE-EPI-002` — Cross-dataset retention

| | |
| --- | --- |
| Statistic | (Δρ on an independent dataset) ÷ (Δρ within the training dataset) |
| Threshold | ≥ **0.50** |
| Applies to | any host with ≥2 independent Ribo-seq datasets |
| On failure | The surface is substantially fitting protocol artefacts. Restrict the artifact to same-pipeline use and mark it prominently. See `EXP-01`. |

### `GATE-BURD-001` — Burden model back-test

| | |
| --- | --- |
| Statistic | Spearman between predicted `g` and measured relative growth cost, on a published burden/expression panel |
| Baseline | `total_protein_output_only` (burden ∝ expression level, ignoring codon composition) |
| Threshold | ρ ≥ **0.50** absolute **and** a positive margin over the baseline with a bootstrap CI excluding 0 |
| On failure | Burden is demoted to `diagnostic` for that capacity tier. See `EXP-04`. |

The absolute component is justified here in a way it is not for `GATE-EPI-001`: growth-rate measurements are far less noisy than codon-level Ribo-seq, and published burden/expression relationships are strong enough that ρ < 0.5 indicates a broken model rather than a hard problem.

### `GATE-UNC-001` — Calibration coverage

| | |
| --- | --- |
| Statistic | Empirical coverage of the nominal 90 % predictive interval on a held-out set |
| Threshold | **[0.87, 0.93]** |
| Applies to | every term declaring `calibrated: true` |
| On failure | `calibrated` is set to `false`, which propagates into probabilistic dominance and automatically depresses affected candidates' robustness scores. Not a build failure — an honest downgrade. |

An **absolute** gate, correctly: a 90 % interval covering 70 % is wrong by definition, with no baseline needed.

### `GATE-KIN-001` — Kinetics term promotion

| | |
| --- | --- |
| Statistic | Coefficient on `pause_concordance` in `solubility ~ concordance + total_expression + (1|family)` |
| Threshold | Significant at α = 0.05 after multiple-comparison correction, correct sign, generalizing to a held-out protein family |
| On failure | Term stays `diagnostic`, or is removed. See `EXP-06`, `AMD-002`. |

---

## 3. Design-quality gates (the comparative benchmarks)

### `GATE-BEN-001` — MFE/CAI frontier vs. LinearDesign / CDSfold

| | |
| --- | --- |
| Setup | Fixed benchmark set of ≥50 proteins. CodonForge runs with exactly two objectives (`structure.mfe_global`, `usage.cai`) to make the comparison apples-to-apples. Competitor outputs enter through the seed-import path so both sides are scored by identical code. |
| Statistic | 2-D hypervolume with a fixed, published reference point |
| Threshold | CodonForge HV ≥ **0.98 ×** the best competitor's HV on ≥90 % of proteins, **and** ≥1 CodonForge point strictly non-dominated by the competitor front on ≥50 % of proteins |
| On failure | Phase 2 does not exit. Not being able to match a specialist two-objective optimizer on its own two objectives means the search or the objective implementation is wrong. |

The `0.98` rather than `1.0` is deliberate: LinearDesign solves a restricted problem with a dynamic program that is exact for its formulation, and a general-purpose MOEA matching an exact DP to within 2 % is the right expectation. Demanding strict dominance everywhere would be demanding that a general tool beat a specialist at the specialist's own game.

### `GATE-BEN-002` — Constraint satisfaction under comparison

| | |
| --- | --- |
| Statistic | Fraction of emitted candidates satisfying every blocking constraint |
| Threshold | **100 %** |
| Note | Reported alongside `GATE-BEN-001` because competitor tools frequently violate the user's assembly constraints, and an MFE/CAI comparison that ignores feasibility is not a fair comparison in either direction. |

### `GATE-BEN-003` — Full-objective hypervolume vs. baselines

| | |
| --- | --- |
| Setup | Full objective set; baselines: CAI-greedy, host-usage sampling, LinearDesign seed, native sequence |
| Statistic | Normalized hypervolume, ≥20 seeds |
| Threshold | Median HV exceeds every baseline with a Wilcoxon signed-rank test at α = 0.01, and the effect size (Cliff's δ) is ≥ 0.5 |
| On failure | Phase 4 does not exit. |

### `GATE-BEN-004` — Ablation attribution

| | |
| --- | --- |
| Setup | Full pipeline vs. the same pipeline with the epistasis layer ablated, and with the burden layer ablated |
| Statistic | Hypervolume difference, and objective-wise difference |
| Threshold | Reported, **not thresholded** |
| Why it exists | A benchmark that only shows "the whole thing beats the baselines" cannot say which layer earned it. Ablation numbers are what make `EXP-10` interpretable, and publishing them is a credibility asset even when a layer contributes little. |

---

## 4. System gates

### `GATE-PERF-001` — Performance budgets

| Case | Tier | Budget |
| --- | --- | --- |
| Single gene ≤2 kb | `fast` | ≤ 30 s |
| Single gene ≤2 kb | `standard` | ≤ 600 s (8 cores + GPU) / ≤ 900 s (8 cores) |
| Multi-cistron ≤15 kb | `fast` | ≤ 120 s |
| Multi-cistron ≤15 kb | `standard` | ≤ 3 600 s (8 cores) |

Measured on the `workstation` reference profile of [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md). Absolute, because a budget is a budget. Failing blocks Phase 4 exit.

### `GATE-PERF-002` — No silent regression

Nightly `TST-PERF-001` fails on a >25 % wall-time regression versus the trailing 7-day median for any phase of either reference job.

### `GATE-DET-001` — Determinism

Bit-identical fronts across 1, 4 and 16 threads, and across two machines of the same architecture, for all three tiers. Absolute; failing blocks any release.

### `GATE-REPRO-001` — Replay fidelity

`codonforge replay --assert-identical` reproduces a stored reference front exactly for every benchmark job. Absolute; failing blocks any release.

---

## 5. Active-learning promotion gates

Evaluated before a recalibrated artifact may be promoted to a channel. Detail in [`../04-api/03-active-learning-ingestion.md`](../04-api/03-active-learning-ingestion.md) §6.

| Gate | Requirement | On failure |
| --- | --- | --- |
| `GATE-ACT-001` | No regression on the tenant's **frozen** holdout beyond a 1-SE tolerance band versus the incumbent | Reject; report which metric regressed |
| `GATE-ACT-002` | No regression on the **public** holdout beyond tolerance | Reject — a tenant adapter must not buy tenant fit with general behaviour |
| `GATE-ACT-003` | Calibration retained: coverage stays in `[0.87, 0.93]` | Reject, or promote with `calibrated: false` |
| `GATE-ACT-004` | Isolation: no foreign-tenant dataset in the manifest | Reject and raise an incident — this is a correctness failure, not a quality one |
| `GATE-ACT-005` | ≥12 observations from ≥2 distinct batches | Defer until more evidence arrives |

`GATE-ACT-005` exists because fitting a capacity curve to a single plate produces a confident, wrong answer, and confident wrong answers are worse than none.

---

## 6. Phase exit criteria

| Phase | Must pass |
| --- | --- |
| 1 | `TST-CORE-*`, `TST-CONS-*`, `TST-SYN-*`, `TST-INGE-*`, `TST-SRCH-002`, `TST-CLI-*`, `GATE-DET-001` |
| 2 | Phase 1 + `TST-STRU-*`, `TST-BURD-*`, `TST-METR-*`, `TST-DEC-*`, `GATE-BEN-001`, `GATE-BEN-002`, `GATE-BURD-001`, `EXP-04` decided |
| 3 | Phase 2 + `TST-EPI-*`, `TST-UNC-*`, `TST-LAND-*`, `GATE-EPI-001`, `GATE-EPI-002`, `GATE-UNC-001`, `EXP-01`/`EXP-08` decided |
| 4 | Phase 3 + `TST-SRCH-008`, `TST-SRCH-009`, `GATE-BEN-003`, `GATE-BEN-004`, `GATE-PERF-001`, `EXP-11` decided |
| 5 | Phase 4 + `TST-ACT-*`, `TST-API-*`, `GATE-ACT-*`, `GATE-REPRO-001`, `EXP-10` running |

`EXP-10` is *running*, not *passed*, at Phase 5 exit — it is a wet-lab experiment whose timeline is set by biology, and blocking a software release on it would be a category error. What Phase 5 must deliver is everything needed to run it.
