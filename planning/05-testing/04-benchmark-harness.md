# Benchmark Harness

Satisfies: the benchmarking half of `BRIEF-2.5`.
Owning module: `cf-bench` (test-only; never a runtime dependency of the shipped library).
Thresholds: [`02-acceptance-thresholds.md`](02-acceptance-thresholds.md) §3.

---

## 1. Design principles

**1. Benchmarks run through the product's own code paths.** Competitor outputs enter via the seed-import path ([ADR-0007](../07-decisions/ADR-0007-external-tool-seed-import.md)) and are scored via `/v1/score`. There is no separate benchmark scoring implementation, because a benchmark that differs from the product is how benchmark claims become wrong without anyone noticing.

**2. Every result is a lockfile plus a front.** A benchmark result that cannot be replayed is an anecdote. `TST-BEN-001` asserts that stored reference results reproduce exactly from their lockfiles.

**3. Comparisons are refused when they would be invalid.** `TST-BEN-002`: the harness will not compare fronts computed at different budget tiers, with different lockfiles, or against different reference points. Hypervolume is not comparable across reference points, and "our hypervolume beat theirs" is otherwise unfalsifiable.

**4. Ablations are first-class.** Every suite reports the full pipeline *and* single-layer ablations. Without them a benchmark can say the system works but not which part earned it — and "which part earned it" is the entire scientific claim.

**5. Baselines include the honest ones.** The native/wild-type sequence and the plain CAI-greedy design are always present. A tool that cannot beat the wild-type sequence should find that out internally, first.

---

## 2. Datasets

| Id | Contents | Used by | Provenance |
| --- | --- | --- | --- |
| `bm_proteins_50` | 50 heterologous proteins, 200–800 codons, spanning GC content, length and hydrophobicity | `bench_mfe_cai_frontier`, `bench_hypervolume_vs_baselines` | Public sequences; manifest with per-record digests. **Not in this build.** The panel that is present is 50 endogenous CDS from `NC_000913.3` (`DEP-BEN-002`). `GATE-BEN-001` is not claimed, because LinearFold's global MFE is not `structure.local_accessibility`. |
| `bm_constructs_20` | 20 multi-cistron constructs, 4–15 kb, 2–6 cistrons | `bench_joint_search`, `perf_*` | Registry parts + published operons |
| `bm_riboseq_ecoli` | ≥2 independent *E. coli* Ribo-seq datasets, reprocessed identically | `bench_epistasis_heldout`, `EXP-01` | GEO/SRA accessions in the manifest |
| `bm_riboseq_yeast` | Same for *S. cerevisiae* | as above | as above |
| `bm_burden_panel` | Published burden/expression panel: expression + growth rate per construct | `bench_burden_backtest`, `EXP-04` | Codon-Harmony-Index-class dataset |
| `bm_solubility_panel` | Synonymous variants with soluble-fraction measurements, ≥3 families | `EXP-06` | Literature compilation |
| `bm_reference_jobs` | The 2 kb and 15 kb reference jobs of [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md) | `TST-PERF-001` | Fixed inputs, fixed seeds |

Each dataset has a manifest with a content digest and per-record provenance. A benchmark result records the dataset digest, so "which version of the benchmark set" is never ambiguous.

---

## 3. Suites

### `bench_mfe_cai_frontier` — head-to-head against LinearDesign and CDSfold

The comparison `BRIEF-2.5` names explicitly, and the one most likely to be read adversarially. Design choices that make it fair in both directions:

- **Restrict DNA Design Suite to two objectives** (`structure.mfe_global`, `usage.cai`). Comparing a 5-objective front against a 2-objective one on 2 axes is not a comparison.
- **Sweep the competitor's trade-off parameter.** LinearDesign's `λ` traces its own frontier; comparing against a single `λ` point is a straw man. The sweep is recorded in `SeedSet.seeds[].source.parameters`.
- **Score both sides with the same code**, via seed import and `/v1/score`.
- **Report constraint feasibility alongside** (`GATE-BEN-002`). Competitor tools frequently violate the user's assembly constraints, and an MFE/CAI comparison that ignores feasibility flatters DNA Design Suite unfairly in one direction and the competitors unfairly in the other. Both numbers are published.
- **Use the same folding engine for both sides**, recorded in the lockfile. A ViennaRNA MFE and a native-kernel MFE are not the same number, and a benchmark that mixes them is measuring the engine.

Threshold: `GATE-BEN-001`.

**Licence note.** LinearDesign and LinearFold are not redistributable (`RISK-DEP-001`, `AMD-009`), so this suite runs in a **separate CI job against user-installed binaries**, or against published output where binaries are unavailable. Its results are published; the binaries are never vendored.

### `bench_hypervolume_vs_baselines` — full objective set

Baselines: CAI-greedy, host-usage sampling, LinearDesign seed, native sequence, and random-feasible. ≥20 seeds per protein. Reports normalized hypervolume with a published reference point derived from each term's `reference_range`, plus per-objective marginals — because a hypervolume win concentrated in one objective means something different from a broad win, and only the marginals show which happened.

Threshold: `GATE-BEN-003`.

### `bench_ablation` — attribution

| Arm | Configuration |
| --- | --- |
| `full` | All layers |
| `no_epistasis` | `epistasis.*` removed from the objective set |
| `no_burden` | `burden.*` removed |
| `no_kinetics` | `kinetics.*` removed |
| `no_joint` | `search.joint_search = false` (multi-cistron only) |
| `constraints_only` | Phase 1 configuration: feasibility + CAI only |

Reported, not thresholded (`GATE-BEN-004`). This is the arm set that makes `EXP-10` interpretable, and publishing it — including where a layer contributes little — is a credibility asset rather than a liability.

### `bench_epistasis_heldout` — model quality

Held-out dwell prediction with by-gene splits and by-gene bootstrap. Reports within-dataset and cross-dataset Δρ against the codon-usage baseline. Thresholds: `GATE-EPI-001`, `GATE-EPI-002`.

### `bench_burden_backtest` — burden model quality

Predicted `g` against measured relative growth cost on `bm_burden_panel`, at each available capacity tier, with the tier-3 path run with abundance data deliberately withheld. Threshold: `GATE-BURD-001`. Feeds `EXP-04`.

### `bench_joint_search` — systems-level attribution

Joint vs. independent per-cistron optimization on `bm_constructs_20` and on the synthetic contention instance `SYN-03`. Reports the hypervolume margin, the decomposition's outer-iteration count, and how often the fully-joint fallback fired. Thresholds: `TST-SRCH-009` (synthetic, exact margin known); reported-only on real constructs.

### `perf_reference_jobs` — performance regression

The 2 kb and 15 kb reference jobs at all three tiers, timed per phase, nightly, on fixed hardware. Fails on >25 % regression versus the trailing 7-day median (`GATE-PERF-002`) or on any budget breach (`GATE-PERF-001`).

---

## 4. Statistical practice

Small choices that decide whether the numbers survive review:

| Practice | Rule |
| --- | --- |
| Seeds | ≥20 per configuration; report median and IQR, never a single run |
| Paired comparison | Same seeds across arms, so tests are paired (Wilcoxon signed-rank) and seed variance does not swamp the effect |
| Effect size | Cliff's δ alongside every p-value; a significant but negligible effect is reported as such |
| Multiple comparisons | Benjamini–Hochberg across the protein set within a suite |
| Resampling unit | By gene for Ribo-seq, by protein for design benchmarks, by batch for wet-lab panels. Never by codon or by replicate. |
| Reference point | Fixed per objective set, derived from `reference_range`, published in the result document |
| Blinding | For any suite feeding a wet-lab experiment, arm labels are withheld from the assay operator |

---

## 5. Result format and publication

Every benchmark run emits:

```
benchmarks/results/<run_id>/
├── manifest.json        # datasets + digests, code version, hardware, seeds
├── lockfiles/           # one per job, so every number is replayable
├── fronts/              # the ParetoFront documents
├── metrics.parquet      # tidy: (suite, arm, protein, seed, metric, value)
├── gates.json           # every gate evaluated, with value and threshold
└── report.html          # rendered comparison, generated from the above
```

Nothing in `report.html` is computed during rendering. Every number comes from the tidy metric table, which comes from replayable jobs. A reviewer can recompute the entire report from the raw artifacts, which is the property that makes the benchmark claims checkable rather than merely stated.

`DEP-BEN-001`. This build writes that table as `metrics.json`, not `metrics.parquet`. The reference point is the nadir end of each published range and is not derived from the front under comparison. A comparison across budget tier, lockfile digest or reference point is refused (`TST-BEN-002`). External baselines run a user-installed executable whose argv the user supplied; the binary is not vendored (`AMD-013`).
