# Complexity Analysis and Performance Budgets

Satisfies: `BRIEF-2.3-g` — "explicit complexity/performance analysis for each algorithm … with target performance budgets for typical use cases (single gene ≤2 kb, multi-cistron construct ≤15 kb)".
Related: [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md) (target hardware and the degradation ladder), `AMD-001` (windowed vs. global folding).

---

## 1. Symbols and reference values

| Symbol | Meaning | Single gene (2 kb) | Multi-cistron (15 kb) |
| --- | --- | --- | --- |
| `L` | editable codons | 666 | 5 000 (across 5 cistrons) |
| `n` | construct length in nt | ~2 400 (incl. flanks) | ~15 000 |
| `S` | mean synonymous-family size | 2.9 | 2.9 |
| `W` | epistasis window (codons) | 10 | 10 |
| `B` | tier-2 block size (codons) | 30 | 30 |
| `b` | beam width | 64 | 64 |
| `M₁` | L1 (delta-capable) terms | 4 | 4 |
| `M₀` | L0 (oracle-only) terms | 1–2 | 1–2 |
| `N` | population size | 212 | 212 |
| `G` | generations | 300 | 300 |
| `E_max` | exact L0 evaluations per generation | 64 | 64 |
| `K` | cistrons | 1 | 5 |
| `\|A\|` | tRNA isoacceptor pools | ~40 | ~40 |

---

## 2. Per-algorithm complexity, consolidated

| # | Algorithm | Time | Memory | Doc |
| --- | --- | --- | --- | --- |
| a | Potts pseudo-likelihood fit (per L-BFGS iteration) | `O(N_codons · W · S)` | `O(64 + 4096·W)` | [01](01-epistatic-surface-fit.md) |
| a | Dwell-head NB GLM fit | `O(iters · N_rows · (2W+8))` | `O(N_rows)` streamed | [01](01-epistatic-surface-fit.md) |
| a | Score full sequence | `O(L · W · S)` | `O(L·S)` | [01](01-epistatic-surface-fit.md) |
| a | Score single substitution (**L1**) | `O(W · S)` | — | [01](01-epistatic-surface-fit.md) |
| a | Full delta table | `O(L · W · S²)` | `O(L·S)` | [01](01-epistatic-surface-fit.md) |
| b | Tier-1 constrained beam DP | `O(L · b · S · W · M₁)` | `O(b · W)` | [02](02-block-local-coordinate-ascent.md) |
| b | Tier-2, one sweep | `O(L · b · S · W · M₁)` | `O(L·S)` | [02](02-block-local-coordinate-ascent.md) |
| b | Tier-2, one instance | `O(sweeps · L · b · S · W · M₁)` + `sweeps` L0 evals | `O(L·S)` | [02](02-block-local-coordinate-ascent.md) |
| b | Seam reconciliation (fallback path only) | `O(K_blocks · 9 · O_ovl · W · M₁)` | `O(K_blocks)` | [02](02-block-local-coordinate-ascent.md) |
| c | Non-dominated sort | `O(M · N²)` | `O(N²)` worst | [03](03-nsga-loop-and-operators.md) |
| c | NSGA-III normalization + association | `O(M²N + M³ + N·H·M)` | `O(N·M)` | [03](03-nsga-loop-and-operators.md) |
| c | Hypervolume (WFG) | `O(N^{M/2})` worst | `O(N·M)` | [03](03-nsga-loop-and-operators.md) |
| c | One generation | `O(N·L·W·(S+M₁)) + E_max · c_L0` | `O(N·L·S)` | [03](03-nsga-loop-and-operators.md) |
| d | Full ledger | `O(L · d̄)`, `d̄`≈1.6 | `O(\|A\|)` | [04](04-systems-level-joint-search.md) |
| d | Ledger delta (**L1**) | `O(1)` | — | [04](04-systems-level-joint-search.md) |
| d | Availability fixed point | `O(iters · \|A\|)` | `O(\|A\|)` | [04](04-systems-level-joint-search.md) |
| d | Joint outer loop | `outer × maxₖ(cistron cost)` | — | [04](04-systems-level-joint-search.md) |
| e | Burden penalty | `O(\|A\|)` | `O(\|A\|)` | [05](05-burden-penalty-and-capacity.md) |
| f | Incremental violation detection | `O(m_max + k + window/3)` | `O(n)` index | [06](06-constraint-repair.md) |
| f | Repair, depth-`d` beam | `O(d · beam · \|P\| · S · detect)` | `O(beam · d)` | [06](06-constraint-repair.md) |
| f | Feasibility analysis | `O(n · m_max + L · b · S)` | `O(n)` | [06](06-constraint-repair.md) |
| — | Structure: local accessibility | `O(n · w²)`, `w`=80 nt | `O(n · w)` | §3 |
| — | Structure: local accessibility, incremental | `O(w³)` per 3-nt edit | — | §3 |
| — | Structure: global MFE | `O(n³)` time, `O(n²)` memory | | §3 |
| — | Kinetics concordance | `O(L)` | `O(L)` | [07](07-folding-kinetics-objective.md) |

---

## 3. The structure objective is the whole problem — and the fix

Global MFE folding inside an evolutionary loop is not affordable, and pretending otherwise would invalidate every budget below.

| Engine | Complexity | `n = 2 400` | `n = 15 000` |
| --- | --- | --- | --- |
| Zuker MFE (ViennaRNA-class) | `O(n³)` time, `O(n²)` memory | `1.4 × 10¹⁰` ops → **~4 s**, ~90 MB | `3.4 × 10¹²` ops → **~20 min**, ~3.5 GB |
| Beam-search linear folding (LinearFold-class) | `O(n · b_f)` | `2.4 × 10⁵` → **~0.15 s** | `1.5 × 10⁶` → **~0.9 s** |
| **Local windowed accessibility (RNAplfold-class), `w = 80`, span 40** | `O(n · w²)` | `1.5 × 10⁷` → **~13 ms** | `9.6 × 10⁷` → **~85 ms** |
| Local windowed, **incremental** after a 3-nt edit | `O(w³)` (≈`w` windows × `O(w²)` each) | `5 × 10⁵` → **~0.4 ms** | same |

At `N·G = 6.4 × 10⁴` evaluations per job, global MFE costs `6.4 × 10⁴ × 4 s ≈ 71 CPU-hours` for a single 2 kb gene. That is a factor of ~400 outside the budget on an 8-core workstation.

**Decision (`AMD-001`): the default structure objective is windowed local accessibility, not global MFE.** Global MFE is available only at `exact` tier and only for final re-ranking of the reduced front (≤20 candidates × 4 s = 80 s, entirely affordable once).

This is not purely a compute concession. The mechanistic claims the structure objective is serving — accessibility of the ribosome binding region, and local unfolding ahead of the elongating ribosome — are *local* claims. A single global MFE number for a 15 kb polycistronic mRNA is a poor proxy for either. The cheap objective is also the better-targeted one; the trade-off is that a design optimized for local accessibility is not optimized for global thermodynamic stability, which matters for mRNA-therapeutic applications and is why `exact` tier exists.

---

## 4. End-to-end budgets

Assumptions: 8 physical cores, ~3 × 10⁹ effective flops/core/s for this workload class, one consumer GPU optional. LM naturalness costs ~15 ms per 2 kb sequence GPU-batched (~120 ms CPU), scaling ~linearly in length.

### 4.1 Single gene ≤ 2 kb, `standard` tier

| Phase | Derivation | Time (8 cores + GPU) | Time (8 cores, no GPU) |
| --- | --- | --- | --- |
| Ingest, feasibility analysis, model resolution | milliseconds each | 0.5 s | 0.5 s |
| Tier 1 | `4.9 × 10⁶` flops | < 0.1 s | < 0.1 s |
| Tier 2 (200 instances) | `7.8 × 10⁹` flops ÷ 8 cores | 0.4 s | 0.4 s |
| Tier 2, L0 evaluations | 1 600 evals, batched 64 | 1.3 s | 12 s |
| Tier 3, in-core work | `G · N · L · W · (S+M₁)` ≈ `4.4 × 10¹⁰` ÷ 8 | 1.8 s | 1.8 s |
| Tier 3, structure | `G · N × 13 ms ÷ 8` | 103 s | 103 s |
| Tier 3, LM (`E_max`=64) | `G · 64` evals batched | 14 s | 130 s |
| Tier 3, sorting + indicators | `G × ~40 ms` | 12 s | 12 s |
| Repair | worst case `G × 1.5 × 10⁸` ÷ 8 | ≤ 2 s | ≤ 2 s |
| Decision, certification, rationale | 20 candidates, full outputs | 8 s | 15 s |
| **Total** | | **≈ 145 s** | **≈ 280 s** |
| **Budget** | | **≤ 600 s** | **≤ 900 s** |

Headroom is ~4× against the GPU budget and ~3× against the CPU-only budget. That margin is intentional: the estimates assume clean batching and no cache misses, and real workloads are messier.

### 4.2 Multi-cistron ≤ 15 kb, 5 cistrons, `standard` tier

| Phase | Derivation | Time (8 cores + GPU) |
| --- | --- | --- |
| Ingest + feasibility | `O(n · m_max)` | 1 s |
| Tier 1, per cistron, parallel | `5 × 4.9 × 10⁶` flops | < 0.2 s |
| Joint outer loop (6 iterations × tier 2, cistrons in parallel) | `6 × 7.8 × 10⁹` ÷ 8 | 20 s |
| Joint outer loop, L0 evaluations | `6 × 1 600` evals | 25 s |
| Joint refinement tier 3, in-core | `G · N · L_tot · W · (S+M₁)` ≈ `3.3 × 10¹¹` ÷ 8 | 14 s |
| Joint refinement, structure | `G · N × 85 ms ÷ 8` | 676 s |
| Joint refinement, LM | `G · 64 × 120 ms`, batched | 110 s |
| Sorting, indicators, repair | scales with `N`, not `L` | 20 s |
| Decision + certification (incl. optional global MFE on 20 candidates) | `20 × 20 min` if `exact` — **not** run at `standard` | 30 s |
| **Total** | | **≈ 900 s (15 min)** |
| **Budget** | | **≤ 7 200 s (2 h)** on 32 cores; ≤ 3 600 s on 8 cores |

The structure term is **75 % of the multi-cistron budget**. Two mitigations are already in the design and are the first things to reach for if the budget is missed: restrict `local_accessibility` to the initiation windows plus a strided sample of the CDS interior (≈5× cheaper), and raise the surrogate's share of structure evaluation so that only `E_max` individuals per generation get an exact fold.

### 4.3 What each tier changes

Full ladder in [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md). Summary of the compute consequences:

| | `fast` | `standard` | `exact` |
| --- | --- | --- | --- |
| L0 terms | none | LM at `E_max`=64/gen | LM on all + global MFE on the front |
| Structure engine | initiation window only | full-CDS windowed | windowed + global MFE re-rank |
| `N` / `G` | 64 / 80 | 212 / 300 | 320 / 600 |
| Uncertainty | point estimates | ensemble of 5 | ensemble of 5 + conformal |
| 2 kb wall time | **~8 s** | ~145 s | ~35 min |
| 15 kb wall time | **~50 s** | ~15 min | ~4 h |

---

## 5. Memory model

| Component | Expression | 2 kb | 15 kb |
| --- | --- | --- | --- |
| Population genotypes | `N · L` bytes | 141 KB | 1.1 MB |
| Cached field tables | `N · L · S · 4` B | 16 MB | 25 MB |
| Structure window matrices | `threads · n · w · 8` B | 12 MB | 77 MB |
| k-mer index | `O(n)` entries | < 1 MB | < 2 MB |
| Potts artifact | `(64 + 4096·W) · 4` B | 164 KB | 164 KB |
| LM checkpoint (GPU or CPU) | model-dependent | 100–500 MB | 100–500 MB |
| Evaluation cache (LRU, capped) | configurable | ≤ 2 GB | ≤ 2 GB |
| **Working set excluding cache** | | **~0.5 GB** | **~0.7 GB** |

Comfortably inside the 32 GB target of [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md). The only component that could break this is a global MFE fold at 15 kb (`O(n²)` ≈ 3.5 GB), which is a further reason `exact` tier restricts global MFE to the reduced front and runs it one candidate at a time.

---

## 6. Scaling behaviour

| Vary | Effect | Note |
| --- | --- | --- |
| `L` (sequence length) | Tier 1–2 linear; tier 3 in-core linear; structure linear (windowed) or **cubic** (global MFE) | Windowed folding is what keeps the whole system linear in `L`. |
| `N` (population) | Linear in evaluation cost; **quadratic** in non-dominated sorting | Sorting is 3–4 orders of magnitude below evaluation at `N ≤ 500`, so the quadratic term does not bind until `N > ~5 000`. |
| `G` (generations) | Linear | Hypervolume-stagnation termination usually stops well before the cap. |
| `M` (objectives) | Linear in sorting, `O(N^{M/2})` in hypervolume, and **exponential in front size** | The binding constraint is presentational, not computational — see `AMD-008`. |
| `K` (cistrons) | Decomposed path: constant wall time (parallel); joint refinement: linear | This is the payoff of the price decomposition. |
| threads | Near-linear to 16 cores; memory-bandwidth-bound beyond | Tier 2 and evaluation are embarrassingly parallel. |

---

## 7. How these numbers are kept honest

Estimates in a planning document decay. Three mechanisms keep them from silently becoming fiction:

1. **`TST-PERF-001` — regression-performance suite.** A nightly benchmark runs the 2 kb and 15 kb reference jobs at all three tiers and records wall time per phase. A >25 % regression versus the trailing 7-day median fails the build.
2. **Every emitted `ParetoFront` carries `summary.wall_seconds` and per-phase timings**, so real-world deviation from these budgets is observable in the field, not just in CI.
3. **`GATE-PERF-001`** in [`../05-testing/02-acceptance-thresholds.md`](../05-testing/02-acceptance-thresholds.md) makes the budgets in §4 an acceptance criterion for Phase 4 completion, not an aspiration.

The largest modelling uncertainty in §4 is the LM cost, which depends on a model choice not yet made. If the selected LM is 5× more expensive than assumed, `standard` tier for the 15 kb case moves from ~15 min to ~24 min — still inside budget, because the structure term dominates. If the *structure* estimate is wrong by 5×, the 15 kb case breaks its 8-core budget, which is why `AMD-001` and the `E_max`-rationing of structure evaluation are the designated first mitigations.
