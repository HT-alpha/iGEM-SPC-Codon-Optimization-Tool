# Algorithm (d) — Systems-Level Joint Search and the Shared Resource Ledger

Satisfies: `BRIEF-2.3-d`.
Owning modules: `cf-search` (`joint/`), `cf-burden` (the ledger itself).
Open question resolved here: [ADR-0004](../07-decisions/ADR-0004-resource-ledger-granularity.md).

---

## 1. A structural fact that constrains everything below

Because the protein sequences are fixed, **the total amino-acid demand of a construct is invariant under synonymous redesign.** A design cannot make the cell synthesize fewer glutamates. What it *can* do is choose which tRNA isoacceptor decodes each of them.

So the resource ledger has exactly two channels through which codon choice acts, and they behave differently:

| Channel | Controllable? | Mechanism |
| --- | --- | --- |
| **tRNA isoacceptor balance** | Yes, but only by *redistribution* within each synonymous family | Shifting demand from a scarce isoacceptor to an abundant one within the same amino-acid family. Total family demand is fixed. |
| **Ribosome occupancy** | Yes, and this one is reducible in absolute terms | Slower elongation means more ribosome-seconds per protein produced, so a slow CDS ties up more of the ribosome pool for the same output. |

This is worth stating at the top of the spec because it bounds what the burden objective can honestly claim. "Reduce host burden" via codon choice means "stop competing for scarce isoacceptors, and stop holding ribosomes hostage" — not "consume fewer resources". A burden model that appears to reduce total amino-acid demand has a bug.

`DEP-BURD-001` — The architecture proposal describes the burden model as reducing resource consumption. This plan narrows the claim: synonymous redesign can rebalance isoacceptor demand and reduce ribosome occupancy, but cannot reduce total amino-acid demand. The objective is specified accordingly.

---

## 2. The ledger

Let `A` be the set of tRNA isoacceptor pools (anticodon-level; ~40 in *E. coli*). Let `w(x, α) ∈ [0,1]` be the share of codon `x`'s decoding events served by pool `α`, from the wobble table `host_profile.trna.wobble_rules_id`, with `Σ_α w(x, α) = 1`.

For each cistron `k` with codon counts `n_k(x)` and expression rate `f_k` (proteins·cell⁻¹·s⁻¹, from `expected_relative_expression` × plasmid copy number × a scale factor):

```
Demand on pool α:      D_α  =  Σ_k  f_k · Σ_x n_k(x) · w(x, α)
Utilization:           u_α  =  D_α / S_α                       # S_α from HostProfile
Availability:          a_α  =  1 / (1 + u_α)                   # M/M/1-style contention
Effective decode rate: r(x) =  r_max · Σ_α w(x, α) · a_α
Ribosome occupancy:    R    =  Σ_k f_k · Σ_i 1 / r(codon_{k,i})   # ribosome-seconds/s
```

`D_α` includes the **native proteome's own demand**, not just the construct's. Otherwise a construct always looks cheap. For tier-1/2 hosts this comes from the measured proteome; for tier-3 hosts it comes from the genome-wide codon usage weighted by a generic expression prior (`ADR-0004`).

**The ledger is incrementally maintainable.** Substituting codon `x → x'` at one position changes `D` only for the ≤3 pools that decode `x` or `x'`:

```
D_α += f_k · (w(x', α) − w(x, α))     for α ∈ decode(x) ∪ decode(x')
```

That is `O(1)` work, which is why `burden.growth_cost` is an L1-capability term and can participate in tier-2 coordinate ascent. Recomputing `a` and `R` after a batch of edits is `O(|A| + L)`.

---

## 3. Why joint search is required

The objective is **not separable across cistrons**. Cistron 1's codon choices change `u_α`, which changes `a_α`, which changes cistron 2's effective elongation rate and therefore cistron 2's dwell profile, ribosome occupancy and burden. Optimizing each cistron alone and concatenating gives the wrong answer whenever two cistrons compete for the same scarce isoacceptor — which is precisely the case a systems-level tool exists to handle.

Optimizing the concatenated genotype directly is correct but expensive: the search space is `∏ₖ ∏ᵢ Syn(a_{k,i})`, and tier-2's block structure no longer has a clean partition across the cistron boundaries.

---

## 4. The chosen method: price-decomposed alternating optimization

The coupling between cistrons is entirely through the low-dimensional vector `D ∈ ℝ^{|A|}` (≈40 numbers). That is a textbook setting for **dual / mean-field decomposition**: summarize the coupling as a price vector, let each cistron optimize independently against the prices, then update the prices.

**Shadow price of pool `α`:**

```
π_α = ∂ Ψ / ∂ D_α          where Ψ = total system cost
                            (burden growth cost + Σ_k weighted dwell cost)
```

Given the closed forms in §2 and the penalty in [`05-burden-penalty-and-capacity.md`](05-burden-penalty-and-capacity.md), `π_α` is available analytically — no numerical differentiation:

```
∂a_α/∂D_α  =  −1 / (S_α (1 + u_α)²)
π_α        =  g'(u) · (1/S_α)                       # burden channel
            + Σ_k f_k Σ_i  w(cᵢ, α) / r(cᵢ)² · r_max · (−∂a_α/∂D_α)   # elongation channel
```

Each cistron then sees a purely **local, additive** per-codon price term `Σ_α w(x, α) · π_α`, and can be optimized with the unmodified tier-1/2/3 machinery.

### 4.1 Pseudocode

```
JOINT_SEARCH(cistrons K, host, constraints, budget) -> front

  x <- { k: initial assignment for cistron k }
  η <- 0.5                                  # damping
  D_hist <- []

  for outer in 1 .. max_outer (default 6):

      # ---- (1) inner fixed point: self-consistent availability ----------
      a <- previous a  or  1 (first pass)
      for inner in 1 .. 50:
          D    <- LEDGER_DEMAND(x, f, w) + NATIVE_DEMAND(host)
          a_hat<- 1 / (1 + D / S)
          a'   <- (1 − η)·a + η·a_hat
          if max|a' − a| < 1e-4: a <- a'; break
          a <- a'
      π <- SHADOW_PRICES(a, D, S, x, f)

      # ---- (2) independent per-cistron optimization ---------------------
      parallel for k in K:
          if not cistrons[k].is_design_target: continue
          x[k] <- OPTIMIZE_CISTRON(k, prices = π, availability = a,
                                   tiers = [tier1, tier2],
                                   constraints = constraints)

      # ---- (3) convergence / oscillation control ------------------------
      D_hist.append(D)
      if RELATIVE_CHANGE(D_hist[-1], D_hist[-2]) < 1e-3: break
      if OSCILLATING(D_hist, window = 3):
          η <- η / 2
          if η < 0.0625:
              return FALLBACK_FULLY_JOINT(x, ...)   # see §4.3
                                                    # recorded as a warning

  # ---- (4) joint refinement: the mean-field assumption is dropped -------
  P0 <- POPULATION_FROM(x, jitter = block_gibbs)
  return TIER3_NSGA3_JOINT(P0, exact_ledger = true, generations = G_joint)
```

### 4.2 Why step (4) is not optional

Steps (1)–(3) are a *mean-field* approximation: each cistron optimizes against an averaged price rather than against the true, assignment-dependent contention. Mean field is exact only in the limit of many small contributors, and a two-cistron construct on a high-copy plasmid is the opposite of that limit. Step (4) runs the ordinary NSGA-III loop on the **concatenated genotype with the ledger evaluated exactly for every individual** (which costs `O(|A|)` extra per evaluation — nothing). Prices got the search into the right neighbourhood cheaply; the joint refinement makes the final front correct.

`DEP-SRCH-005` — The architecture proposal describes a "systems-level joint search" without specifying how joint optimization stays tractable. This plan uses price decomposition for the bulk of the work and reserves genuinely joint (concatenated-genotype) optimization for a final refinement phase, with an explicit fallback when the decomposition fails to converge.

### 4.3 Failure mode and fallback, stated up front

Damped fixed-point iteration on a supply–demand system can oscillate: cistron 1 abandons a pool, the price collapses, cistron 1 returns, and so on. The design does not assume this away.

- **Detection:** `‖D^t − D^{t−1}‖` fails to decrease over 3 consecutive outer iterations.
- **First response:** halve `η`. Heavier damping converges more slowly but more reliably.
- **Give-up threshold:** `η < 0.0625`.
- **Fallback:** run the fully joint search from the start (concatenated genotype, exact ledger, no prices). Correct, roughly `|K|×` more expensive per generation because the genotype is `|K|×` longer.
- **Visibility:** a `CF-JOINT-001` warning appears in `job_status.warnings`, and `ParetoFront.summary` records that the decomposition was abandoned. A silent fallback would make the performance budget unreproducible.

`TST-SRCH-008` constructs a synthetic two-cistron instance designed to oscillate (two cistrons of equal size competing for a single scarce isoacceptor with a sharp penalty) and asserts that the detector fires and the fallback produces a solution at least as good as the decomposed path.

---

## 5. Complexity

Let `K` = number of cistrons, `L_k` = length of cistron `k`, `L_tot = Σ L_k`, `|A| ≈ 40`.

| Operation | Cost | Note |
| --- | --- | --- |
| Full ledger demand | `O(L_tot · d̄)`, `d̄` ≈ 1.6 pools per codon | ~8 000 flops for a 15 kb construct |
| Incremental ledger update (1 codon) | `O(1)` (≤3 pools) | This is what makes `burden.*` an L1 term |
| Inner fixed point | `O(iters · |A|)` ≤ `50 × 40` | 2 000 flops — free |
| Shadow prices | `O(|A| + L_tot)` | once per outer iteration |
| Per-cistron optimization (tiers 1–2) | `O(sweeps · L_k · b · S · W · M)` | §[`02`](02-block-local-coordinate-ascent.md) |
| **Decomposed path, wall time** | `outer × maxₖ(per-cistron cost)` (cistrons run in parallel) | `6 × ` single-cistron tier-2 |
| **Joint refinement** | `G_joint × N · (L_tot·W·M + |A| + L0 cost)` | dominated by L0 evaluation, as always |
| **Fully joint fallback** | `gens × N · (L_tot·W·M + L0)` | `K×` the per-generation cost of a single cistron |

The decomposition's payoff is in the tier-2 phase: `K` cistrons optimized concurrently against a shared 40-number price vector cost the same wall time as one cistron, whereas the concatenated genotype costs `K×`. For a 5-cistron, 15 kb construct that is the difference between ~30 s and ~150 s of tier-2 time, before the (dominant) L0 costs are counted.

---

## 6. Granularity — what the ledger does and does not track

Fully specified in [ADR-0004](../07-decisions/ADR-0004-resource-ledger-granularity.md). Summary:

| Dimension | v1 decision |
| --- | --- |
| tRNA resolution | **Per-isoacceptor pool (anticodon level)**, with wobble shares. Not per-gene, not per-charged/uncharged. |
| Amino-acid pools | Not tracked. Total amino-acid demand is invariant under synonymous redesign (§1), so an amino-acid ledger would be a constant. |
| Ribosomes | Tracked as aggregate occupancy `R`, one number. |
| Time resolution | **Single steady state.** Two-phase (pre/post-induction) is an opt-in extension gated on `EXP-03`. |
| Energy / ATP / GTP | Not tracked in v1: proportional to total amino-acid demand, hence invariant (§1). |
| Chaperone / secretion capacity | Not tracked. Named as a v2 candidate in the risk register (`RISK-BURD-003`). |

---

## 7. Tests

| Test | Assertion |
| --- | --- |
| `TST-BURD-001` | Ledger conservation: total amino-acid demand is bit-identical before and after 10⁵ random synonymous substitutions (§1). A violation means the wobble table or the demand accumulator is wrong. |
| `TST-BURD-002` | Incremental vs. full ledger agree to `1e-12` over 10⁵ random edits. |
| `TST-SRCH-008` | Oscillation detection and fallback (§4.3). |
| `TST-SRCH-009` | Joint beats independent: on the synthetic two-cistron contention instance of `TST-SYN-003`, joint search must beat independent per-cistron optimization by at least the analytically known margin. This is the test that justifies the whole layer's existence. |

---

## 8. What this build runs

`DEP-SRCH-019` — `design_operon` designs each cistron with the existing single-CDS job, then keeps a synonymous swap only when `burden.growth_cost` of the caller-supplied flux-weighted codon counts falls and the sequence stays feasible. Shadow prices, native-proteome demand, and the joint NSGA-III refinement in §4 are not run. The caller supplies every flux. A default expression rate would be an invented number.
