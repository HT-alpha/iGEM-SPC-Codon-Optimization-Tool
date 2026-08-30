# Algorithm (c) — The NSGA-II/III Loop and Synonymous-Respecting Operators

Satisfies: `BRIEF-2.3-c`.
Owning modules: `cf-search` (loop, operators), `cf-metrics` (sorting, reference directions, indicators).
Related: [ADR-0005](../07-decisions/ADR-0005-uncertainty-and-pareto-ranking.md) (uncertainty in ranking), `AMD-008` (many-objective front explosion).

---

## 1. The property that makes this easy, and the one that does not

**The synonymous design space is a product space.** An individual is an assignment `x ∈ ∏ᵢ Syn(aᵢ)`. Every position is chosen independently, and every choice preserves the protein. Two consequences:

- **Every mutation and every crossover at codon boundaries produces a valid individual.** There is no repair, no penalty and no rejection sampling needed for the protein constraint — the one constraint that matters most. Most MOEA applications to sequence design spend their complexity budget here; CodonForge gets it free by choosing the right genotype.
- **What *can* break is the hard sequence constraints** (restriction sites, repeats, GC windows). These are handled by constraint-domination plus the repair operator of [`06-constraint-repair.md`](06-constraint-repair.md), never by discarding individuals.

Genotype representation: `Vec<u8>` of codon indices, length `L`, plus a cached 64-bit hash for duplicate suppression and a cached field table for `O(W·S)` deltas.

---

## 2. NSGA-II or NSGA-III?

| Objectives `M` | Algorithm | Why |
| --- | --- | --- |
| `M ≤ 3` | **NSGA-II** | Crowding distance is a good, cheap diversity proxy in low dimension, and NSGA-II has no reference-direction bookkeeping. |
| `M ≥ 4` | **NSGA-III** (default) | Crowding distance degrades badly as `M` grows: nearly every individual becomes non-dominated, so selection pressure comes almost entirely from the diversity operator, and crowding distance is a poor one in high dimension. Reference-direction niching supplies structured, predictable coverage instead. |

CodonForge's default objective set (`epistasis.dwell_deviation`, `structure.local_accessibility`, `burden.growth_cost`, `usage.cai`, optionally `kinetics.pause_concordance` and `lm.naturalness`) is `M = 4…6`, so **NSGA-III is the default**, and the plan is written around it.

`DEP-SRCH-003` — The architecture proposal writes "NSGA-II/III" without choosing. This plan chooses NSGA-III for `M ≥ 4` and keeps NSGA-II only for `M ≤ 3`, with the switch automatic and recorded in the job status.

`M ≥ 4` also triggers the front-size problem addressed by `AMD-008`: the *search* handles many objectives fine, but the *presentation* does not, which is why `job_request.objectives[].role` allows demoting a term to an `epsilon_constraint` or a `diagnostic`.

---

## 3. Reference directions and normalization

Reference directions come from Das–Dennis with `p` divisions, `H = C(M+p−1, p)`, `p` chosen so `H` is the smallest value `≥ N_pop`; a two-layer construction (boundary + inner layer at 0.5 scale) is used for `M ≥ 6`, where single-layer Das–Dennis puts almost all directions on the simplex boundary.

| `M` | `p` | `H` | `N_pop` default |
| --- | --- | --- | --- |
| 3 | 12 | 91 | 92 |
| 4 | 7 | 120 | 120 |
| 5 | 6 | 210 | 212 |
| 6 | 4 + 2 (two-layer) | 126 + 21 = 147 | 148 |

Objectives are normalized every generation by the standard NSGA-III procedure: update the ideal point, find extreme points by minimizing an achievement scalarizing function per axis, solve for the hyperplane intercepts, and translate/scale. Two guards against the well-known degeneracies:

- If the extreme-point matrix is singular or an intercept is non-positive, fall back to the per-axis max of the current front (a degenerate front, e.g. when two objectives are perfectly correlated, otherwise produces a division by ~0).
- Terms declared `maximize` are negated once, at ingress to `cf-metrics`, so all internal math is minimization. The sign convention lives in exactly one place.

**Niching selection is deterministic.** `DEP-SRCH-009`. Deb's algorithm makes two random choices: which least-occupied reference point to fill when several are tied, and which associated member to take when the chosen niche is already occupied. Both are replaced by deterministic rules — the lowest reference index, and the minimum perpendicular distance. Determinism by construction rather than by seeding: `GATE-DET-001` requires bit-identical output across thread counts, and a rule with no randomness in it cannot depend on scheduling, whereas a correctly-seeded one still has to be proved not to. It also avoids adding two sub-seeding sites under tier 3's reserved root. The cost is a mild preference for members lying closest to a direction, where Deb spreads randomly among a niche's members; that trades a little diversity for a guarantee the release gate treats as absolute.

---

## 4. Ranking: constraint-domination with an explicit `unknown`

```
COMPARE(u, v) -> {u better, v better, incomparable}

  if u.feasible and not v.feasible:  return u better
  if v.feasible and not u.feasible:  return v better

  if not u.feasible and not v.feasible:
      # normalized total violation, lower is better
      return compare(CV(u), CV(v))          # scalar, so never incomparable

  # both feasible
  if u has any UNKNOWN term or v has any UNKNOWN term:
      return incomparable                    # never impute; see below
  return PARETO_COMPARE(u.f, v.f)
```

`CV(x) = Σ_c w_c · normalized_excess_c(x)` over blocking constraints, with each constraint's excess normalized by its own scale so that "one extra BsaI site" and "GC 3 % over the window limit" are commensurable.

**Only *actionable* violations count.** `DEP-SRCH-008`, reconciling this section with `DEP-CONS-002`/`DEP-CONS-004`. A violation lying wholly inside immutable sequence — a site in the vector backbone, a repeat spanning only locked codons — cannot be removed by any synonymous assignment, so charging it to the individual is charging it for something it cannot change. At best every individual carries the same constant and the term is noise; at worst the whole population is classed infeasible forever, and since constraint-domination makes any feasible individual beat any infeasible one, the entire ranking is then driven by a quantity no individual can influence. `ConstraintChecker` already reports actionability per violation, and every other layer — repair, the feasibility analyser, tier 1's pruning and post-check, tier 2's accept step — draws the line the same way.

**Infeasible individuals are kept, not discarded.** Deb's constraint-domination makes any feasible individual beat any infeasible one, so feasibility is reached under selection pressure, while infeasible individuals survive as stepping stones across a constraint barrier. Discarding them is the classic way to make a tightly-constrained codon design problem return nothing.

**`UNKNOWN` is a first-class value, never imputed.** When a provider returns `status: partial` (missed its latency budget), the missing term is `UNKNOWN`. Such an individual is placed at the back of its feasibility class and put at the head of the next generation's exact-evaluation queue. Imputing a mean or a worst-case value would bias the front in a way that is invisible in the output.

Non-dominated sorting is the standard fast method: `O(M · N²)`. At `M=5, N=212` that is `≈2.2 × 10⁵` comparisons per generation — three to four orders of magnitude below evaluation cost, so no effort is spent on the asymptotically better `O(N logᴹ⁻¹ N)` algorithms.

---

## 5. Variation operators

All operators respect `locked_codons` unconditionally: locked positions are excluded from every index set before sampling.

### 5.1 Mutation

| Operator | Description | Default weight | Cost |
| --- | --- | --- | --- |
| `uniform_synonymous` | Each editable position mutates independently with `p_m = 1/L_edit`; the replacement is drawn from `Syn(aᵢ)\{cᵢ}` weighted by host relative adaptiveness. | 0.35 | `O(L)` |
| `potts_guided` | Same positions, but the replacement is drawn from the Potts conditional `P(x | context)` at temperature `T = 1.5`. Produces context-plausible moves. | 0.35 | `O(L·W·S)` |
| `block_gibbs` | Pick a start uniformly and a length `~ Geometric(1/12)`; do one Gibbs sweep through the block sampling each codon from its Potts conditional. A large, coherent move. | 0.20 | `O(len·W·S)` |
| `violation_directed` | For infeasible individuals only: mutate positions inside a current violation span. Lets an individual escape infeasibility on its own instead of waiting for the repair operator. | 0.10 | `O(span)` |

`potts_guided` is deliberately not the sole mutation operator. It samples from the host's context distribution, so used alone it drags the population toward host-typical sequence and erodes exactly the diversity the multi-objective search needs. `T = 1.5` flattens the conditional for the same reason.

### 5.2 Crossover — block-uniform by default

```
BLOCK_UNIFORM_CROSSOVER(parent₁, parent₂, mean_block = 30) -> child₁, child₂
  i <- 0
  take_from_1 <- coin_flip()
  while i < L:
      len <- 1 + Geometric(1 / mean_block)
      j   <- min(i + len, L)
      child₁[i:j] <- (take_from_1 ? parent₁ : parent₂)[i:j]
      child₂[i:j] <- (take_from_1 ? parent₂ : parent₁)[i:j]
      take_from_1 <- not take_from_1
      i <- j
  # locked positions are restored from the parents afterwards; they are
  # identical in both parents by construction, so this is a no-op in practice
  # and exists only as a defensive invariant check.
```

**Why block-uniform and not per-codon uniform.** The entire premise of the epistatic layer is that adjacent codons interact within a `W ≈ 10` window. Per-codon uniform crossover breaks, in expectation, half of all adjacent pairs — it systematically destroys the co-adapted context that tier 2 spent its whole budget building. Drawing block lengths with mean ≈ `B` keeps most windows intact while still mixing. One-point and two-point crossover remain available (`search.crossover` config) and are the right choice when `W = 1`, i.e. when the epistatic layer is disabled.

`DEP-SRCH-004` — The architecture asks for "synonymous-codon-respecting mutation/crossover operators". This plan observes that *every* codon-boundary operator is synonymous-respecting because the space is a product space, and therefore spends the design effort on a different question — preserving epistatic context — which is where operator choice actually matters.

### 5.3 Duplicate suppression and seed injection

Offspring are hashed on their codon array; exact duplicates of any current population member or of an already-generated sibling are rejected and re-generated (up to 8 attempts, then force an extra mutation). Duplicate-free populations matter more here than in continuous MOEAs because the space is finite and tier 2 deliberately converges many instances to similar local optima.

Imported seeds ([ADR-0007](../07-decisions/ADR-0007-external-tool-seed-import.md)) are injected at initialization only, capped at `seed_set.injection_policy.max_fraction_of_population` (default 0.2), each accompanied by two mutated copies so crossover has material immediately.

---

## 6. Controlling the cost of L0 terms: surrogate-assisted evaluation

Every offspring's L1 terms are computed exactly in-core — they are nearly free. The L0 terms (LM naturalness, global MFE at `exact` tier) are the entire cost of tier 3, so they are rationed:

```
EVALUATE_GENERATION(offspring, surrogate ĝ, budget E_max) -> evaluated offspring

  for x in offspring:
      x.f_L1 <- EXACT_L1(x)                       # in-core, exact
      x.f_L0_hat <- ĝ(features(x))                # ridge regression, ~microseconds

  chosen <- NONDOMINATED(offspring using f_L1 ∪ f_L0_hat)      # (a)
          ∪ RANDOM_SAMPLE(offspring, 0.10·|offspring|)         # (b)
          ∪ TOP_K_BY_EXPECTED_HV_GAIN(offspring, k)            # (c)
          ∪ DEFERRED_FROM_LAST_GENERATION                      # (d)
  chosen <- TRUNCATE(chosen, E_max)                            # default 64

  results <- WORKER.SCORE_BATCH(chosen)           # ONE rust->python crossing
  for x in offspring \ chosen:
      x.f_L0 <- UNKNOWN                           # not imputed
  refit ĝ every G_s = 10 generations on all (x, true f_L0) pairs seen
```

Component (b) — the uniform random sample — is not optional: without it the surrogate is only ever trained on points the surrogate already likes, and it drifts. It is the same failure mode as the active-learning feedback loop in `AMD-004`, at a smaller scale.

**Certification.** No candidate may appear in a `ParetoFront` with a surrogate-derived objective. After termination, `cf-search` exactly evaluates every member of the reduced front and re-sorts; if certification changes the front, the certified front is what is emitted and `summary.termination_reason` is unchanged but a `warnings` entry records the correction. `TST-SRCH-005` asserts that every emitted candidate's provenance shows an exact evaluation.

---

## 7. The generation loop

```
TIER3_NSGA3(P₀, ref_dirs, constraints, budget) -> front

  P <- P₀                                          # from tier 2, size N
  EVALUATE_GENERATION(P, ĝ, E_max)
  hv_history <- []

  for gen in 1 .. max_generations:
      Q <- []
      while |Q| < N:
          p₁, p₂ <- TOURNAMENT_SELECT(P, k=2, comparator=COMPARE)
          c₁, c₂ <- BLOCK_UNIFORM_CROSSOVER(p₁, p₂)   with prob p_c = 0.9
          c₁ <- MUTATE(c₁);  c₂ <- MUTATE(c₂)
          c₁ <- REPAIR(c₁, constraints, delta_table)  # see 06-constraint-repair
          c₂ <- REPAIR(c₂, constraints, delta_table)
          Q <- Q ∪ DEDUP(c₁, c₂)

      EVALUATE_GENERATION(Q, ĝ, E_max)
      R <- P ∪ Q
      fronts <- FAST_NONDOMINATED_SORT(R, COMPARE)
      P <- SELECT_BY_REFERENCE_DIRECTION_NICHING(fronts, ref_dirs, N)

      hv_history.append(HYPERVOLUME(P, reference_point))
      TABU_TICK(constraints)                       # decays repair tabu tenure

      if STAGNATED(hv_history, window=30, rel_eps=1e-3): break
      if OVER_BUDGET(): break

  return CERTIFY(P)
```

Termination is on hypervolume stagnation by default because it is the only signal that tracks the thing the user cares about. Generation and wall-clock caps are backstops. `hv_history` is streamed into `job_status.current_front_summary.hypervolume_history` so a UI can show it live.

---

## 8. Complexity per generation

| Step | Cost | At `M=5, N=212, L=666` |
| --- | --- | --- |
| Variation (crossover + mutation) | `O(N · L · W · S)` for `potts_guided` | `≈ 4.1 × 10⁶` flops |
| Repair | `O(N · V · beam · depth · W)`, `V` = violations per individual (typically 0–2) | `≈ 10⁵–10⁸` flops (typical–worst) |
| Exact L1 evaluation | `O(N · L · W · M₁)` | `≈ 5.6 × 10⁶` flops |
| Surrogate L0 prediction | `O(N · d)` , `d` ≈ 200 features | `≈ 4 × 10⁴` flops |
| Exact L0 evaluation | `E_max` × per-item cost | 64 × ~15 ms = **≈ 1 s** (batched, GPU) |
| Non-dominated sort | `O(M · N²)` | `≈ 2.2 × 10⁵` comparisons |
| Normalization + association | `O(M²N + M³ + N·H·M)` | `≈ 2.3 × 10⁵` flops |
| Hypervolume (WFG, incremental) | `O(N^{M/2})` worst case; front is capped at `N` and `M ≤ 6` | ~10–40 ms |

**Conclusion: exact L0 evaluation is ~95 % of tier-3 wall time.** Every optimization that matters is therefore about reducing `E_max` or making an L0 term cheaper — which is what `AMD-001` (windowed instead of global folding) and the surrogate rationing above are for. Optimizing the sorting or the operators would be wasted effort.

---

## 9. Tests

| Test | Assertion |
| --- | --- |
| `TST-METR-001` | Non-dominated sorting and hypervolume agree with a brute-force implementation on 10⁴ random small fronts. |
| `TST-METR-002` | Das–Dennis directions sum to 1, are unique, and their count equals `C(M+p−1,p)`; two-layer construction covers the interior. |
| `TST-SRCH-005` | Every candidate in an emitted front carries an exact (non-surrogate) evaluation. |
| `TST-SRCH-006` | Protein invariance: after 10⁶ randomized operator applications, every individual still translates to the input protein and every locked span is unchanged. |
| `TST-SRCH-007` | Constraint-domination ordering matches Deb's definition on an enumerated truth table including `UNKNOWN` cases. |
