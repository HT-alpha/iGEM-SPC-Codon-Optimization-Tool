# Algorithm (b) — Three-Tier Search, Block-Local Coordinate Ascent, and Seam Reconciliation

Satisfies: `BRIEF-2.3-b`.
Owning module: `cf-search`.
Related: [ADR-0002](../07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md) (why exact deltas, not gradients), [`08-complexity-and-performance-budgets.md`](08-complexity-and-performance-budgets.md).

---

## 1. The three tiers and what each is for

| Tier | Role | Method | Cost class | Output |
| --- | --- | --- | --- | --- |
| **1 — seed** | Get from "nothing" to "feasible and decent" fast | Beam-search DP over the decomposable objective, on the product with a constraint automaton | `O(L · b · S)` | A handful of feasible starting points |
| **2 — block CA** | Exploit the decomposable structure hard, cheaply, and produce a *diverse* population | Red-black block coordinate ascent, one instance per reference direction | `O(sweeps · L · b · S · W)` per instance | `N_pop` locally-optimal, spread-out individuals |
| **3 — MOEA** | Explore the true multi-objective front, including the expensive non-decomposable terms | NSGA-III ([`03-nsga-loop-and-operators.md`](03-nsga-loop-and-operators.md)) | `O(gens · N_pop · full evaluation)` | The Pareto front |

The division of labour is dictated by [ADR-0002](../07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md): tiers 1 and 2 touch only L1-capable (delta-providing) terms and therefore never leave the Rust core; tier 3 is the only tier that pays for L0 terms, and it pays once per generation via the batch-yield protocol.

---

## 2. Tier 1 — constrained beam DP seeding

The decomposable part of the objective is a chain with interaction span `W`. Exact DP over that chain needs state = the last `W` codons, i.e. `|S|^W` states — `6¹⁰ ≈ 6 × 10⁷` at the default window. Exact DP is therefore used only when `W ≤ 3` (`≤216` states, used by the synthetic test suite where it provides ground truth); otherwise beam search.

Hard constraints that are *regular languages* — forbidden motifs, restriction sites, homopolymer runs, `custom_regex` — are compiled into an Aho-Corasick automaton and carried in the DP state as a product. A state that would enter an accepting (violating) state is simply never generated, so tier 1 output is feasible with respect to those constraints **by construction**.

Constraints that are not regular — `kmer_repeat` (a counting constraint over the whole construct) and `gc_window` — are *not* enforced in the DP. They are checked after and handed to the repair operator. This is a deliberate scope limit: encoding a global occurrence count in the DP state would multiply the state space by the count bound for every monitored k-mer.

```
TIER1_SEED(protein, host, constraints, beam b, weights λ) -> [seeds]

  A <- AHO_CORASICK(regular constraints)          # motif DFA
  beam <- [ (dfa_state = A.start, prefix = "", ctx = [], score = 0) ]

  for i in 0 .. L-1:
      next <- []
      for st in beam:
          for x in Syn(a_i):
              if i in locked: x <- locked_codon(i)
              s' <- A.step(st.dfa_state, x)
              if A.is_violating(s'): continue      # infeasible, prune
              Δ <- Σ_m λ_m · delta_m(x | st.ctx)   # O(W·M), exact, L1 terms only
              next.append( state(s', st.prefix+x, shift(st.ctx, x), st.score + Δ) )
      beam <- TOP_B_BY_SCORE_WITH_DIVERSITY(next, b)
  return DISTINCT(beam)
```

`TOP_B_BY_SCORE_WITH_DIVERSITY` keeps the best `b` states but caps how many may share the same trailing `W`-codon context, which prevents the beam collapsing onto one lineage — the standard beam-search failure.

**Cost.** `O(L · b · S · W · M)`. At `L=666, b=64, S=2.9, W=10, M=4`: `≈ 4.9 × 10⁶` flops, i.e. a few milliseconds. Tier 1 is free relative to everything else.

---

## 3. Tier 2 — red-black block coordinate ascent

### 3.1 The partition, and why it removes the reconciliation problem

Partition the **codon positions** into contiguous blocks of size `B`, and 2-colour them:

```
 codon   0        B       2B       3B       4B       5B
       [0,B)   [B,2B)  [2B,3B) [3B,4B) [4B,5B) [5B,6B)
 colour  RED     BLACK    RED     BLACK    RED     BLACK
```

`DEP-SRCH-007` — Codon positions, not editable positions. An earlier wording of this section said "editable positions" while the diagram showed codon ranges, and the two differ on a sequence with locked codons. Partitioning codon positions makes the independence argument direct: two positions in different same-colour blocks are then at least `B` codons apart by construction, with no reasoning about how many locked positions happen to lie between them. A block containing only locked codons is a no-op, which costs nothing. The alternative reading is not *wrong* — `B` editable positions apart implies at least `B` codons apart — but it makes the argument indirect for no gain, and the implementation follows the diagram.

**Choose `B ≥ R`,** where `R` is the largest interaction radius among the scored terms. Then any two positions in *different red blocks* are at least `B ≥ R` codons apart, so no coupling term connects them. `W` is the right value for `R` only when the epistatic window is the widest coupling in play; §4 exists precisely because a term can reach further, so the condition is stated against what the terms declare rather than against `W`. `TermProvider::interaction_radius` is what they declare it with, and tier 2 refuses a term whose radius exceeds `B` rather than assuming independence it does not have. Holding all black blocks fixed, the red blocks are conditionally independent: they can be optimized **simultaneously, in parallel, exactly, with no interaction and therefore nothing to reconcile.** Then swap colours and repeat.

This is the central design decision of tier 2, and it is a departure from the architecture proposal:

`DEP-SRCH-001` — The architecture proposal specifies overlapping blocks plus an overlap-reconciliation step. With a red-black partition and `B ≥ W`, overlap is unnecessary and reconciliation is vacuous: each half-sweep is an exact block coordinate-ascent step on a set of conditionally independent variables. The overlapping-block scheme is retained only for the residual case in §4, where it genuinely cannot be avoided.

**Monotonicity.** A half-sweep is a coordinate-ascent step over a partition of the variables. If each block update never returns something worse than the incumbent, the scalarized objective is non-decreasing across the whole sweep, and since the objective is bounded above, tier 2 converges. Block updates are made non-worsening by *seeding the beam with the incumbent assignment* (§3.3) — a one-line implementation detail that buys a convergence guarantee.

### 3.2 Multi-objective handling: one instance per reference direction

Coordinate ascent is single-objective, but tier 2's job is to seed a *diverse* tier-3 population. Run `N_pop` independent instances, instance `r` scalarizing with the augmented Tchebycheff achievement scalarizing function against reference direction `λ^(r)` (Das–Dennis, shared with NSGA-III via `cf-metrics`):

```
ASF(f; λ, z*, z^nad) = max_m [ λ_m · (f_m − z*_m) / (z^nad_m − z*_m) ]
                     + ρ · Σ_m  λ_m · (f_m − z*_m) / (z^nad_m − z*_m),     ρ = 1e-4
```

`ASF` is not additive over positions, so it cannot drive a chain DP directly. Tier 2 therefore uses a two-step accept:

1. **Linearize.** At the incumbent, the active `max` term identifies a weight vector; take the weighted-sum linearization `Σ_m w_m f_m` with `w` from the active term plus the `ρ` augmentation. This *is* additive, so the DP works.
2. **Accept on the truth.** Evaluate the true `ASF` of the block proposal; accept it only if it does not worsen. Because the incumbent is in the beam, "do nothing" is always available, so acceptance is always possible.

Monotonicity therefore holds in the *true* scalarized objective, not merely in the linearization.

### 3.3 Block update

```
UPDATE_BLOCK(block [p, q), fixed context, weights w, beam b) -> assignment

  # state = trailing min(W, q-p) codons of the partial assignment
  beam <- [ (ctx = left_context(p), assign = [], score = 0) ]
  incumbent_path <- current assignment on [p, q)

  for i in p .. q-1:
      next <- []
      for st in beam:
          for x in Syn(a_i):
              if i is locked: x <- locked codon
              if MOTIF_DFA_WOULD_VIOLATE(st, x): continue
              Δ <- Σ_m w_m · delta_m(x | st.ctx)      # O(W·M) exact deltas
              next.append(extend(st, x, Δ))
      # -- the line that buys monotonicity --
      next <- next ∪ { incumbent_path_prefix(i) }
      beam <- TOP_B(next, b)

  best <- argmax over beam of (score + right_boundary_correction)
  return best.assign
```

`right_boundary_correction` accounts for couplings from the block's last `W` codons into the fixed context on the right; without it the DP is myopic at the block's trailing edge and produces a systematic seam artefact.

After the DP, a **1-opt polish** pass over the block (try every synonym at every position, accept improvements, repeat until no change) removes any residual beam-truncation loss. Polish is `O(B · S · W · M)` per pass and typically converges in 2–3 passes.

### 3.4 Full tier-2 loop

```
TIER2_BLOCK_CA(seed, λ, B, b, max_sweeps) -> individual

  assert B >= W                                    # partition validity
  x <- seed
  f <- FULL_EVALUATE(x)                            # includes L0 terms, once
  best <- (x, f)

  for sweep in 1 .. max_sweeps:
      for colour in [RED, BLACK]:
          blocks <- blocks_of(colour)
          parallel for blk in blocks:               # conditionally independent
              w <- LINEARIZE_ASF(f, λ)
              prop <- UPDATE_BLOCK(blk, context(x), w, b)
              prop <- POLISH(prop)
              if ASF_local(prop) >= ASF_local(current(blk)):
                  commit(x, blk, prop)
          # deltas are exact, so the scalarized surrogate value after the
          # sweep is known without re-evaluating; only the L0 terms need a
          # real evaluation.
      f <- FULL_EVALUATE(x)                        # 1 L0 evaluation per sweep
      if ASF(f) > ASF(best.f): best <- (x, f)
      if no block changed: break                   # block-wise local optimum
  return best
```

**Non-decomposable (L0) terms.** Tier 2 does *not* evaluate the LM naturalness or global MFE per move — that would cost a Python round-trip per move and violate the once-per-generation rule from [ADR-0001](../07-decisions/ADR-0001-language-and-runtime-split.md). They are evaluated once per sweep, and the best full-objective incumbent is retained (elitism). The honest statement of the guarantee is therefore:

> Tier 2 is monotone in the decomposable surrogate. Against the full objective it is elitist, not monotone: it returns the best full evaluation it ever saw, and that is at least as good as its input.

---

## 4. Overlap and seam reconciliation — the residual case

Red-black removes reconciliation whenever `B ≥ W`. Two situations break that condition and are the only places where overlapping blocks and reconciliation are still needed:

- **`L < 2W`** — the sequence is shorter than two windows, so no useful partition exists. Handled by optimizing the whole sequence as one block; no reconciliation needed either.
- **A term with an effective interaction radius larger than `B`.** Long-range RNA base pairing is the real instance: two codons 400 nt apart can pair. Enlarging `B` to cover it would defeat the point of blocking.

For that second case the design keeps a fallback: overlapping blocks plus an **exact chain reconciliation**.

Blocks `Bₖ = [k(B−O), k(B−O)+B)` with overlap `O ≥ W`. Each block is optimized independently from the same base, producing proposal `Pₖ`. Consecutive proposals disagree on the overlap region `Sₖ = Bₖ ∩ Bₖ₊₁`.

```
RECONCILE_SEAMS(proposals P₀..P_{K−1}, incumbent x) -> reconciled assignment

  # Choice at each seam: take the left block's proposal, the right block's
  # proposal, or the incumbent. Three options per seam, and the objective
  # contribution of a seam depends only on its own choice and its two
  # neighbours' choices, because interactions span at most O.
  # => a chain DP over seams with 3 states each.

  for k in 0 .. K-1:
      for choice c in {LEFT, RIGHT, INCUMBENT}:
          cost[k][c] <- EXACT_OBJECTIVE_CONTRIBUTION(seam k under c, neighbours)

  best <- CHAIN_DP(cost)                     # O(K · 3² · seam_eval_cost)
  y <- assemble(best)
  y <- POLISH_SEAM_WINDOWS(y, radius = W)    # 1-opt within ±W of each seam
  return y
```

Two properties matter:

- **The reconciled result is never worse than the incumbent**, because `INCUMBENT` is always among the choices at every seam and the DP is exact over those choices.
- **Cost is negligible**: `O(K · 9)` seam evaluations, each `O(O · W · M)`. For `L=666, B=60, O=20`: `K ≈ 13`, so ~117 seam evaluations of ~800 flops each.

`DEP-SRCH-002` — The architecture's reconciliation step is specified here as an exact 3-state chain DP with the incumbent always available, rather than as a heuristic merge. The heuristic version (take whichever proposal scores better locally) can decrease the global objective when two adjacent seams interact, which would break tier 2's convergence argument.

---

## 5. Complexity

Let `L` = editable codons, `S` = mean synonym count (≈2.9 for a typical protein), `W` = window, `B` = block size, `b` = beam width, `M` = number of L1 terms, `N_pop` = tier-2 instance count.

| Quantity | Expression | Default values (`L=666, W=10, B=30, b=64, S=2.9, M=4`) |
| --- | --- | --- |
| One block DP | `O(B · b · S · W · M)` | `30·64·2.9·10·4 ≈ 2.2 × 10⁵` flops |
| One polish pass over a block | `O(B · S · W · M)` | `≈ 3.5 × 10³` flops |
| Blocks per sweep | `L / B` | 22 |
| One full sweep (both colours) | `O(L · b · S · W · M)` | `≈ 4.9 × 10⁶` flops |
| Sweeps to convergence | 5–12 empirically; hard cap 20 | 8 typical |
| **One tier-2 instance** | `O(sweeps · L · b · S · W · M)` | `≈ 3.9 × 10⁷` flops + 8 L0 evaluations |
| **Full tier 2** | `× N_pop` | `200 × 4.0 × 10⁷ ≈ 7.9 × 10⁹` flops + 1 600 L0 evaluations |

`7.9 × 10⁹` flops of simple integer-indexed float work is roughly 3 s single-threaded and ~0.4 s on 8 cores. The **1 600 L0 evaluations dominate**: at ~15 ms each for the LM term that is 24 s serial, or ~1.3 s batched 64-wide on a GPU. Tier 2's wall time is therefore set by the L0 terms, not by the coordinate ascent — which is why `fast` tier disables them entirely (see [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md)).

**Memory.** `O(L · S)` for the cached field table (`666 × 6 × 4 B ≈ 16 KB` per individual) plus the shared model artifacts (~164 KB). Tier 2 is not memory-bound at any realistic construct length.

**Parallelism.** Within a sweep, `L/2B ≈ 11` blocks per colour run concurrently; across instances, all `N_pop` tier-2 runs are independent. Effective parallelism is therefore ≥ `N_pop`, far above any workstation core count, so tier 2 scales linearly with cores up to memory bandwidth.

---

## 6. Determinism

Tier 2 is deterministic given `LandscapeLockfile.rng`. Two requirements, both cheap and both non-negotiable for published results:

1. RNG streams are assigned by **block index and instance index**, never by worker id, so the result does not depend on thread count.
2. The per-block score accumulation uses a **fixed-order binary tree reduction**, because floating-point addition is not associative and a work-stealing reduction order would make identical inputs produce different beams.

`TST-SRCH-004` asserts bit-identical output across 1, 4 and 16 threads.

---

## 7. Tests

| Test | Assertion |
| --- | --- |
| `TST-SRCH-001` | Monotonicity: the scalarized surrogate is non-decreasing across every half-sweep, over 500 randomized runs. |
| `TST-SRCH-002` | Optimality on a brute-forceable instance: on the toy space of `TST-SYN-001` (8 codons, ≤3 synonyms each, 6 561 sequences), tier 1 + tier 2 must find the exact optimum of the decomposable objective in ≥99% of seeds. |
| `TST-SRCH-003` | Seam reconciliation never worsens: for 10⁴ randomized overlapping-block configurations, `RECONCILE_SEAMS` output scores `≥` the incumbent. |
| `TST-SRCH-004` | Thread-count determinism (above). |
