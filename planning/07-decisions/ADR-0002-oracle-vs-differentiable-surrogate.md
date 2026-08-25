# ADR-0002 — Landscape Model as Oracle vs. Differentiable Surrogate

- **Status:** Accepted
- **Resolves:** open question 3.2 (`BRIEF-3.2`)
- **Related:** [ADR-0001](ADR-0001-language-and-runtime-split.md), [`../03-algorithms/02-block-local-coordinate-ascent.md`](../03-algorithms/02-block-local-coordinate-ascent.md), [`../02-schemas/landscape_query.schema.json`](../02-schemas/landscape_query.schema.json)

## Question as posed

> Should the search engine treat the landscape model as a black-box scoring oracle (favouring modularity, per the architecture's swappable-layer goal) or should parts of it (e.g. the structure predictor) expose gradients to allow faster local search? What is the performance cost of each choice at construct-level sequence lengths?

## Decision

**The contract is a black-box oracle. Gradients are an optional capability that nothing depends on. The thing that actually makes local search fast is not gradients — it is exact deltas.**

Three capability levels, declared per term in `ModelManifest.capability` and echoed in every `LandscapeResponse`:

| Level | Operation | Who implements it |
| --- | --- | --- |
| **L0 Scalar** | score a full sequence | mandatory for every predictor |
| **L1 Delta** | exact `Δterm` for one synonymous substitution, and a full delta table | every decomposable term: CAI, GC, k-mer, Potts context, dwell head, tRNA/burden ledger, windowed structure accessibility |
| **L2 Gradient** | `∂term/∂(relaxed codon logits)` | optional; currently only the LM naturalness term could offer it |

The search engine **requires** L0 and **exploits** L1. It never requires L2.

## Why gradients are the wrong tool here

The design space is discrete. A codon is a categorical choice among 1–6 synonyms; there is no continuous path between `CGT` and `AGA`. A gradient therefore only exists with respect to a continuous relaxation, and once you have committed to a relaxation you are optimizing a different problem and must round back at the end.

For a decomposable term, the **exact** finite change is available for the same asymptotic cost as a directional derivative — and it is exact rather than a linearization. There is no regime in which the gradient is better information.

| | Gradient (L2) | Exact delta (L1) |
| --- | --- | --- |
| What it says | rate of change under an infinitesimal, fictitious relaxation | the true change from the actual move |
| Accuracy for the move you will make | approximate; error grows with the categorical gap | exact |
| Cost for the Potts term | one backward pass over the window: `O(W·S)` | `O(W·S)` |
| Requires an autograd runtime | yes | no |
| Requires all predictors in one framework | yes, if used end-to-end | no |

The last row is decisive against a differentiable landscape. Making the whole stack differentiable means every predictor lives in one autograd runtime — the structure kernel, the Potts model, the tRNA ledger and the LM would all have to be, say, PyTorch modules. That directly contradicts the swappable-layer goal: a contributor could no longer drop in a C-library folding engine or a hand-written Rust ledger.

## Performance cost, measured against the plan's own budgets

At `L = 666` codons (2 kb), `W = 10`, `S = 2.9`, `M₁ = 4` terms:

| Operation | Cost | Ratio |
| --- | --- | --- |
| Full rescore (L0 only, decomposable terms) | `O(L·W·S·M₁)` ≈ `7.7 × 10⁴` flops | 1× |
| Exact delta for one substitution (L1) | `O(W·S·M₁)` ≈ `120` flops | **~640× cheaper** |
| Full delta table over all positions (L1) | `O(L·W·S²·M₁)` ≈ `2.2 × 10⁵` flops | one table ≈ 3 rescores, and it answers ~1 900 move queries |

Tier 2 issues roughly `1.5 × 10⁴` delta queries per instance per restart. With L1 that is `1.8 × 10⁶` flops; with L0-only rescoring it is `1.2 × 10⁹` — a **660×** difference, which turns a sub-second tier-2 phase into a ten-minute one.

At construct scale (`L = 5 000`, 15 kb) the ratio grows linearly with `L`, to roughly `4 800×`. This is the whole argument: **the oracle/delta distinction matters far more than the oracle/gradient distinction, and it scales in the direction that hurts.**

For an L0-only term (the LM), tier 2 simply does not evaluate it per move; it is evaluated once per sweep and once per generation. The cost of a term being L0-only is therefore not a slowdown, it is a *reduction in the fidelity of tier 2's search direction* — which is acceptable because tier 3 optimizes the full objective and tier 2 is only a seeder.

## What the design gives up

- **Gradient-based local search over a relaxation** (e.g. optimizing continuous codon logits then rounding) is not available as a first-class strategy. It could be added later behind L2 without changing any contract, and `EXP-09` is registered to test whether it beats block coordinate ascent on the benchmark set. Current expectation: it will not, because the rounding step discards exactly the constraint-feasibility structure the discrete search maintains.
- **End-to-end differentiable design** (backpropagating from an objective to a sequence, as some mRNA-design work does) is out of scope for v1.

## Consequences

- `LandscapeQuery.target` has three shapes (`full`, `delta`, `delta_table`) rather than one. This is the main complexity cost of the decision, and it is confined to one schema.
- Every L1 predictor must pass a delta-correctness test asserting `delta == score(after) − score(before)` to `1e-9`. This is the highest-value test in the package: a wrong delta corrupts every tier-2 result while leaving tier-3 results looking plausible, so it fails silently without this test.
- Predictors may be written in any language and shipped as any artifact format. The swappable-layer goal survives intact.
- `LandscapeResponse.gradients` remains in the schema, fully specified, unused by the engine. Specifying it now costs nothing and prevents a future contributor from inventing an incompatible gradient convention.

## Departures

`DEP-LAND-001` — The architecture proposal frames the choice as oracle vs. differentiable surrogate. This plan introduces a third option — an exact-delta capability — and shows it dominates the gradient option for every decomposable term, which is most of them.
