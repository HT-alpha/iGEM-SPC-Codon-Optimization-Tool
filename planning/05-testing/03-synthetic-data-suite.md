# Synthetic-Data Test Suite

Satisfies: `BRIEF-2.5-b` — "a synthetic-data test suite (small toy proteins/hosts with known ground truth) so core algorithms can be validated before any real biological training data is available".
Owning module: `cf-testkit`.

---

## 1. What synthetic data is for here

Real biological validation cannot start until Ribo-seq is ingested, models are fit and experiments run — which is Phase 3 at the earliest. But the algorithms that need validating most urgently are built in Phase 1: the constraint machinery, the coordinate ascent, the delta computations. Waiting for biology to test them would mean building four phases on unverified foundations.

Synthetic data gives **ground truth without biology**. Each generator below is constructed so that the correct answer is known exactly, by construction or by exhaustive enumeration, so a failure is unambiguous.

The suite is also what makes the four "no reference implementation" layers reviewable: a claim like "the joint search beats independent optimization" is a hypothesis on real data and a *theorem* on `SYN-03`, where the margin is computable in closed form.

---

## 2. The generators

### `SYN-01` — Toy exhaustive design space

**Construction.** An 8-residue protein whose amino acids each have ≤3 synonymous codons (e.g. `M-K-R-I-L-D-E-*` restricted to a reduced codon table). The design space is `3⁸ = 6 561` sequences — small enough to enumerate completely.

**Ground truth.** The exact Pareto front, computed by brute force over all 6 561 sequences under any objective set.

**Validates.** `TST-SYN-001`, `TST-SRCH-002`. This is the only place in the system where a search result can be compared against the *true* optimum rather than against another algorithm. Every search change is regression-tested here first: IGD+ to the true front must be 0 for tier 1 + tier 2 on the decomposable objective, and the tier-3 front must contain every true extreme point.

**Why it is small on purpose.** A space of 6 561 is enumerable in milliseconds, so this test can run on every commit.

### `SYN-02` — Synthetic Potts host

**Construction.** Draw `h` and `J_d` (`d = 1..W`, `W = 3`) from a known distribution with the zero-sum gauge imposed. Generate `N_genes` synthetic CDSs by Gibbs sampling from the resulting synonymous-conditional distribution, given random amino-acid backbones.

**Ground truth.** The generating parameters.

**Validates.** `TST-SYN-002`, `TST-EPI-001` (parameter recovery: Pearson `r ≥ 0.95` on `J` in zero-sum gauge). Also, by sweeping `N_genes`, it gives the **sample-size curve** that sets the data tiering thresholds in [ADR-0003](../07-decisions/ADR-0003-epistasis-window-and-order.md) §2 — "how much data do you need before the couplings are recoverable" is answerable synthetically, before any real Ribo-seq is available.

Extended to generate synthetic A-site counts by drawing `n_i ~ NegBinomial(exp(α_g + Σψ_d(...)), φ)`, which gives the dwell head the same treatment.

### `SYN-03` — Two-cistron contention instance

**Construction.** Two cistrons of equal length whose amino-acid composition forces both to draw on a single scarce isoacceptor. The synonymous families are constructed so that each cistron alone would take the scarce pool, but the joint optimum splits demand.

**Ground truth.** The optimal split, and the exact margin by which the joint optimum beats the best independent-optimization result, computable in closed form because the ledger is a linear function of codon counts and the penalty is a known scalar function.

**Validates.** `TST-SYN-003`, `TST-SRCH-009`. **This is the test that justifies the systems-level layer's existence.** If joint search cannot beat independent optimization on an instance constructed so that it must, the layer is not doing what it claims. A parameterized version (varying scarcity and cistron size) also produces the oscillation instance for `TST-SRCH-008`.

### `SYN-04` — Planted constraint violations

**Construction.** Take a feasible sequence; inject a violation of a chosen class (restriction site, homopolymer run, k-mer repeat, GC-window excursion) at a chosen position, in a context engineered so that the minimal repair distance is known — either a unique single-codon fix, or a deliberately-constructed set of `k` equally-minimal fixes.

**Ground truth.** The minimal edit distance and the full set of minimal fixes.

**Validates.** `TST-SYN-004`, `TST-CONS-002` (minimality **and** the diversity property: across 1 000 seeds, more than one of the `k` minimal fixes must be returned). The multi-fix construction is what makes the randomized tie-breaking testable rather than merely asserted.

### `SYN-05` — Infeasible constraint sets

**Construction.** One purpose-built input per feasibility check: a blocking violation inside immutable backbone; a k-mer limit violating the pigeonhole bound; a position with no admissible synonym; a motif set that closes off every DP continuation at a known position; a GC band outside the attainable range.

**Ground truth.** Which check should fire, at which position, naming which constraint.

**Validates.** `TST-CONS-006`. Also serves as the fixture set for the `/v1/feasibility` API tests, so the user-facing message for each failure is exercised, not just the detection.

### `SYN-06` — Synthetic host profiles across the data tiers

**Construction.** Three variants of one synthetic host: tier 1 (with a fabricated burden/expression panel), tier 2 (ribosome counts and elongation rate, no panel), tier 3 (genome-derived codon usage only). Identical underlying "true" biology; different visible data.

**Ground truth.** The true capacity, known because the host was generated from it.

**Validates.** `TST-BURD-005`, and — more valuably — it makes the *degradation* measurable: how far does the tier-3 estimate drift from the truth the tier-1 fit recovers? That number is the honest expectation to publish for tier-3 hosts, and it is obtainable synthetically well before `EXP-04` can be run on real data.

### `SYN-07` — Property-based strategies

Hypothesis (Python) and proptest (Rust) strategies rather than fixed fixtures:

| Strategy | Generates |
| --- | --- |
| `arbitrary_protein` | Random proteins with configurable length and composition |
| `arbitrary_cds` | A CDS in frame, with a valid stop, optionally with flanks |
| `arbitrary_substitution` | A synonymous substitution at an unlocked position |
| `arbitrary_constraint_set` | A set of constraints, optionally guaranteed satisfiable |
| `arbitrary_ledger` | A tRNA supply/demand state within physiological bounds |
| `arbitrary_front` | An `M`-objective candidate set with a known non-dominated subset |

These drive the invariant tests — protein invariance, delta correctness, ledger conservation, repair termination — which are the tests most likely to catch a subtle regression during refactoring.

---

## 3. What synthetic data cannot do

Stated so the suite is not over-trusted:

- It cannot validate that the model class is right for biology. `SYN-02` proves the fitter recovers a Potts surface; it says nothing about whether real ribosome density is Potts-like. Only `EXP-01` can.
- It cannot validate the objective *definitions*. `SYN-03` proves the joint search optimizes the ledger; it cannot say whether the ledger reflects real contention.
- It cannot calibrate uncertainty against reality. Synthetic coverage is coverage of the generating process, which is optimistic by construction.

The suite's job is to guarantee that when a real-data experiment fails, the failure is in the science and not in the code. That is a narrower claim than "the system is validated", and it is the claim the suite can actually support.

---

## 4. Cost

Every generator is `O(milliseconds)` and every fixture fits in a few kilobytes, so the entire synthetic suite runs on every commit. This is deliberate: a ground-truth suite that only runs nightly stops being the thing developers check against.
