# Algorithm (a) — Fitting the Epistatic Surface by Pseudo-Likelihood

Satisfies: `BRIEF-2.3-a`.
Owning module: `cf-epistasis` (scoring in Rust, fitting in `codonforge.epistasis.train`).
Open question resolved here: [ADR-0003](../07-decisions/ADR-0003-epistasis-window-and-order.md) (window size and interaction order).
Experiments that must run before Phase 3 ships: `EXP-01`, `EXP-02`, `EXP-08`.

---

## 1. What is being modelled, precisely

The architecture proposal says "epistatic fitness surface" without committing to a statistical object. Planning has to commit, so here it is.

**Claim being encoded.** The translational cost of a codon is not a property of the codon; it is a property of the codon *in its local context* — the neighbours occupying the ribosome's E, P and A sites and the surrounding footprint. Codon-pair effects and context-dependent dwell are the empirical basis for this claim; the modelling task is to estimate, per host, a function that scores a codon window rather than a codon.

**Two components, fit separately, exposed as two terms.**

| Component | Data | Statistical object | Exposed term |
| --- | --- | --- | --- |
| **A. Context model** | The host's native CDS set (no Ribo-seq required) | Synonymous-restricted Potts model over codon identities, fit by pseudo-likelihood | `epistasis.context_naturalness` |
| **B. Dwell head** | Ribosome-profiling A-site counts | Negative-binomial GLM over the same window features, plus nuisance covariates | `epistasis.dwell_deviation` |

Splitting them is a departure from the brief's single "epistatic surface" and the reason is data availability: component A needs only an annotated genome, which every host has, while component B needs Ribo-seq, which perhaps a few dozen organisms have. A host without Ribo-seq keeps a real, useful context model instead of losing the layer entirely (`AMD-006`).

`DEP-EPI-001` — The architecture treats the epistatic surface as one model. This plan splits it into a data-cheap context prior and a data-expensive dwell head so that the layer degrades gracefully instead of being all-or-nothing.

---

## 2. Component A: synonymous-restricted Potts model

### 2.1 Parameterization

Let a gene `g` have codons `c₁ … c_L` and amino acids `a₁ … a_L`. Let `Syn(a)` be the synonymous codon set of amino acid `a` (size 1–6). Let `W` be the window half-width in codons.

The conditional distribution of the codon at position `i` given everything else:

```
                exp( F_i(x) )
P(cᵢ = x | c₋ᵢ) = ──────────────────────,      x ∈ Syn(aᵢ)
                  Σ    exp( F_i(x') )
                x'∈Syn(aᵢ)

F_i(x) = h(x) + Σ    [ J_d(x, c_{i+d}) + J_d(c_{i−d}, x) ]
                d=1..W
```

- `h : 64 → ℝ` — codon fields.
- `J_d : 64 × 64 → ℝ` — offset-specific couplings, one matrix per offset `d ∈ {1..W}`. Offset-specific rather than distance-independent because the ribosome's E/P/A geometry is not translation-invariant.
- Parameter count `= 64 + 4096·W`. For `W = 10` that is **41 024** parameters — small enough to fit on a laptop and small enough that regularization is a real constraint rather than a formality.

### 2.2 The one modification that makes this work

Standard DCA/Potts normalizes over the full alphabet. **Here the partition function runs over the synonymous set only.** Three consequences follow, and they are the reason this formulation was chosen over a plain codon-pair-bias table or a plain DCA:

1. The model is *directly* a model of the quantity CodonForge manipulates — synonymous codon choice — rather than a model of sequence that must be reinterpreted.
2. The amino-acid sequence is conditioned on rather than modelled, so the parameters transfer across proteins. This is the generalization bet that `RISK-EPI-001` names and `EXP-01` tests.
3. Amino-acid composition cannot leak into the objective. A model normalized over all 64 codons would assign higher likelihood to sequences rich in common amino acids, which is meaningless when the protein is fixed.

`DEP-EPI-002` — The architecture proposal describes DCA/Potts epistasis modelling generically. This plan narrows it to a *conditional* Potts model normalized over synonymous sets, which is not standard DCA.

### 2.3 Gauge fixing

Potts models are invariant under `J_d(x,y) → J_d(x,y) + u(y)` with a compensating shift in `h`, so the raw parameters are unidentifiable and uninterpretable. Fix the **zero-sum (Ising) gauge within synonymous families**:

```
Σ_{x∈Syn(a)}  J_d(x, y) = 0     ∀ a, y, d
Σ_{y∈Syn(b)}  J_d(x, y) = 0     ∀ b, x, d
Σ_{x∈Syn(a)}  h(x)      = 0     ∀ a
```

Projection onto this gauge is a linear operation applied after each optimizer step, and it is what makes `J_d` entries comparable across offsets — required for the rationale tag `epistasis_favourable_context` to mean anything.

### 2.4 Objective and its one good property

```
ℓ(θ) = Σ_g Σ_i log P(cᵢ^g | c₋ᵢ^g)  −  λ_h‖h‖²  −  λ_J Σ_d ‖J_d‖²  −  μ Σ_d ‖J_d‖_{2,1}
```

The log-conditional is concave in `θ` (log-sum-exp composed with an affine map), the `L2` terms are strictly concave, and the group-`L1` term is concave-but-nonsmooth. Therefore **the fit has a unique global optimum** and there is no seed-dependent model selection to audit. This is worth the parameterization constraints on its own: an epistatic surface whose values depend on a random initialization could not be versioned meaningfully under `ADR-0009`.

- Smooth case (`μ = 0`): L-BFGS.
- Group-sparse case (`μ > 0`, groups = per-offset coupling matrices, used to *select the window*): OWL-QN or proximal L-BFGS. The selected offsets are the empirical answer to "how big should the window be" and feed `ADR-0003`.

### 2.5 Pseudocode

```
FIT_CONTEXT_MODEL(genes G, window W, reg λ_h, λ_J, μ) -> (h, J)

  # ---- feature layout -------------------------------------------------
  # Precompute, once, the neighbour index lists. Positions within W of a
  # gene boundary use a shorter window rather than padding: padding with a
  # sentinel codon injects a fictitious context that the model then learns.
  for g in G:
      for i in 0..L_g-1:
          nbr[g][i] <- { (d, i+d) : 1<=d<=W, i+d < L_g }
                     ∪ { (-d, i-d) : 1<=d<=W, i-d >= 0 }

  θ <- 0                                     # concave objective: 0 is fine
  repeat until converged (L-BFGS / OWL-QN):

      ℓ <- 0 ;  ∇ <- 0

      for g in G:                            # parallel over genes (rayon)
        for i in 0..L_g-1:
            S <- Syn(a[g][i])
            if |S| == 1: continue            # M, W: no information, skip

            # field for each candidate state
            for x in S:
                F[x] <- h[x]
                for (d, j) in nbr[g][i]:
                    F[x] += (d > 0) ? J[ d][x][c[g][j]]
                                    : J[-d][c[g][j]][x]

            p <- softmax(F over S)
            y <- c[g][i]
            ℓ += log p[y]

            # gradient: observed minus expected, the usual softmax form
            for x in S:
                r <- (x == y ? 1 : 0) - p[x]
                ∇h[x] += r
                for (d, j) in nbr[g][i]:
                    if d > 0: ∇J[ d][x][c[g][j]] += r
                    else:     ∇J[-d][c[g][j]][x] += r

      ∇ -= 2λ * θ                            # L2
      θ <- LBFGS_STEP(θ, ℓ, ∇)
      if μ > 0: θ <- PROX_GROUP_L1(θ, μ)     # per-offset group shrinkage
      θ <- PROJECT_ZERO_SUM_GAUGE(θ)

  return (h, J)
```

### 2.6 The exposed term

```
epistasis.context_naturalness(s) = (1/|E|) · Σ_{i∈E} log P(cᵢ | c₋ᵢ)      [nats/codon]
```
where `E` is the set of *editable* positions (locked codons contribute a constant and are excluded so the term is comparable across constructs). Direction: maximize. Its `reference_range` is the p05–p95 of this quantity over the host's native highly-expressed CDS set, which is what makes a raw nats/codon figure interpretable to a user.

---

## 3. Component B: the dwell head

### 3.1 Model

For gene `g`, position `i`, replicate `r`, let `n_{g,i,r}` be the A-site-assigned read count.

```
n_{g,i,r} ~ NegBinomial( mean = μ_{g,i,r},  dispersion = φ )

log μ_{g,i,r} = α_{g,r}                      # gene x replicate offset
              + Σ_{d=1..W} ψ_d(cᵢ, c_{i+d})  # context couplings (the epistasis)
              + ω(cᵢ)                        # codon main effect
              + f_ramp(i / L_g)              # position spline, 5' ramp
              + β₁ · dG_local(i)             # folding energy of the window ahead
              + β₂ · charge_tunnel(i)        # nascent-peptide charge, codons i-30..i-10
              + β₃ · n_pro(i)                # proline in P/A sites
              + β₄ · SD_affinity(i)          # internal Shine-Dalgarno (bacteria only)
              + batch_r
```

**`α_{g,r}` is not optional.** Without a gene-level offset, transcript abundance and initiation rate — which have nothing to do with elongation context — are absorbed into codon coefficients, and the model learns that "codons found in highly expressed genes are fast". That single omission would invalidate the entire layer.

The nuisance covariates `β₁…β₄` exist for the same reason: they are the known non-context drivers of local dwell. Leaving them out does not make the model simpler, it makes the couplings wrong.

### 3.2 Coupling to component A

`ψ_d` and `J_d` describe the same physical window. Three options were considered; the plan picks the third and makes the choice ablatable:

| Option | Verdict |
| --- | --- |
| Fit `ψ` freely | Wastes the (much larger) native-CDS data that constrains `J`. |
| Force `ψ_d = γ_d · J_d` | Too strong: what is *frequent* in native genes is not identical to what is *fast*, and conflating them re-introduces the codon-usage fallacy CodonForge exists to move past. |
| **Empirical-Bayes shrinkage: `ψ_d ~ N(γ_d·J_d, τ_d²)`, with `γ_d, τ_d` fit** | Chosen. Data-poor hosts shrink toward the context model; data-rich hosts pull away from it. `γ_d` is itself an interesting diagnostic — how much "frequent" and "fast" actually agree in this host. |

`EXP-02` measures whether the shrinkage helps on held-out genes versus free fitting; if it does not, the prior is dropped and `ψ` is fit freely.

### 3.3 Fitting

Penalized maximum likelihood via IRLS-within-L-BFGS on the NB likelihood, with `φ` estimated by profile likelihood on a held-out gene split. Splits are **by gene, never by codon** — adjacent codons within a gene are strongly dependent, so a random-codon split leaks and produces an inflated held-out correlation. `ModelManifest.training.split_strategy` records this explicitly so the mistake is visible in review.

### 3.4 The exposed term, and what "better" means

This is where the plan must make a decision the architecture leaves open. `dwell_deviation` is **not** "minimize total predicted dwell". Uniformly fast translation is exactly the design error the architecture argues against, so an objective that rewards it would undo the layer's purpose.

```
epistasis.dwell_deviation(s) =
      w₁ · RMS_i ( z_i(s) − z_ref(i / L) )          # rhythm mismatch
    + w₂ · Σ_i  softplus( z_i(s) − z_pause )        # unintended strong pauses
    − w₃ · (pause mass inside user- or predictor-designated pause windows)
```

- `z_i(s)` is the predicted dwell for the candidate, standardized against the host's native per-codon dwell distribution.
- `z_ref(·)` is the host's mean native profile as a function of *relative* position — the canonical 5′ ramp and its decay.
- `z_pause` is the host's 99th native percentile; the `softplus` makes isolated extreme pauses expensive without penalizing normal variation.
- The third term is what allows `cf-biophys-kinetics` to *ask* for a pause at a domain boundary rather than fight for one.

Direction: minimize. Unit: `z_score`.

`DEP-EPI-003` — The architecture does not state what the epistatic surface should optimize toward. This plan defines it as concordance with the host's native translational rhythm plus suppression of unintended strong pauses, not as dwell minimization.

---

## 4. Delta computation (the L1 capability)

This is the part that makes the whole search affordable, so it is specified rather than left to the implementer.

Maintain, for every position `j` and every candidate state `y ∈ Syn(a_j)`, the cached field `F_j(y)` from §2.1 and the cached normalizer `Z_j`.

Substituting `c_i → x'` changes `F_j(·)` only for `|i − j| ≤ W`, and for each such `j` the change is a single coupling term:

```
DELTA_SINGLE(i, x') -> Δ(context_naturalness), Δ(dwell terms)

  Δ <- 0
  # (1) the substituted position itself
  Δ += [F_i(x') − F_i(cᵢ)] − [log Z_i(x') − log Z_i(cᵢ)]     # Z_i unchanged: same Syn set

  # (2) each neighbour's conditional, whose context just changed
  for d in 1..W:
      for j in {i−d, i+d} ∩ [0, L):
          for y in Syn(a_j):
              F'_j(y) = F_j(y) + ΔJ(y, cᵢ → x', offset)
          Δ += [F'_j(c_j) − log Σ_y exp F'_j(y)]
             − [F_j(c_j)  − log Σ_y exp F_j(y)]
  return Δ
```

**Cost `O(W · |Syn|)` ≈ 60 flops.** Compare with `O(L · W)` ≈ 6 600 flops for a full rescore of a 666-codon CDS: a **~110×** saving per move, which is the concrete performance argument in [ADR-0002](../07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md).

Full delta table over all positions and alternatives: `O(L · W · |Syn|²)` ≈ `666 × 10 × 36 ≈ 2.4 × 10⁵` flops — sub-millisecond, and it is computed once per coordinate-ascent sweep rather than once per move.

The same caching scheme applies to the dwell head, with the addition that `dG_local` and `charge_tunnel` covariates must be recomputed over their own windows; `dG_local` is the expensive one and is served by the incremental structure cache described in [`08-complexity-and-performance-budgets.md`](08-complexity-and-performance-budgets.md).

---

## 5. Higher-order terms

Full third-order couplings would need `64³ × C(W,2) ≈ 1.2 × 10⁷` parameters against `~10⁶` codon observations. That is not a modelling choice, it is a category error. The plan therefore uses **gated, selected triplets**:

```
SELECT_TRIPLETS(pairwise model, data, K, δ) -> accepted triplets

  1. Fit the pairwise model; compute Pearson residuals r_{g,i}.
  2. For each candidate triplet pattern t = (x, y, z, d₁, d₂) with d₁<d₂<=W,
     accumulate the score statistic  U(t) = Σ r_{g,i} · 1[pattern t occurs at i]
     and its variance under the pairwise null.
  3. Rank by U²/Var(U); keep the top K (default K = 2000).
  4. Refit with those K triplet parameters included.
  5. ACCEPT iff held-out log-likelihood per codon improves by >= δ
     AND the bootstrap 95% CI of the improvement excludes 0.
     Otherwise discard all triplets and ship the pairwise model.
```

Step 5 is a hard gate, not a preference: a triplet extension that cannot beat its own held-out baseline is an overfit, and shipping it would silently degrade every design produced for that host. `EXP-08` runs this procedure for E. coli and yeast and reports whether step 5 ever passes.

---

## 6. Complexity and resource envelope

| Operation | Time | Memory | Notes |
| --- | --- | --- | --- |
| Component A fit, one gradient evaluation | `O(N · W · S)` | `O(64 + 4096·W)` params | `N` = total codons, `S` = mean synonym count (≈2.9). E. coli: `1.3×10⁶ × 10 × 2.9 ≈ 3.8×10⁷` flops. |
| Component A full fit | `≈ 200` iterations × above ≈ `7.6×10⁹` flops | < 200 MB | Minutes on 32 vCPU, vectorized. |
| Component B fit (NB GLM) | `O(I · N · P_eff)` | `O(N)` rows streamed | `P_eff` ≈ active features per row ≈ `2W + 8`. |
| Triplet selection | `O(N · W² · S)` scan + refit | `O(K)` | The scan is the expensive part: `1.3×10⁶ × 100 × 2.9 ≈ 3.8×10⁸`. |
| Score full sequence | `O(L · W · S)` | `O(L · S)` cache | 666 codons: ≈ 2×10⁴ flops. |
| Score single substitution | `O(W · S)` | — | ≈ 60 flops. |
| Full delta table | `O(L · W · S²)` | `O(L · S)` | ≈ 2.4×10⁵ flops. |

Artifact size: fields + couplings at `float32` = `(64 + 4096·10) × 4 B ≈ 164 KB`. Small enough that a host's epistatic surface is cheap to version, cheap to ship and cheap to keep several versions of — which is what `ADR-0009` assumes.

---

## 7. Validation (see [`../05-testing/01-test-plan.md`](../05-testing/01-test-plan.md))

| Test | What it establishes |
| --- | --- |
| `TST-EPI-001` | Parameter recovery: generate sequences from a known synthetic Potts surface, refit, require Pearson `r ≥ 0.95` on `J` entries in zero-sum gauge. Catches gauge, indexing and gradient bugs without any biology. |
| `TST-EPI-002` | Delta correctness: `DELTA_SINGLE` must equal `full_score(after) − full_score(before)` to within `1e-9` for 10⁴ random substitutions. This is the single most important test in the package — an incorrect delta silently corrupts every tier-2 result while leaving tier-3 results looking plausible. |
| `TST-EPI-003` | Held-out dwell prediction beats a codon-usage-only baseline by `≥ 0.05` Spearman with a bootstrap CI excluding zero (`GATE-EPI-001`). |
| `EXP-01` | Cross-dataset and cross-protein generalization — the falsification test for `RISK-EPI-001`. |
