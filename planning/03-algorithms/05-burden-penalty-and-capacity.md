# Algorithm (e) — The Saturating Burden Penalty and Host-Capacity Estimation

Satisfies: `BRIEF-2.3-e`.
Owning module: `cf-burden` (penalty in Rust; capacity calibration in `codonforge.burden.calibrate`).
Related: [ADR-0004](../07-decisions/ADR-0004-resource-ledger-granularity.md), [`04-systems-level-joint-search.md`](04-systems-level-joint-search.md) §2 for the ledger this consumes.

---

## 1. Two penalty functions, because one cannot do both jobs

The brief asks for "the saturating burden-penalty function". Saturation is biologically correct — growth-rate cost cannot exceed 100 %, and a cell 10× over its translational capacity is not 5× sicker than one 2× over. But a saturating function is a bad search signal: once utilization is well past the knee, the function is flat, its delta is ~0, and the optimizer cannot tell a bad design from a catastrophic one. Coordinate ascent in that regime does nothing.

The design therefore keeps **two** functions with different jobs, and says so explicitly rather than letting an implementer discover the flat-gradient problem during Phase 2.

### 1.1 The reported objective — bounded and saturating

```
Per-resource utilization:      u_α  = D_α / S_α                       (α ∈ pools ∪ {ribosome})
Soft-max aggregation:          u_eff = ( Σ_α ( u_α / θ_α )^p )^(1/p),   p = 4
Growth-cost penalty:           g(u_eff) = g_max · u_eff^n / (1 + u_eff^n),   n = 2 (default)
```

- **`p = 4` soft-max aggregation** rather than a plain mean or a hard `max`. Burden is a bottleneck phenomenon: one exhausted isoacceptor throttles everything, so a mean under-reacts. A hard `max` over-reacts and is non-differentiable, which destroys the delta table. `p = 4` sits between them and is smooth. `p` is a host parameter, fit at capacity tier 1. The formula line above is an L_p norm, which is not between a mean and a max. The implementation uses the power mean `(mean of (u/θ)^p)^(1/p)`. `DEP-BURD-005`.
- **`g_max ∈ [0, 1]`** is the asymptotic relative growth-rate cost. Bounding it is what makes the penalty saturating.
- Reported as `burden.growth_cost`, unit `relative_growth_rate`, direction `minimize`.

### 1.2 The internal search signal — unbounded soft hinge

```
h(u_α) = (1/κ) · log( 1 + exp( κ · (u_α/θ_α − 1) ) ),      κ = 8 (default)
H      = Σ_α  w_α · h(u_α)
```

`h` is ≈0 below capacity, ≈linear above it, and smooth everywhere. `H` is used:

- as the price basis for the shadow-price computation in the joint search;
- as an `advisory` constraint-violation contribution when a design is grossly over capacity;
- as a tie-break inside tier-2 coordinate ascent when `g` has saturated.

`H` is never reported to the user as an objective value: it has no interpretable unit. It appears in `LandscapeResponse.per_position` diagnostics only.

`DEP-BURD-002` — The architecture proposal specifies one saturating burden penalty. This plan ships a bounded saturating function for reporting and an unbounded soft hinge for search, because a saturating objective alone provides no gradient in the regime where it matters most.

---

## 2. Estimating the capacity threshold `θ`

Three tiers, selected automatically from `HostProfile.data_completeness`. The tier is recorded in `HostProfile.capacity.estimator_tier`, propagated into `LandscapeResponse.terms[...].provenance`, and surfaced to the user — because a tier-3 burden number means something categorically different from a tier-1 one.

### 2.1 Tier 1 — fit to measured burden/expression data

Requires a panel of constructs measured in the same host with both an expression readout and a growth-rate readout (≥20 constructs spanning ≥1 order of magnitude in expression).

```
FIT_CAPACITY_TIER1(observations O, host) -> (θ, n, g_max, p, uncertainty)

  # each observation gives (predicted demand vector D_i, measured relative growth γ_i)
  for obs in O:
      D_i    <- LEDGER_DEMAND(obs.sequence, obs.f_estimate, host.w)
      γ_i    <- obs.growth_rate / obs.control_growth_rate

  minimize over (θ, n, g_max, p):
      Σ_i  ( (1 − γ_i)  −  g(u_eff(D_i; θ, p); n, g_max) )² / σ_i²
      subject to  0 <= g_max <= 1,  n in [0.5, 6],  p in [1, 16]

  uncertainty <- PARAMETRIC_BOOTSTRAP(2000 resamples)
  return parameters with per-parameter CI
```

This is the only tier that produces a genuine *growth-rate prediction*. Everything else is a calibrated proxy.

### 2.2 Tier 2 — derived from absolute physiology

Requires `ribosomes_per_cell`, `elongation_rate_codons_per_s` and a doubling time — available for the handful of model organisms with published growth-law measurements.

```
FIT_CAPACITY_TIER2(host) -> (θ, uncertainty)

  C_total   <- ribosomes_per_cell · active_fraction · elongation_rate     # aa/s per cell
  D_native  <- (protein_mass_per_cell / mean_aa_mass) · ln2 / doubling_time_s
  θ_global  <- C_total / D_native                                        # ~1.1-1.3 typically

  # per-pool split
  if has_trna_abundance:  S_α <- abundance_α · charging_turnover_α
  else:                   S_α <- (gcn_α / Σ gcn) · C_total

  # the dominant uncertainties are ribosomes_per_cell and active_fraction;
  # propagate them rather than reporting a point estimate
  uncertainty <- DELTA_METHOD(σ_ribosomes, σ_active_fraction, σ_elongation)
  return θ_global, S, uncertainty
```

`g_max` is not identifiable at this tier and is fixed at the literature default `0.6` with a wide prior; this is recorded in the manifest so it is not mistaken for a fit value.

### 2.3 Tier 3 — relative-to-native normalization (the no-data case)

This is the tier that most hosts land in, and getting it right is the difference between CodonForge being usable outside model organisms and not.

**Assumption:** a host's tRNA pools are, on average, matched to its own native codon demand. This is the same empirical regularity that underlies tAI — tRNA gene copy number tracks genomic codon usage — so tier 3 is not a fabrication, it is a documented approximation.

```
FIT_CAPACITY_TIER3(host) -> (S, θ = 1, semantics = "profile_mismatch")

  # expression-weighted native demand profile from the genome alone
  for each annotated CDS gene:
      weight_g <- CAI(gene, highly_expressed_reference)^ρ      # ρ = 2, expression proxy
  D_native_α <- Σ_g weight_g · Σ_x n_g(x) · w(x, α)

  S_α  <- D_native_α / Σ_β D_native_β          # supply assumed proportional to native demand
  θ_α  <- 1                                     # by construction
  return S, θ
```

At tier 3 the quantity `u_α = D_α/S_α` is **how much the construct's demand profile deviates from the host's native profile** — a mismatch score, not a utilization ratio. The design makes three consequences explicit:

1. `HostProfile.capacity.theta_unit` is `fraction_of_native_demand`, not an absolute unit.
2. `g_max` is set to 1 and `burden.growth_cost` is documented as a dimensionless mismatch score in `[0,1)`, **not** a predicted growth-rate cost. The API returns unit `dimensionless` at this tier.
3. `Uncertainty.calibrated = false` and the decision layer attaches the `low_confidence_host` rationale tag.

`DEP-BURD-003` — The architecture proposal treats the burden model as predicting host burden. This plan states that for hosts without proteomics or tRNA abundance data (tier 3) the burden term is a demand-profile mismatch score and must not be presented as a growth-rate prediction. Reporting the same unit at all three tiers would be the single most misleading thing CodonForge could do.

### 2.4 User configuration

`θ`, `κ`, `n`, `g_max` and `p` are all overridable per job. A user with their own burden measurements for their own strain has better information than any tier here, and the design should not make them fight the defaults. Overrides are recorded in the lockfile.

---

## 3. Incremental evaluation

```
DELTA_BURDEN(position i, cᵢ -> x', cistron k) -> Δg, Δh

  for α in decode(cᵢ) ∪ decode(x'):                    # <= 3 pools
      D_α += f_k · ( w(x', α) − w(cᵢ, α) )
      u_α  = D_α / S_α
  recompute u_eff over the <=3 touched pools using cached Σ (u_α/θ_α)^p
  Δg = g(u_eff_new) − g(u_eff_old)
  Δh = Σ_touched w_α · ( h(u_α_new) − h(u_α_old) )

  # ribosome channel: elongation rate at position i changed, and so did the
  # rate at EVERY position decoded by a touched pool. That second part is
  # O(#positions using a touched pool), so it is deferred: the ribosome
  # occupancy R is recomputed exactly once per sweep, not per move, and the
  # per-move delta uses the frozen availability vector.
  return Δg, Δh
```

The last comment is a real, load-bearing approximation and is called out rather than hidden: within one tier-2 sweep the availability vector `a` is held fixed, making per-move burden deltas exact with respect to a frozen `a` and approximate with respect to the true self-consistent `a`. `a` is refreshed at every sweep boundary and exactly at every tier-3 evaluation. `TST-BURD-003` bounds the resulting error at `< 1 %` of the term's reference range for realistic constructs, and fails the build if a change makes it worse.

---

## 4. Complexity

| Operation | Cost |
| --- | --- |
| Full penalty from a ledger | `O(|A|)` ≈ 40 flops |
| Per-move delta (frozen `a`) | `O(1)` — at most 3 pools |
| Availability refresh (fixed point) | `O(iters · |A|)` ≤ 2 000 flops |
| Ribosome occupancy recompute | `O(L_tot)` |
| Tier-1 capacity fit | `O(iters · |O| · L̄)`, `|O|` ≈ 20–200 observations — seconds |
| Tier-3 native profile | `O(genome codons)` ≈ `1.4 × 10⁶` — once per host, cached in the artifact |

The burden layer is computationally free relative to everything else in the system. Its cost is entirely in *data*, which is why §2 is three times longer than §1.

---

## 5. Tests and gates

| Id | Assertion |
| --- | --- |
| `TST-BURD-003` | Frozen-availability delta error bound (§3). |
| `TST-BURD-004` | Monotonicity and boundedness: `g` is non-decreasing in each `u_α`, and `0 ≤ g ≤ g_max` over 10⁶ random ledgers. |
| `TST-BURD-005` | Tier-3 semantics: a tier-3 host's response carries unit `dimensionless`, `calibrated = false`, and the `low_confidence_host` tag. Asserting this in a test is the only way to keep §2.3's honesty requirement from eroding. |
| `GATE-BURD-001` | Back-test against a published burden/expression dataset: Spearman `≥ 0.5` between predicted `g` and measured relative growth cost, and a positive margin over a "total protein output only" baseline with a bootstrap CI excluding zero. |
| `EXP-04` | Does tier-3 profile mismatch correlate with measured burden at all? If not, the tier-3 term is demoted from a Pareto axis to a diagnostic for tier-3 hosts. |
