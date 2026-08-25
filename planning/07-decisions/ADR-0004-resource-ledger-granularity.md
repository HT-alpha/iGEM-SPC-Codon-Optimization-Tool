# ADR-0004 — Resource-Ledger Granularity and Host-Capacity Estimation

- **Status:** Accepted; time resolution deferred to `EXP-03`
- **Resolves:** open question 3.3 (`BRIEF-3.3`)
- **Related:** [`../03-algorithms/04-systems-level-joint-search.md`](../03-algorithms/04-systems-level-joint-search.md), [`../03-algorithms/05-burden-penalty-and-capacity.md`](../03-algorithms/05-burden-penalty-and-capacity.md)

## Question as posed

> At what resolution should tRNA/ribosome demand be tracked — per-isoacceptor pool only, or also per-timepoint/expression-phase? How is host-capacity estimated for organisms lacking published tRNA abundance or proteomics data?

## Decision

### Granularity

| Dimension | v1 | Rationale |
| --- | --- | --- |
| tRNA | **Per-isoacceptor pool (anticodon level)**, with wobble decoding shares | This is the level at which synonymous choice acts. Finer (per tRNA gene) adds parameters with no additional decision-relevant information, since isoacceptor genes are not separately regulatable by codon choice. Coarser (per amino acid) makes the ledger *constant* — see below. |
| Amino-acid pools | **Not tracked** | Protein sequences are fixed, so total amino-acid demand is invariant under synonymous redesign. An amino-acid ledger would be a constant column. |
| Energy (ATP/GTP) | **Not tracked** | Proportional to total amino-acid demand, hence also invariant. |
| Ribosomes | **Tracked, as one aggregate occupancy number** `R` | This one *is* reducible in absolute terms: slower elongation means more ribosome-seconds per protein. It is the only channel through which codon choice reduces total resource consumption rather than redistributing it. |
| Time | **Single steady state** in v1; a two-phase (pre/post-induction) ledger is an opt-in extension gated on `EXP-03` | See below. |
| Charged vs. uncharged tRNA | Not tracked separately; `charged_fraction` scales supply when known | Charging dynamics need kinetic parameters essentially no host has. |
| Chaperones, secretion, membrane insertion | Not tracked | Real burden channels, but not ones synonymous codon choice controls directly. Registered as `RISK-BURD-003` for v2. |

The invariance argument deserves to be stated once, plainly, because it bounds what the layer can honestly claim: **synonymous redesign cannot make a cell synthesize fewer amino acids.** It can only (a) move demand between isoacceptors within a synonymous family and (b) change how long a ribosome is occupied. Anything the burden model appears to save beyond those two channels is a bug.

### Why not per-timepoint

Per-timepoint tracking requires an induction model: how fast the promoter turns on, how transcript levels rise, how the tRNA pools re-equilibrate. That needs kinetic parameters that even *E. coli* only partially has, and no non-model host has at all. Introducing it in v1 would mean either (i) fabricating those parameters, or (ii) restricting the whole systems-level layer to one or two organisms.

Instead, `EXP-03` answers the empirical question that actually matters: **does timepoint resolution change the design?** It runs the two-phase ledger against the steady-state ledger on the benchmark construct set and measures the divergence between the resulting Pareto fronts. If designs agree, the extension is not built. If they diverge materially, the extension is built and gated to hosts that have the parameters.

This is the pattern the brief asks for: rather than picking arbitrarily, name the experiment that resolves it.

## Host-capacity estimation for data-poor organisms

Three tiers, selected automatically from `HostProfile.data_completeness`, fully specified in [`../03-algorithms/05-burden-penalty-and-capacity.md`](../03-algorithms/05-burden-penalty-and-capacity.md) §2:

| Tier | Requires | Produces | Reported as |
| --- | --- | --- | --- |
| 1 | ≥20 constructs with paired expression + growth measurements | Fitted `θ, n, g_max, p` with bootstrap CIs | A genuine predicted relative growth cost |
| 2 | Ribosome count, elongation rate, doubling time | `θ` from `C_total / D_native`; `g_max` fixed at a literature default | A calibrated proxy for growth cost |
| 3 | **An annotated genome. Nothing else.** | Supply assumed proportional to the host's own expression-weighted native codon demand; `θ ≡ 1` by construction | A **demand-profile mismatch score**, explicitly *not* a growth-rate prediction |

Tier 3 is the answer to "organisms lacking published tRNA abundance or proteomics data", and it is the tier most hosts will use. Its assumption — that a host's tRNA pools are on average matched to that host's own native codon demand — is the same empirical regularity that underlies tAI, so it is a documented approximation rather than an invention.

**The critical design decision is what tier 3 *reports*.** It would be easy, and wrong, to emit the same `relative_growth_rate` unit at all three tiers. This plan requires instead:

- `theta_unit = "fraction_of_native_demand"` and the term's unit becomes `dimensionless`;
- `Uncertainty.calibrated = false`;
- the `low_confidence_host` rationale tag is attached to every candidate;
- `TST-BURD-005` asserts all of the above, so the honesty requirement cannot erode through refactoring.

`EXP-04` tests whether the tier-3 mismatch score correlates with measured burden at all. If it does not, the term is demoted from a Pareto axis to a `diagnostic` for tier-3 hosts — which the `role` field in the job request already supports, so no redesign is needed to act on that result.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Track per-tRNA-gene demand | No decision-relevant resolution gained; codon choice cannot address individual genes. |
| Track amino-acid and energy pools | Constant under the design operator (§Decision). Would present the user with objectives that cannot be improved. |
| Require tRNA abundance data, refuse hosts without it | Restricts CodonForge to ~5 organisms and forfeits the open-source, any-host goal. |
| Use a generic "typical bacterium" tRNA profile for unknown hosts | Worse than tier 3, and dishonest: it presents another organism's biology as this one's. Tier 3 at least uses the host's *own* genome. |
| Full two-phase ledger in v1 | Needs kinetic parameters that do not exist for most hosts; `EXP-03` first. |

## Consequences

- The ledger is `|A| ≈ 40` numbers, so it is essentially free to compute and `O(1)` to update per substitution. This is what lets `burden.growth_cost` be an L1 term and participate in tier-2 coordinate ascent.
- The `capacity.estimator_tier` field propagates from `HostProfile` into every `LandscapeResponse.terms[...].provenance` and out to the user. The tier is part of the answer, not an internal detail.
- A two-phase extension can be added later without a schema change: `GrowthCondition` already carries `growth_phase`, and the ledger would become a per-phase vector behind the same interface.

## Departures

`DEP-BURD-004` — The architecture proposal presents the burden model as predicting host burden generally. This plan specifies three estimation tiers with *different output semantics*, and requires that the tier-3 output be labelled a mismatch score rather than a growth-rate prediction.
