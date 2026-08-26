# Glossary and Notation

---

## 1. Terms that mean something specific in this package

| Term | Meaning here |
| --- | --- |
| **Epistasis / epistatic surface** | **Context-dependence of synonymous codon effects on translation** — codon-pair effects and ribosome-footprint interactions. *Not* protein-level or organismal epistasis, and not DCA's structural contacts. `AMD-005` recommends using "context-dependent codon effects" in user-facing material for exactly this reason; the `epistasis.*` namespace is retained internally to avoid churning every schema. |
| **Landscape** | The composed set of scoring terms the search queries. Not a physical fitness landscape; a queryable model ensemble. |
| **Term** | One named, namespaced scalar the landscape can produce (`usage.cai`, `burden.growth_cost`). The unit of composition, versioning and objective selection. |
| **Capability (L0/L1/L2)** | What a predictor can do for a term: score a full sequence (L0), give an exact single-substitution delta (L1), give a gradient w.r.t. relaxed codon logits (L2). See [ADR-0002](../07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md). |
| **Delta / delta table** | The exact change in a term from one synonymous substitution / from every alternative at a set of positions. The mechanism that makes tier-2 search affordable. |
| **Tier (search)** | One of the three search stages: `tier1_seed`, `tier2_block_ca`, `tier3_moea`. |
| **Tier (budget)** | One of `fast` / `standard` / `exact`, which changes *algorithmic choices*, not just parameters. [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md). |
| **Tier (capacity)** | 1, 2 or 3, describing how host capacity was estimated and therefore what the burden number *means*. [ADR-0004](../07-decisions/ADR-0004-resource-ledger-granularity.md). |
| **Lockfile** | The sealed record of every artifact digest, engine version, input digest and RNG parameter used by one job. The reproducibility primitive, and the mechanism that makes recalibration non-blocking. |
| **Artifact** | An immutable, content-addressed trained model or parameter set plus its manifest. |
| **Channel** | A mutable pointer from `(scope, kind, host)` to one artifact digest. Promotion is an atomic pointer move. |
| **Ledger** | The tRNA-isoacceptor supply/demand vector plus ribosome occupancy, shared across cistrons. |
| **Blocking vs. advisory constraint** | Blocking makes the construct unbuildable and drives constraint-domination and repair; advisory is a preference that participates in the objective vector. Conflating them is what makes a codon optimizer return nothing. |
| **Red-black partition** | The 2-colouring of contiguous codon blocks with `B ≥ W` that makes same-colour blocks conditionally independent, so they can be optimized in parallel with nothing to reconcile. |
| **Certification** | Exact (non-surrogate) re-evaluation of the reduced front before emission. No candidate is ever reported with a surrogate score. |
| **Robustness** | The probability a candidate stays non-dominated under resampling from the joint predictive distribution. Reported per candidate. |
| **Tenant** | An isolation boundary for private wet-lab data and the models derived from it. Derived from the credential, never from a request body. |

---

## 2. Mathematical notation

| Symbol | Meaning | Typical value |
| --- | --- | --- |
| `L` | Editable codons in the design | 666 (2 kb) / 5 000 (15 kb) |
| `n` | Construct length in nucleotides | ~2 400 / ~15 000 |
| `S` | Mean synonymous-family size at the designed positions | ≈2.9 |
| `Syn(a)` | Synonymous codon set of amino acid `a` | 1–8 codons (max 6 in tables 1 and 11; 8 in table 14) |
| `W` | Epistasis window half-width, in codons | 10 (default) |
| `B` | Tier-2 block size, in codons; must satisfy `B ≥ W` | 30 |
| `b` | Beam width | 64 |
| `M` | Number of objectives | 4–6 |
| `M₁` / `M₀` | L1-capable / L0-only terms | 4 / 1–2 |
| `N` | Population size | 212 |
| `G` | Generations | 300 |
| `H` | Reference directions, `C(M+p−1, p)` | 120–210 |
| `E_max` | Exact L0 evaluations per generation | 64 |
| `K` | Cistrons | 1–6 |
| `A` | tRNA isoacceptor pools | ≈40 |
| `h(x)`, `J_d(x,y)` | Potts fields and offset-`d` couplings | 64 and 64×64×`W` parameters |
| `ψ_d`, `ω`, `α_g` | Dwell-head couplings, codon main effect, gene offset | — |
| `w(x, α)` | Share of codon `x`'s decoding served by pool `α` | `Σ_α w = 1` |
| `D_α`, `S_α`, `u_α`, `a_α` | Demand, supply, utilization, availability for pool `α` | — |
| `θ`, `κ`, `n`, `g_max`, `p` | Capacity threshold, hinge sharpness, Hill exponent, asymptotic growth cost, soft-max order | 1.0, 8, 2, 0.6, 4 |
| `π_α` | Shadow price of pool `α` in the joint search | — |
| `z(i)`, `z_ref(·)`, `z_pause` | Predicted dwell z-score, native reference profile, pause threshold (99th pct) | — |

**Indexing convention.** Everything downstream of ingest is **codon-indexed**. Nucleotide indices appear only in ingest-facing schemas and are always named `*_nt`. Spans are half-open `[start, end)`. Codon index = `16b₀ + 4b₁ + b₂` with `A=0, C=1, G=2, T=3`.

**Sign convention.** All internal optimization is minimization. Terms declared `maximize` are negated once, at ingress to `cf-metrics`. The convention lives in exactly one place.

---

## 3. Identifier prefixes used throughout

| Prefix | Meaning | Defined in |
| --- | --- | --- |
| `BRIEF-*` | A requirement from the assignment brief | [`traceability.yaml`](traceability.yaml) |
| `ADR-*` | An architecture decision record | [`../07-decisions/`](../07-decisions/ADR-0001-language-and-runtime-split.md) |
| `AMD-*` | A proposed amendment to the architecture proposal | [`03-architecture-amendments.md`](03-architecture-amendments.md) |
| `DEP-*` | A local departure from the architecture, flagged in place | [`04-departure-index.md`](04-departure-index.md) |
| `EXP-*` | A named experiment with a pre-registered decision rule | [`../08-experiments/experiment-register.md`](../08-experiments/experiment-register.md) |
| `RISK-*` | A risk with a bet, falsifier, mitigation and fallback | [`../06-build/06-risk-register.md`](../06-build/06-risk-register.md) |
| `TST-*` | A specified automated test | [`../05-testing/01-test-plan.md`](../05-testing/01-test-plan.md) |
| `GATE-*` | A pass/fail acceptance threshold | [`../05-testing/02-acceptance-thresholds.md`](../05-testing/02-acceptance-thresholds.md) |
| `SYN-*` | A synthetic-data generator | [`../05-testing/03-synthetic-data-suite.md`](../05-testing/03-synthetic-data-suite.md) |
| `P{1..5}-*` | A build ticket in the given phase | [`../06-build/`](../06-build/00-sequencing-overview.md) |
| `CF-*-###` | A runtime error code returned by the API and CLI | [`../04-api/01-rest-api.md`](../04-api/01-rest-api.md) §2.6 |

Cross-references between these are checked mechanically by `validate_planning_package.py --check refs`: an `ADR-`, `AMD-`, `EXP-`, `RISK-`, `GATE-`, `TST-` or `DEP-` identifier cited anywhere must be defined in its owning document.
