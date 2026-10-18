# ADR-0001 — Implementation Language and Runtime Split

- **Status:** Accepted
- **Satisfies:** `BRIEF-2.1-c`
- **Date:** 2026-08
- **Supersedes:** —

## Context

The architecture proposal requires (a) a search engine that evaluates a very large number of candidate synonymous variants, (b) neural predictors that need PyTorch, (c) an open-source, contributor-friendly codebase with independently swappable layers, and (d) reproducible, auditable results.

Quantifying (a) settles most of the question. For a 2 kb CDS (666 codons, mean 2.9 synonymous options per codon) the search performs, per job:

| Tier | Move count (default budget) |
| --- | --- |
| Tier 2 block-local coordinate ascent | `666 positions × 2.9 alternatives × ~8 sweeps ≈ 1.5 × 10⁴` **objective deltas per restart**, × 32 restarts ≈ `5 × 10⁵` |
| Tier 3 NSGA-III | `200 individuals × 300 generations = 6 × 10⁴` **full evaluations**, each internally touching all 666 positions |

So the hot path is ~10⁶ decomposable-term delta evaluations and ~10⁵ full-sequence scorings. In CPython, a single delta evaluation costing ~40 interpreted bytecode operations lands at roughly 2–4 µs, giving 2–4 s per 10⁶ deltas *before* any real work — but the realistic figure is far worse, because a pure-Python delta implementation is closer to 30–60 µs, i.e. 30–60 s per 10⁶ deltas, per restart. NumPy vectorization does not rescue this: coordinate ascent is inherently sequential in its accept/reject decisions, so the arrays are length-3 to length-6 and vectorization overhead dominates.

## Decision

**A two-language split: a Rust compiled core for everything on the per-move hot path, and Python for orchestration, training, I/O and presentation. The two communicate through PyO3 with a hard rule that no call crosses the boundary more often than once per generation.**

### Runtime assignment

| Runtime | Modules | Justification |
| --- | --- | --- |
| `rust-core` (Rust 2021, compiled into `dna_design_suite._core`, abi3 wheel) | `cf-types`, `cf-metrics`, `cf-landscape-api`, `cf-constraints`, `cf-epistasis` (scoring), `cf-burden`, `cf-biophys-structure`, `cf-biophys-kinetics`, `cf-landscape`, `cf-search` | Per-move hot path, or a type/trait definition that the hot path needs. |
| `python-orchestration` (CPython ≥3.11) | `cf-ingest`, `cf-modelreg`, `cf-decision`, `cf-orchestrator`, `cf-active`, `cf-api`, `cf-cli`, `cf-bench`, `cf-testkit`, all `*/train/` targets | Called `10⁰–10²` times per job; benefits from Biopython, pandas, PyTorch, FastAPI, matplotlib. |
| `python-worker` (CPython ≥3.11 + PyTorch) | `cf-biophys-lm`, GPU batch scoring | Needs torch; isolated so a CPU-only install never imports it. |

### The batch-yield evaluation protocol

The structural consequence of the once-per-generation rule is that the search engine cannot *call* the expensive predictors. Instead `cf-search` is written as a resumable state machine that yields evaluation batches:

```
loop:
    batch  = engine.step()          # rust; returns Vec<Genotype> needing L0-only terms
    if batch.is_empty(): break
    scores = worker.score(batch)    # python/torch, one crossing, Arrow IPC
    engine.absorb(scores)           # rust
```

In-core terms (CAI, GC, k-mer, Potts epistasis, tAI/burden ledger, windowed accessibility) never leave Rust and are evaluated per move. Out-of-core terms (LM naturalness, global MFE when requested) are evaluated once per generation on a whole population.

## Alternatives considered

| Option | Why rejected |
| --- | --- |
| **Pure Python + NumPy/Numba** | Numba would cover the delta kernels but not the search engine's control flow, the constraint automaton, or the cache. Numba's `nopython` subset makes the trait-object plugin seam (required for swappability) impossible. Debugging Numba compile failures is a worse contributor experience than reading Rust. |
| **Pure Python, accept the slowdown** | Puts a 2 kb single-gene job at 6–20 h, versus a 10 min budget. Fails `BRIEF-2.3-g` outright. |
| **C++17 + pybind11** | Closest competitor, and the incumbent in this field (ViennaRNA, LinearFold, CDSfold are all C++). Rejected on: (1) build/distribution — `maturin` produces manylinux/macOS/Windows wheels from one command, whereas a C++ project needs bespoke CI per platform; (2) memory safety in a codebase whose core data structure is an index-heavy sliding window over a byte array; (3) `rayon` gives correct-by-construction data parallelism for the red-black block sweep, where a hand-rolled C++ thread pool is a recurring source of bugs. Rust FFI to the existing C/C++ engines is a solved problem (`bindgen`), so the incumbent-ecosystem argument does not force C++. |
| **Rust for everything, Python only as a CLI shell** | Loses PyTorch, Biopython, pandas, statsmodels and the entire scientific plotting ecosystem. Model *training* would have to be reimplemented; that is a much larger risk than an FFI boundary. |
| **Julia** | Solves the two-language problem in principle, but the contributor pool in synthetic biology is overwhelmingly Python, and the package's headline dependencies (PyTorch checkpoints for CodonTransformer/CodonBERT-class models) are Python-native. |
| **JAX end-to-end with a differentiable landscape** | Would make ADR-0002's L2 capability universal, but requires every predictor to live inside one autograd runtime, which directly contradicts the swappable-layer goal. See ADR-0002. |

## Consequences

**Positive**

- The performance budget in [`08-complexity-and-performance-budgets.md`](../03-algorithms/08-complexity-and-performance-budgets.md) becomes achievable on a workstation, which is what ADR-0008 targets.
- The rust/python seam coincides with the *train vs. infer* seam, so a contributor who only wants to improve a model never touches Rust.
- `cargo tree`/`cargo deny` give a machine-checkable answer to "did anyone violate the layering", which is how the anti-cycle rules in [`02-module-dependency-graph.md`](../01-modules/02-module-dependency-graph.md) are enforced.

**Negative**

- Two toolchains in CI, and contributors touching `cf-search` need Rust. Mitigation: `maturin develop` in the dev container; prebuilt wheels for all supported platforms so a *user* never needs a Rust toolchain.
- Types are defined twice (Rust structs + `.pyi` stubs). Mitigation: JSON Schema in `planning/02-schemas/` is normative; both sides are generated/checked against it by `contract_*` tests, so drift fails CI rather than surfacing at runtime.
- Debugging across the FFI boundary is harder. Mitigation: the boundary is narrow — four functions (`step`, `absorb`, `snapshot`, `resume`) — and every batch is serializable, so any failure is reproducible from a dumped batch.

**Licensing note.** The project ships **Apache-2.0** (patent grant matters here: see `RISK-DEP-002` on the LinearDesign patent filing). Neither ViennaRNA nor LinearFold is redistributable under Apache-2.0 terms, so both are *optional, user-installed, feature-gated* back ends behind the `StructureEngine` seam, never vendored. See `RISK-DEP-001` in the [risk register](../06-build/06-risk-register.md).

## Departures from the architecture proposal

| ID | Departure | Rationale |
| --- | --- | --- |
| `DEP-RT-001` | The proposal suggests "a compiled/accelerated core for the search engine **and folding predictors**". This plan additionally moves the epistatic surface *scoring*, the burden ledger and the constraint machinery into the compiled core, while leaving all *fitting* in Python. | The split that matters is per-move vs. per-job, not model vs. search. The Potts score and the tRNA ledger are evaluated as often as the search operators themselves; the pseudo-likelihood fitter runs once per host per month. |
| `DEP-RT-002` | The LM naturalness predictor is *not* in the compiled core, contrary to a naive reading of "biophysical predictors in the accelerated core". | A transformer forward pass is already a compiled kernel; wrapping it in Rust adds nothing and loses the HuggingFace ecosystem. It is instead moved out of the per-move path entirely by the batch-yield protocol. |
