# ADR-0007 — Importing External Tool Output as Seed Candidates

- **Status:** Accepted
- **Resolves:** open question 3.6 (`BRIEF-3.6`)
- **Related:** [`../02-schemas/seed_set.schema.json`](../02-schemas/seed_set.schema.json), [`../05-testing/04-benchmark-harness.md`](../05-testing/04-benchmark-harness.md)

## Question as posed

> Should DNA Design Suite be able to import a LinearDesign/CDSfold or CodonTransformer output as a seed candidate for its population search, to make migration easier for existing users? If so, what is the interface for that?

## Decision

**Yes, and it should be built in Phase 1 rather than deferred — because the same interface is what the benchmark harness needs.**

Seed import is unusually cheap relative to its value. The whole mechanism is: parse a sequence, check it translates to the target protein, repair it to feasibility, inject it into the initial population, and record where it came from. There is no model, no training, and no algorithmic risk.

### Interface

- **Schema:** [`seed_set.schema.json`](../02-schemas/seed_set.schema.json).
- **REST:** `POST /v1/seeds:import` accepts either a structured body (`tool`, `version`, `parameters`, `payload`) or raw FASTA; returns a `SeedSet` with a per-seed acceptance verdict.
- **Job:** `job_request.seeds` carries a `SeedSet` directly.
- **CLI:** `codonforge seeds import --tool lineardesign --target-protein P --input out.fa`.
- **Library:** `cf.seeds.import_(...)`.

### Acceptance is strict, and rejections are explained

A seed is admitted only if it translates to **exactly** the target protein under the target translation table. Rejection reasons are enumerated in the schema (`protein_mismatch`, `length_mismatch`, `invalid_alphabet`, `not_whole_codons`, `internal_stop`, `duplicate`, `locked_span_conflict`), and a `protein_mismatch` carries a **positional diff**:

```
{ "accepted": false, "reason": "protein_mismatch",
  "protein_diff": [ { "position": 2, "expected": "R", "found": "H" } ] }
```

A user who pastes the wrong file must be told which residue disagrees, not that "the import failed". This is the single most likely user error in the whole feature.

### Injection policy

| Parameter | Default | Rationale |
| --- | --- | --- |
| `max_fraction_of_population` | 0.20 | A population dominated by near-identical seeds from one tool collapses diversity within a few generations and the search never explores. |
| `deduplicate` | true | Two tools frequently emit the same CAI-maximal sequence. |
| `repair_before_inject` | true | An external tool knows nothing about the user's assembly standard, so seeds routinely carry restriction sites. Repair happens *before* injection so the population starts feasible. |
| `also_inject_mutated_copies` | true | Gives crossover material to work with from generation 1. |
| `protect_for_generations` | 0 | No elitism grace period. A seed that cannot survive on merit should not be protected — and a seed that *is* good will survive anyway, which is itself the answer to "was the external tool's output any good". |

### Provenance survives to the output

`ParetoFront.candidates[].provenance.seed_source` records the originating tool, and the `seed_derived` rationale tag renders as *"Derived from a LinearDesign seed after 128 generations"*. This is what turns a migration feature into an evidence feature: a user can see directly whether DNA Design Suite improved on the tool they were using.

## Why this earns its place in Phase 1

Three uses, one mechanism:

1. **Migration.** An existing user keeps their current workflow and adds DNA Design Suite on top rather than replacing it.
2. **Benchmarking.** `BRIEF-2.5` requires MFE/CAI trade-off comparisons against LinearDesign and CDSfold. Placing a competitor's output on DNA Design Suite's own axes *is* seed import plus `/v1/score`. Without it, the benchmark needs a parallel, unvalidated ingestion path — and benchmark code that differs from product code is how benchmark claims become wrong.
3. **Warm starts.** A user re-running a design after changing one constraint can seed from the previous job's front (`tool: previous_job`), which is both faster and produces a more comparable result.

## Adapters and what they must not do

| Tool | Input format | Notes |
| --- | --- | --- |
| LinearDesign | FASTA / stdout | Record `lambda` in `source.parameters` — a point on an MFE/CAI trade-off curve is meaningless without the trade-off parameter that produced it. |
| CDSfold | FASTA | — |
| CodonTransformer / CodonBERT | FASTA or plain sequence | Record the organism token and whether sampling was deterministic. |
| DNA Chisel | GenBank | May carry annotations; only the CDS is imported. |
| IDT / GenScript / JCat | FASTA | Vendor tools; often already constraint-filtered for a different standard. |
| Native gene | FASTA / GenBank | The most important "tool" of all: the wild-type sequence is the baseline every result should be compared against. |

**Adapters operate purely at the nucleotide level.** No adapter imports a tokenizer, a model or a scoring function from the originating tool. This keeps `cf-ingest` free of a dependency on `cf-biophys-lm` and keeps the module graph acyclic; it also means adding an adapter is ~30 lines and no new dependency.

`reported_metrics` stores the originating tool's own claimed metrics *separately* from DNA Design Suite's recomputation, so a discrepancy is visible rather than hidden. Different tools compute CAI against different reference sets, and silently overwriting one with the other would produce comparison tables that neither tool's authors would recognize.

## Consequences

- `ParetoFront` candidates carry a `seed_import` origin, so "did we beat LinearDesign" is answerable from a single job's output.
- Seeds are counted against the population budget, so a large seed set reduces random exploration — hence the 20 % cap.
- No new module: seed import lives in `cf-ingest`, which already owns format parsing.

## Departures

None. The architecture proposal leaves this open; this plan answers yes and specifies the interface.
