# CLI and Library Surfaces

Satisfies: `BRIEF-2.4-a` (the CLI and library halves of "REST/CLI/library").
Owning modules: `cf-cli`, `cf-orchestrator`.

---

## 1. One implementation, three faces

```
                 ┌─────────────┐
   REST  ───────▶│             │
   CLI   ───────▶│ cf-orchestr.│──▶ cf-search / cf-landscape / cf-decision
   library ─────▶│             │
                 └─────────────┘
```

The library entry point takes a `DesignJobRequest` — the *same* object the REST body validates against and the same object the CLI builds from its flags. There is no second validation path, no second defaulting rule and no second set of objective presets. When the three surfaces drift, it is always because someone added a code path; this structure makes that structurally hard rather than merely discouraged.

`TST-CLI-001` asserts the equivalence directly: for a corpus of requests, `codonforge design --json <req>`, `POST /v1/jobs`, and `dna_design_suite.design(**req)` produce identical lockfiles and identical fronts under a fixed seed.

**Equivalence is necessary and not sufficient, and `TST-CLI-001` alone is not enough.** All three surfaces share one orchestrator, so three surfaces that all skip a tier agree perfectly and the test passes. `TST-CLI-003` therefore asserts separately that the pipeline *ran* what the budget planned: a job configured for tier 1 and tier 2 must report them `done` rather than `skipped`, and the `CF-ORCH-006` warnings naming them must be absent. Their disappearance is the evidence the wiring is real; their presence in `P1-29` was the evidence it was not.

**What a Phase 1 front is.** `front` and `explain` presuppose a populated Pareto front, and the population loop that produces one is `P2-02`/`P2-03`. Phase 1 emits a schema-valid `pareto_front` document with the candidates the assembled tier-1/tier-2 pipeline produced and `reduction: "none"`. The document shape, the objective order and the objective senses are real and end-to-end; the multi-objective content is not, and a single-objective job yields a front of one. `DEP-CLI-003` records this so the verb is not mistaken for the full capability.

---

## 2. Library

```python
import dna_design_suite as cf

result = cf.design(
    host="escherichia_coli:562/K-12_MG1655",
    protein="MKRIL*",                       # or cds=..., or construct=...
    objectives=["epistasis.dwell_deviation",
                "structure.local_accessibility",
                "burden.growth_cost",
                "usage.cai"],
    constraints=cf.constraints.golden_gate(enzyme="BsaI") + cf.constraints.vendor("twist"),
    budget="standard",
    seed=20260825,
)

result.front              # ParetoFront, schema-validated
result.lockfile           # LandscapeLockfile
result.best(by="burden.growth_cost")
result.to_fasta("out.fa")
result.explain("cand_a1b2c3d4")   # rendered per-position rationale
```

Three properties the library must have and the REST surface cannot provide:

- **Synchronous by default.** A researcher in a notebook wants a blocking call. `design(..., background=True)` returns a handle (`job_id`, `done()`, `result(timeout)`) and the job runs on a worker pool (`DEP-CLI-006`). The call does not run the job on the caller's thread. `DEP-CLI-004` records the Phase 1 refusal. Seed import is `DEP-API-003`. An unknown route stays `CF-CLI-003`.
- **In-process execution.** No server, no network, no daemon. `cf.design` drives the same Rust core directly through `dna_design_suite._core`.
- **Every intermediate is reachable.** `cf.landscape.score(...)`, `cf.constraints.analyse(...)`, `cf.search.tier2(...)` are public. A methods paper needs to run one layer in isolation, and forcing that through a REST API would be absurd.

---

## 3. CLI

```
codonforge design      --host H (--protein P | --cds F | --construct F) [options]
                       [--seed N] [--budget fast|standard|exact]
                       [--checkpoint-dir D] [--resume]
codonforge score       --host H --cds F --terms T,...
codonforge feasibility --construct F --constraints F
codonforge seeds import --tool lineardesign --target-protein P --input F
codonforge front       --job JOB [--reduction none|knee|clustering] [--format json|fasta|csv]
codonforge explain     --job JOB --candidate CAND
codonforge replay      --lockfile F
codonforge observe submit --input F
codonforge hosts list [--supports-term T]
codonforge models list [--host H] [--kind K]
codonforge models show ARTIFACT_ID
codonforge bench run   --suite SUITE [--baseline lineardesign,cdsfold]
codonforge serve       [--port 8080]
```

### 3.1 Conventions

| Convention | Rationale |
| --- | --- |
| `--json <file>` accepted by every verb, and equals the REST body | Makes "reproduce what the UI did" a copy-paste, and makes the CLI scriptable without learning a second flag grammar. |
| `--output json` on every verb | The default human output is a table; `--output json` emits the schema-validated document. Pipelines never parse the table. |
| Exit codes: `0` success, `1` design failed, `2` invalid input, `3` infeasible constraints, `4` quota/budget, `5` internal | An infeasible constraint set is not a crash and must be distinguishable in a script. |
| Every error prints the same `CF-*` code the API returns | One vocabulary across surfaces; searching the code finds the same documentation. |
| `--seed` and `--lockfile` on `design` | Reproducibility is a flag away, not a research project. |
| `--checkpoint-dir` and `--resume` on `design` | `DEP-CLI-001`. `P1-29` delivers checkpoint/resume and this document had no way to reach it, so the capability existed with no user-facing surface. `--resume` continues from the checkpoint in `--checkpoint-dir` for the same job id, refusing a checkpoint written under a different lockfile rather than silently changing which artifacts produced the numbers. |
| `SIGINT` cancels at the next segment boundary | `DEP-CLI-002`. `P1-29`'s cancellation is cooperative and checked between segments; a blocking CLI call had no specified way to request it. The first `SIGINT` sets the flag and the job stops cleanly with a checkpoint, exiting `0` with `state: cancelled`. A second `SIGINT` aborts immediately, because a user pressing it twice has stopped asking politely — the checkpoint from the last completed segment is still on disk, so nothing is lost that was not already unwritten. |

### 3.2 `codonforge replay`

```
codonforge replay --lockfile job_xyz.lock.json --assert-identical
```

Re-runs a job from its lockfile and, with `--assert-identical`, diffs the resulting front against the original, exiting non-zero on any difference. This is both a user-facing reproducibility tool and the mechanism behind `TST-ORCH-003`.

### 3.3 `codonforge explain`

Renders a candidate's per-position rationale as an aligned text or HTML view: the codon sequence, the diff against the input, and each rationale annotation positioned under its span. It is the CLI face of the closed tag vocabulary in [`../02-schemas/rationale_tags.json`](../02-schemas/rationale_tags.json).

```
  pos   0    1    2    3    4    5
        ATG  AAA  CGC  ATC  CTG  TAA
   in   ATG  AAA  CGT  ATT  CTG  TAA
             ^    ^^^  ^^^
                  │    └── rare_codon_replaced (usage.cai +0.08, conf 0.95)
                  └─────── restriction_site_removed (BsaI, golden_gate_bsai)
  [0,6) structure_relaxed_initiation  (structure.local_accessibility +0.11, conf 0.81)
```

### 3.4 SBOL 3 export

`DEP-SBOL-001`. `TST-SBOL-001`.

`dna-design-suite design --output sbol` writes sorted N-Triples. The same bytes are `DesignResult.to_sbol()` and `GET /v1/jobs/{id}/front?format=sbol`. All three call `dna_design_suite.export.sbol.to_sbol` on the front, the sequence document, the host profile, and the constraint set. The writer does not translate the CDS again and does not choose a candidate other than `front.recommended`.

The file contains the recommended coding sequence (5' flank, CDS, 3' flank), sequence features for that CDS, the flanks, `locked_codons`, and rationale annotations whose codon spans lie inside the CDS, plus the sequence document's protein string when it is present. The host profile and the constraint set are annotations in `https://dna_design_suite.org/sbol`. SBOL 3 has no codon-usage class. `sbol:Constraint` means a topological relation between features, so the export creates none.

pySBOL3 (distribution name `sbol3`) is the optional extra `sbol`. A missing install is `CF-SBOL-001`. `Document.validate()` must report no errors and no warnings or the export raises `CF-SBOL-002` and returns no file. There is no SBOL import and no round trip. FASTA and CSV still list every candidate on the front.

---

## 4. Configuration precedence

```
CLI flag  >  --json file  >  CODONFORGE_* env var  >  ~/.config/codonforge/config.toml  >  built-in default
```

Whatever resolution occurs is echoed into the lockfile's `environment` block, so a result never depends on an unrecorded machine-local setting.

---

## 5. Distribution

- `pip install codonforge` — pure wheels with the compiled core, manylinux/macOS/Windows, no Rust toolchain required for users.
- `pip install codonforge[gpu]` — adds the torch stack for the LM term.
- `pip install codonforge[nonfree]` — **does not exist by design.** ViennaRNA and LinearFold are not redistributable under Apache-2.0 terms (`RISK-DEP-001`), so the alternative structure engines are enabled by the user installing the upstream package themselves and DNA Design Suite detecting it at runtime. Shipping an extra that pulls them in would make DNA Design Suite's own wheel non-redistributable.
- A container image with the default (permissive) engines only.
