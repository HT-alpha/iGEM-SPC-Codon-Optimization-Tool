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

`TST-CLI-001` asserts the equivalence directly: for a corpus of requests, `codonforge design --json <req>`, `POST /v1/jobs`, and `codonforge.design(**req)` produce identical lockfiles and identical fronts under a fixed seed.

---

## 2. Library

```python
import codonforge as cf

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

- **Synchronous by default.** A researcher in a notebook wants a blocking call. `cf.design(..., background=True)` returns a handle for the async case.
- **In-process execution.** No server, no network, no daemon. `cf.design` drives the same Rust core directly through `codonforge._core`.
- **Every intermediate is reachable.** `cf.landscape.score(...)`, `cf.constraints.analyse(...)`, `cf.search.tier2(...)` are public. A methods paper needs to run one layer in isolation, and forcing that through a REST API would be absurd.

---

## 3. CLI

```
codonforge design      --host H (--protein P | --cds F | --construct F) [options]
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
- `pip install codonforge[nonfree]` — **does not exist by design.** ViennaRNA and LinearFold are not redistributable under Apache-2.0 terms (`RISK-DEP-001`), so the alternative structure engines are enabled by the user installing the upstream package themselves and CodonForge detecting it at runtime. Shipping an extra that pulls them in would make CodonForge's own wheel non-redistributable.
- A container image with the default (permissive) engines only.
