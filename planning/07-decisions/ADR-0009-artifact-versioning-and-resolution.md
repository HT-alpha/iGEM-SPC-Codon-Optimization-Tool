# ADR-0009 — Model Artifact Versioning and Query-Time Resolution

- **Status:** Accepted
- **Satisfies:** `BRIEF-2.2-d` (versioning scheme for host-specific trained artifacts) and `BRIEF-2.4-c` (query-time version selection)
- **Operational spec:** [`../04-api/04-model-resolution.md`](../04-api/04-model-resolution.md)
- **Schemas:** [`model_manifest.schema.json`](../02-schemas/model_manifest.schema.json), [`landscape_lockfile.schema.json`](../02-schemas/landscape_lockfile.schema.json)

## Context

Four requirements land on the same construct, and they are usually solved separately and incompatibly:

1. Results must be reproducible and auditable (`BRIEF-2.2-d`).
2. Host- and dataset-specific model versions must be selected at query time (`BRIEF-2.4-c`).
3. A wet-lab submission must trigger recalibration without blocking other users' jobs (`BRIEF-2.4-b`).
4. Layers must remain independently swappable (`BRIEF-2.1`).

## Decision

**Artifacts are immutable and content-addressed. Coordinates (`artifact_id`) are for humans; digests are identity. Resolution runs once per job and its full result is sealed into a lockfile.**

### 1. Coordinate

```
cf:<kind>:<host_id | host_agnostic>:<MAJOR.MINOR.PATCH>
```

| Component | Bumped when |
| --- | --- |
| **MAJOR** | The artifact's *output changes meaning*: a new or renamed term, a changed unit, a changed sign convention, or a bumped `contract_major`. Never auto-upgraded. |
| **MINOR** | Retrained on new/additional data, or a changed architecture, with the same output meaning. Auto-upgraded on a channel; **changes results**. |
| **PATCH** | Same data, same architecture: a bug fix, a recalibration, a re-export. Auto-upgraded; may change results only within the calibration tolerance. |

**The coordinate is not an identity.** A channel can be re-pointed and a MINOR can be rebuilt from the same inputs with a different toolchain. Identity is `weights[].digest` and the manifest digest. Anything claiming reproducibility cites the digest; the lockfile records both.

`DEP-VER-001` — The architecture proposal asks for a versioning scheme so that results are reproducible. This plan makes the human-readable coordinate explicitly insufficient for reproducibility and requires a content digest everywhere the coordinate appears, because coordinate-only pinning is the standard way "reproducible" systems quietly stop being reproducible.

### 2. Manifest

Every artifact carries a [`ModelManifest`](../02-schemas/model_manifest.schema.json). The fields that exist specifically to catch known failures:

| Field | Failure it catches |
| --- | --- |
| `training.split_strategy` | `random_codon` splits leak across adjacent Ribo-seq positions and inflate held-out metrics. Recording the strategy makes the mistake visible in review rather than invisible in a number. |
| `code_version.dirty` | An artifact built from an uncommitted tree cannot be rebuilt. Public promotion refuses `dirty: true`. |
| `evaluation.baseline_comparison` | An absolute metric is not evidence. Every promotion states which baseline it beat and by how much, with a CI. |
| `isolation.scope` + `training.datasets[].tenant_id` | The mechanically-checked half of [ADR-0006](ADR-0006-active-learning-data-isolation.md). |
| `parent_artifact` | Adapter/fine-tune lineage, so a tenant adapter's base is auditable. |
| `calibration.empirical_coverage` | An uncalibrated `sigma` entering probabilistic dominance ([ADR-0005](ADR-0005-uncertainty-and-pareto-ranking.md)) corrupts rankings silently. |
| `lifecycle.status` | `yanked` artifacts stay resolvable **by digest forever** — a published paper must remain reproducible — but are never selected by a channel. |

### 3. Channels

A channel is a small document mapping `(scope, kind, host_id) → (artifact_id, digest)`, with names `stable` / `candidate` / `nightly`. Promotion is an atomic write, gated on every required `evaluation.gates` entry passing and on `code_version.dirty == false` for public channels. Channel history is retained so "which artifact was `stable` on date D" is answerable.

### 4. Resolution ladder, run once

`pin → tenant_channel → host_channel → nearest_host → base → builtin`, with a compatibility filter at every rung (matching `contract_major`, non-yanked, provides the term, visible to the tenant, sufficient capability level). Full algorithm, host-distance metric and admissibility rules in [`../04-api/04-model-resolution.md`](../04-api/04-model-resolution.md).

Resolution runs **once, when the job leaves `queued`** — not per query. A tier-2 sweep issues ~10⁶ landscape queries; resolving per query would be slow and, worse, non-deterministic if a channel moved mid-job.

### 5. The lockfile

The sealed resolution result, plus engine versions, input digests, RNG plan and environment. It is the single construct that satisfies all four requirements from the Context:

- **Reproducibility:** `codonforge replay --lockfile F --assert-identical` re-runs and diffs.
- **Query-time selection:** the ladder ran at pin time; every query carries its outcome.
- **Non-blocking recalibration:** a new artifact appearing mid-job is *invisible* to that job, because the job reads pinned digests.
- **Swappability:** the lockfile names engines and artifacts by interface + digest, so swapping an implementation is a resolution-time change, not a code change.

A replay that cannot find a pinned digest **fails loudly** (`CF-RESOLVE-006`). It never substitutes a newer artifact — a silent substitution would turn "reproduce my result" into "produce a different result without telling me".

## Alternatives considered

| Option | Why not |
| --- | --- |
| Semantic version pinning only (`^1.4`) | Cannot reproduce a result once a MINOR is republished. Familiar from package managers, and package managers added lockfiles for exactly this reason. |
| Content addressing only, no human coordinate | Reproducible but unusable: nobody can read a changelog of digests, and "which is the current E. coli model" becomes unanswerable. |
| A global "model set" version (one number for everything) | Couples unrelated artifacts: a burden-model patch would force a new epistasis version. Also makes host-specific evolution impossible. |
| Re-resolve per query with a cache | The failure mode is a job whose early generations used one model and whose later generations used another, producing an incoherent front. |

## Consequences

- Every response, front, and job status carries provenance down to the digest. This is verbose and it is the point.
- The registry must retain yanked artifacts indefinitely. Artifacts are small (the Potts surface is 164 KB; the LM checkpoint is the only large one), so retention is cheap.
- Adding a term is additive and bumps nothing, which is why `requested_terms` has no wildcard: a wildcard would make adding a term a breaking change for every existing client and every existing lockfile.
