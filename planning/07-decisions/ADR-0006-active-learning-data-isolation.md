# ADR-0006 — Active-Learning Data Isolation and the Community Pool

- **Status:** Accepted
- **Resolves:** open question 3.5 (`BRIEF-3.5`)
- **Related:** [`../04-api/03-active-learning-ingestion.md`](../04-api/03-active-learning-ingestion.md), [`../02-schemas/model_manifest.schema.json`](../02-schemas/model_manifest.schema.json)

## Question as posed

> How is a lab's private wet-lab data kept from leaking into another user's host/protein-family model when both share the same base architecture, while still allowing an opt-in community-pooled baseline?

## Decision

**Tenant data never flows into a shared weight set. Isolation is structural — separate artifacts, separate stores, adapter-only training — and it is verified mechanically by a test that reads manifests, not by policy.**

### 1. Model lineage

```
   public base artifact                     (trained on public data only)
            │
            ├── tenant A adapter            (trained on tenant A data only)
            ├── tenant B adapter            (trained on tenant B data only)
            └── community-pooled artifact   (trained on explicitly-donated data)
```

- A **tenant adapter** is a separate artifact: a low-rank/residual correction plus tenant intercepts and batch effects. It contains no base weights and is never merged into the base.
- The base artifact is **never** retrained on tenant data. Not with differential privacy, not with aggregation — not at all in v1. This is a stronger guarantee than "privately trained", and it is far easier to verify.
- Community pooling is a **separate lineage**, trained only on observations whose `sharing.scope == "community"` with a timestamped consent record.

### 2. Four enforcement layers

| Layer | Mechanism | Failure it prevents |
| --- | --- | --- |
| **API** | Tenant is derived from the credential; no endpoint accepts a tenant id in a body or path | Cross-tenant addressing is not expressible |
| **Storage** | Tenant artifacts and observation datasets live in per-tenant, separately-credentialed prefixes | A training job cannot read another tenant's data even if it tried |
| **Manifest** | `isolation.scope` + `training.datasets[].tenant_id` on every artifact | Provenance is machine-readable, not documentary |
| **Test** | `GATE-ACT-004` / `TST-ACT-002`: for every artifact, `isolation.scope == public ⟹ no dataset has a tenant_id`, and `scope == tenant ⟹ every dataset is public or belongs to that tenant` | Erosion over time — the invariant is re-checked on every promotion |

The fourth layer is the one that matters in practice. Access controls decay as a system grows; a promotion-time assertion over machine-readable provenance does not.

### 3. Resolution consequences

- A tenant adapter sits at rung 2 of the resolution ladder, above the public host channel. A tenant with recalibrated models automatically gets them; other tenants cannot see or pin them.
- A lockfile pinning a tenant-A artifact cannot be replayed by tenant B: replay fails with `CF-RESOLVE-006`. Reproducibility does not silently degrade into cross-tenant access.
- Attempting to pin a foreign artifact returns `CF-AUTH-003`, not `404`. The difference between "does not exist" and "not yours" is itself information, and a `404` would leak the existence of other tenants' artifacts by timing and code.

### 4. The community pool

Opt-in per **observation**, never per account, and never retroactive.

```
sharing: { scope: "community", license: "...", consent_recorded_at: ..., consent_actor: ... }
```

The contribution path:

1. The submitter sets `scope: community` on a specific observation and names a license.
2. Consent is recorded with an actor and a timestamp, immutably.
3. The observation enters a review queue: sequences are screened (biosecurity screening against controlled-sequence lists; a design tool that pools user sequences must do this), and metadata is checked for identifiers.
4. Accepted observations join a public dataset with its own `dataset_id` and digest.
5. The community-pooled artifact is trained on that public dataset and is `isolation.scope: public`.

**Revocation.** A tenant may revoke consent. Revocation removes the observation from the public dataset going forward and marks derived artifacts `deprecated` with a successor scheduled. It does **not** retroactively unpublish an already-released artifact — that is not technically possible, and saying otherwise would be a false promise. The consent text states this plainly.

### 5. What about membership inference?

A tenant adapter is trained on few examples and could, in principle, memorize them. Two facts make this a bounded concern in v1: tenant adapters are never served to other tenants, and never published. If tenant model sharing is ever added (a plausible v2 feature — a consortium wanting a shared model), differential privacy becomes necessary rather than optional, and the plan should not pretend otherwise. Registered as `RISK-SEC-001`.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Federated learning across tenants | Solves a problem CodonForge does not have. Federation is for when a shared model is *required*; here, per-tenant adapters over a public base are simpler, verifiable, and give better per-tenant fit. |
| DP-SGD on a shared model trained on all tenant data | Adds a privacy budget to reason about, degrades utility, and — decisively — still requires all tenant data to reach one training process. The isolation guarantee becomes statistical rather than structural. |
| Full per-tenant fine-tunes of the base model | Better fit in principle, but a full copy of the weights per tenant makes the manifest-level isolation check much weaker (the artifact contains base weights, so provenance is harder to attest) and multiplies storage. |
| One shared model, contractual isolation only | Not a technical control. The question asks how leakage is *prevented*, not how it is prohibited. |

## Consequences

- Tenants with little data get most of their quality from the public base, which is the right default: a 12-observation adapter should not be able to overwrite a model fit on a million codons.
- A tenant cannot benefit from another tenant's private data. That is the point, and it is a real cost: the community pool exists to give an opt-in path to the shared benefit.
- Recalibration runs in a per-tenant worker context, which also makes the non-blocking property of [`../04-api/03-active-learning-ingestion.md`](../04-api/03-active-learning-ingestion.md) §2 easy to guarantee — separate queues fall out of separate credentials.

## Departures

`DEP-ACT-002` — The architecture proposal describes an active-learning loop and a community baseline without specifying the isolation boundary. This plan forbids tenant data from entering any shared weight set in v1, and makes the community pool an explicitly-consented, separately-reviewed, separately-trained lineage rather than an aggregation of private updates.
