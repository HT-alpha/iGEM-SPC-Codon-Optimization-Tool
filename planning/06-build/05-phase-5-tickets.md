# Phase 5 — The Active-Learning Loop

**What Phase 5 delivers:** wet-lab feedback ingestion, tenant-isolated scheduled recalibration, the community pool, and multi-tenant operation. It also delivers everything `EXP-10` — the project's capstone validation — needs in order to run.

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P5-01` | Observation ingestion | `cf-active` | `POST /v1/observations`; all eight validation checks; append-only tenant dataset with content digests | `TST-ACT-001` | Phase 4 | M |
| `P5-02` | Tenant storage partitioning | `cf-active`, `cf-modelreg` | Per-tenant, separately-credentialed artifact and dataset prefixes; the manifest isolation fields | `TST-ACT-002`, `GATE-ACT-004` | `P5-01` | M |
| `P5-03` | Recalibration scheduler | `cf-active` | Count/time/manual triggers; 1-hour debounce; per-`(tenant, kind, host)` collapsing; a separate worker queue | `TST-ACT-003` | `P5-02` | M |
| `P5-04` | Dwell-head adapters | `cf-epistasis` | Low-rank/residual correction plus tenant intercepts and batch random effects; adapter artifact packaging with `parent_artifact` lineage | Adapter improves tenant holdout fit without touching base weights | `P5-03` | L |
| `P5-05` | Capacity tier-1 refit | `cf-burden` | Fit `θ, n, g_max, p` from tenant burden/expression observations with censoring handled | A synthetic tenant panel recovers the generating capacity | `P5-03` | M |
| `P5-06` | Acceptance gating and promotion | `cf-modelreg` | All five `GATE-ACT-*` gates; gate reports; promotion only on full pass; rejection reasons surfaced | `TST-ACT-004`, `GATE-ACT-001`…`005` | `P5-04`, `P5-05` | M |
| `P5-07` | Exploration arms | `cf-active`, `cf-decision` | Arm assignment surfaced in the UI/CLI; selection-rule covariate and inverse-propensity weighting in the fitter; the frozen never-trained holdout | `GATE-ACT-001` evaluated on the frozen holdout | `P5-03` | L |
| `P5-08` | Community pool | `cf-active` | Per-observation consent with actor and timestamp; review queue; biosecurity screening against controlled-sequence lists; public dataset assembly; revocation semantics | An end-to-end donation reaches a public dataset only after review; revocation removes it going forward | `P5-02` | L |
| `P5-09` | Multi-tenant API | `cf-api` | Bearer-token tenancy; per-tenant quotas in job-seconds; `X-CodonForge-Quota-*`; `CF-AUTH-*` semantics including `403`-not-`404` for foreign artifacts | `TST-API-003` | `P5-02` | L |
| `P5-10` | `GET /v1/recalibrations` | `cf-api` | Recalibration state, gate reports, rejection reasons | `contract_*` against the endpoint schema | `P5-06` | S |
| `P5-11` | `EXP-06` support | `cf-bench` | Solubility-panel ingestion; the pre-registered mixed model; the kinetics promotion decision path | `EXP-06` decided; `GATE-KIN-001` evaluated | `P5-01` | M |
| `P5-12` | `EXP-10` support | `cf-bench` | Blinded arm generation including the ablation arm; randomized plate layout export; result ingestion; the pre-registered analysis | `EXP-10` running with its analysis pre-registered | `P5-07` | L |
| `P5-13` | Reproducibility hardening | `cf-orchestrator` | `GATE-REPRO-001` across the whole benchmark corpus; long-term artifact retention policy; lockfile format stability guarantee | `GATE-REPRO-001` green on every stored benchmark job | Phase 4 | M |

## Phase 5 exit checklist

- [ ] Phase 4 checklist still green
- [ ] `TST-ACT-*`, `TST-API-*` green
- [ ] `GATE-ACT-001`…`005`, `GATE-REPRO-001` pass
- [ ] `TST-ACT-003` demonstrates a bit-identical front under a concurrent artifact promotion — the executable proof of `BRIEF-2.4-b`
- [ ] `EXP-06` decided; the kinetics term is promoted, kept as a diagnostic, or removed, per `AMD-002`
- [ ] `EXP-10` running with a pre-registered analysis plan and a blinded ablation arm
- [ ] Community-pool consent, review and revocation paths exercised end-to-end

---

## After Phase 5

Not planned in detail here, but recorded so the design does not accidentally foreclose them:

| Candidate | Prerequisite |
| --- | --- |
| Two-phase resource ledger | `EXP-03` showing material design divergence |
| Co-translational folding simulation (kinetics v2) | `EXP-06` promoting the proxy *and* large systematic residuals remaining |
| Attention-based context model replacing the Potts surface | `RISK-EPI-002` materializing; would forfeit the concave objective and `O(W·S)` deltas, so it needs its own performance plan |
| Island-model parallel search | A cluster deployment with a validated convergence comparison |
| Cross-tenant shared models | Differential privacy, because `RISK-SEC-001` stops being bounded once tenant models are shared |
| Chaperone and secretion capacity in the ledger | `RISK-BURD-003`; needs a data source that does not currently exist |
