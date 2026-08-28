"""cf-orchestrator — job lifecycle.

Ticket P1-29. Contract: `schemas/v1/job_status.schema.json`.

Owns the order that makes the guarantees hold — seal the landscape lockfile,
plan the tiers, run in segments, emit a status — and the three things a job adds
to a bare search: a compute budget, an RNG plan taken from the sealed lockfile,
checkpoint/resume, and cooperative cancellation.

It is also where three earlier tickets' warnings finally reach a document a user
can read: `P1-23`'s operator availability, `P1-18b`'s collection of them, and
`P1-25`'s tier-2 capability precheck all land in `job_status.warnings` tagged
with the `stage` and `subject` of `DEP-ORCH-001`.

Declared dependencies: cf-types, cf-ingest, cf-landscape, cf-search,
cf-constraints, cf-decision, cf-modelreg.
"""

from .budget import ComputeBudget, TierPlan, TierPlanError, Warning_
from .job import Cancelled, Job, JobInputs

__all__ = [
    "Cancelled",
    "ComputeBudget",
    "Job",
    "JobInputs",
    "TierPlan",
    "TierPlanError",
    "Warning_",
]
