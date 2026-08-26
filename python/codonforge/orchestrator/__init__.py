"""cf-orchestrator — support layer.

Job lifecycle, landscape lockfile pinning at submission time, worker pool management,
checkpoint/resume, provenance recording, budget tier selection.

Declared dependencies: cf-types, cf-ingest, cf-landscape, cf-search, cf-constraints, cf-
decision, cf-modelreg.
Contract: planning/02-schemas/job_request.schema.json
"""
