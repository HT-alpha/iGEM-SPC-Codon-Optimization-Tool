"""cf-active — active-learning layer.

Wet-lab observation ingestion, validation, tenant dataset append, scheduled
recalibration jobs, acceptance gating, exploration-arm bookkeeping. Communicates with
predictors only by writing new artifacts into `cf-modelreg`; it never imports `cf-
landscape`.

Declared dependencies: cf-types, cf-ingest, cf-modelreg.
Contract: planning/02-schemas/observation.schema.json
"""
