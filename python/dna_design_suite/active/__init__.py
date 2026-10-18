"""cf-active — active-learning layer.

`POST /v1/observations` validates the eight `CF-OBS-*` checks and appends the
document. The outcome model is fit from stored explicit CDS and a numeric
point readout, or left `not_fitted`. Recalibration is not scheduled: this
package does not import a predictor and does not write a model artifact.

Declared dependencies: cf-types, cf-ingest, cf-modelreg.
Contract: planning/02-schemas/observation.schema.json
"""

from .observations import ObservationError, ObservationLog, accept
from .outcome import fit_outcome

__all__ = [
    "ObservationError",
    "ObservationLog",
    "accept",
    "fit_outcome",
]
