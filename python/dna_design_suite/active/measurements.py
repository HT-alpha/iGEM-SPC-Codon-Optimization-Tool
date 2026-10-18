"""A short lab-measurement form, stored as a full observation.

Every field that is filled is present in the stored document. The OLS line
of one readout on GC fraction, when it fits, does not change codon choices
and does not train a sequence model. `DEP-ACT-004`.
"""

from __future__ import annotations

import secrets
from datetime import UTC, datetime
from typing import Any

from dna_design_suite import _core
from dna_design_suite.active.observations import (
    _CENSORING,
    _CONDITION,
    _HOST,
    _MEASURES,
    ObservationError,
)

_CROCKFORD = "0123456789ABCDEFGHJKMNPQRSTVWXYZ"
_ASSAYS = {
    "fluorescence",
    "luminescence",
    "western_blot",
    "rt_qpcr",
    "riboseq",
    "ms_proteomics",
    "growth_rate",
    "od600_auc",
    "soluble_fraction",
    "specific_activity",
    "titer",
    "flow_cytometry",
}
_ROLES = {"negative", "positive", "reference_construct", "wild_type", "empty_vector"}
_NOT_A_CODON_MODEL = (
    "Stored as an observation. The ordinary-least-squares line of one "
    "uncensored readout on GC fraction, when it fits, does not change codon "
    "choices. No sequence model is trained."
)


def codon_model_statement() -> str:
    return _NOT_A_CODON_MODEL


def build_observation(tenant: str, body: dict[str, Any]) -> dict[str, Any]:
    """Turn the short form into an observation `accept` can store."""
    if not isinstance(body, dict):
        raise ObservationError("CF-OBS-001", "measurement must be a JSON object")
    cds = str(body.get("cds") or "").upper().replace("U", "T").replace(" ", "")
    if not cds:
        raise ObservationError("CF-OBS-001", "cds is required")
    host = str(body.get("host_id") or body.get("host") or "")
    if _HOST.match(host) is None:
        raise ObservationError(
            "CF-OBS-007",
            f"host_id {host!r} must look like genus_species:taxid",
        )
    measures = str(body.get("measures") or "")
    if measures not in _MEASURES:
        raise ObservationError(
            "CF-OBS-005",
            f"measures must be one of {', '.join(sorted(_MEASURES))}",
        )
    try:
        value = float(body["value"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ObservationError("CF-OBS-001", "value must be a number") from exc
    unit = str(body.get("unit") or "")
    if not unit:
        raise ObservationError("CF-OBS-001", "unit is required")
    assay = str(body.get("assay") or "fluorescence")
    if assay not in _ASSAYS:
        raise ObservationError("CF-OBS-001", f"assay {assay!r} is not a known assay")
    condition = str(body.get("condition_id") or "cond_unspecified")
    if _CONDITION.match(condition) is None:
        raise ObservationError(
            "CF-OBS-007",
            "condition_id must match cond_ and then 1 to 48 letters, digits, or . _ -",
        )
    role = str(body.get("control_role") or "negative")
    if role not in _ROLES:
        raise ObservationError("CF-OBS-003", f"control_role {role!r} is not a known control")
    censoring = str(body.get("censoring") or "none")
    if censoring not in _CENSORING:
        raise ObservationError(
            "CF-OBS-004",
            f"censoring {censoring!r} is not none, left, right, or interval",
        )
    try:
        table = int(body.get("translation_table") or 11)
    except (TypeError, ValueError) as exc:
        raise ObservationError("CF-OBS-002", "translation_table must be an integer") from exc
    try:
        sequence = _core.sequence_document(cds, table, None, None, None)
    except ValueError as exc:
        raise ObservationError("CF-OBS-002", str(exc)) from exc
    import json

    document = json.loads(sequence)
    measured = str(body.get("measured_on") or _now())
    submitted = _now()
    notes = str(body.get("notes") or body.get("assay_note") or "")
    return {
        "schema_version": "1.0.0",
        "observation_id": _observation_id(),
        "submitted_at": submitted,
        "tenant_id": tenant,
        "subject": {"kind": "explicit", "cds": cds, "digest": document["digest"]},
        "host_id": host,
        "condition": {"condition_id": condition},
        "readouts": [
            {
                "readout_id": "ro_measurement",
                "assay": assay,
                "measures": measures,
                "value": value,
                "unit": unit,
                "n_replicates": 1,
                "normalization": {"method": "none"},
                "dynamic_range": {"censoring": censoring},
            }
        ],
        "controls": [{"control_id": "ctl_measurement", "role": role}],
        "exploration": {"is_exploration_arm": False, "selection_rule": "user_chosen"},
        "protocol": {"measured_on": measured, "notes": notes},
        "sharing": {"scope": "private"},
    }


def observations_csv(documents: list[dict[str, Any]]) -> str:
    """One row per stored observation. Fields that were filled are columns."""
    header = (
        "observation_id,cds,host_id,measures,value,unit,assay,condition_id,"
        "control_role,measured_on,censoring,notes"
    )
    rows = [header]
    for document in documents:
        subject = document.get("subject") or {}
        readout = (document.get("readouts") or [{}])[0]
        control = (document.get("controls") or [{}])[0]
        protocol = document.get("protocol") or {}
        dynamic = readout.get("dynamic_range") or {}
        rows.append(
            ",".join(
                _cell(value)
                for value in (
                    document.get("observation_id", ""),
                    subject.get("cds", ""),
                    document.get("host_id", ""),
                    readout.get("measures", ""),
                    readout.get("value", ""),
                    readout.get("unit", ""),
                    readout.get("assay", ""),
                    (document.get("condition") or {}).get("condition_id", ""),
                    control.get("role", ""),
                    protocol.get("measured_on", ""),
                    dynamic.get("censoring", ""),
                    protocol.get("notes", ""),
                )
            )
        )
    return "\n".join(rows) + "\n"


def _observation_id() -> str:
    return "obs_" + "".join(_CROCKFORD[secrets.randbelow(32)] for _ in range(26))


def _now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _cell(value: object) -> str:
    text = str(value).replace('"', '""')
    if any(mark in text for mark in ',\n\r"'):
        return f'"{text}"'
    return text
