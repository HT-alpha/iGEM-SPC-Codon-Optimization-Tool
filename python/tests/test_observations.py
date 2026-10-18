"""TST-ACT-001 for the eight observation checks, and the outcome model.

A check fires on a submission built to fail it. A submission with two
explicit point readouts and two GC fractions is fit by least squares on that
fraction. One row, a constant GC, or a censored assay stays `not_fitted`
and carries no coefficient.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_design_suite import _core
from dna_design_suite.active.observations import ObservationError, accept
from dna_design_suite.active.outcome import fit_outcome
from dna_design_suite.api import serve
from dna_design_suite.modelreg import ArtifactStore
from dna_design_suite.modelreg.manifest import digest_of

TENANT = "t_localdev1"
OTHER = "t_otherdev1"
_IDS = (
    "obs_01JBQ8Z4K7M2N3P5R6T7V8W9XY",
    "obs_01JBQ8Z4K7M2N3P5R6T7V8W9XX",
    "obs_01JBQ8Z4K7M2N3P5R6T7V8W9XW",
)


def _explicit(cds: str, value: float, observation_id: str, **readout_overrides: object) -> dict:
    sequence = json.loads(_core.sequence_document(cds, 11, None, None, None))
    readout = {
        "readout_id": "ro_1",
        "assay": "fluorescence",
        "measures": "protein_abundance",
        "value": value,
        "unit": "RFU",
        "n_replicates": 1,
        "normalization": {"method": "none"},
        "dynamic_range": {"censoring": "none"},
    }
    readout.update(readout_overrides)
    return {
        "schema_version": "1.0.0",
        "observation_id": observation_id,
        "submitted_at": "2026-09-30T16:20:00Z",
        "tenant_id": TENANT,
        "subject": {"kind": "explicit", "cds": cds, "digest": sequence["digest"]},
        "host_id": "escherichia_coli:511145/BL21",
        "condition": {"condition_id": "cond_lb_37"},
        "readouts": [readout],
        "controls": [{"control_id": "ctl", "role": "negative"}],
        "exploration": {"is_exploration_arm": True, "selection_rule": "random_feasible"},
        "protocol": {"measured_on": "2026-09-28T14:05:00Z"},
        "sharing": {"scope": "private"},
    }


def _code(document: dict, **kwargs: object) -> str:
    with pytest.raises(ObservationError) as caught:
        accept(document, **kwargs)
    return caught.value.code


def test_each_check_fires_and_a_good_submission_is_stored() -> None:
    good = _explicit("ATGAAATAA", 1.0, _IDS[0])
    stored = accept(good)
    assert stored["ingest_status"]["state"] == "validated"
    assert stored["digest"] == digest_of(stored)

    bad_id = _explicit("ATGAAATAA", 1.0, _IDS[0])
    bad_id["observation_id"] = "obs_not_a_ulid"
    assert _code(bad_id) == "CF-OBS-001"

    lied = _explicit("ATGAAATAA", 1.0, _IDS[0])
    lied["digest"] = "sha256:" + "ab" * 32
    assert _code(lied) == "CF-OBS-001"

    mismatched = _explicit("ATGAAATAA", 1.0, _IDS[0])
    mismatched["subject"]["digest"] = "sha256:" + "cd" * 32
    assert _code(mismatched) == "CF-OBS-002"

    missing_candidate = _explicit("ATGAAATAA", 1.0, _IDS[0])
    missing_candidate["subject"] = {
        "kind": "design",
        "job_id": "job_01JBQ8Z4K7M2N3P5R6T7V8W9XY",
        "candidate_id": "cand_a1b2c3d4",
        "as_built_digest": "sha256:" + "ab" * 32,
    }
    assert _code(missing_candidate, candidates={}) == "CF-OBS-002"

    uncontrolled = _explicit("ATGAAATAA", 1.0, _IDS[0])
    uncontrolled["controls"] = []
    assert _code(uncontrolled) == "CF-OBS-003"

    uncensor = _explicit("ATGAAATAA", 1.0, _IDS[0])
    del uncensor["readouts"][0]["dynamic_range"]
    assert _code(uncensor) == "CF-OBS-004"

    unmeasured = _explicit("ATGAAATAA", 1.0, _IDS[0])
    del unmeasured["readouts"][0]["measures"]
    assert _code(unmeasured) == "CF-OBS-005"

    unexplored = _explicit("ATGAAATAA", 1.0, _IDS[0])
    del unexplored["exploration"]
    assert _code(unexplored) == "CF-OBS-006"

    unknown_host = _explicit("ATGAAATAA", 1.0, _IDS[0])
    unknown_host["host_id"] = "Escherichia coli"
    assert _code(unknown_host) == "CF-OBS-007"

    short_replicates = _explicit("ATGAAATAA", 1.0, _IDS[0], n_replicates=2, replicate_values=[1.0])
    assert _code(short_replicates) == "CF-OBS-008"


def test_the_shipped_example_passes_when_the_candidate_is_known() -> None:
    path = (
        Path(__file__).resolve().parents[2]
        / "planning/02-schemas/examples/observation.example.json"
    )
    document = json.loads(path.read_text(encoding="utf-8"))
    subject = document["subject"]
    stored = accept(
        document,
        candidates={(subject["job_id"], subject["candidate_id"]): subject["as_built_digest"]},
    )
    assert stored["digest"] == digest_of(stored)
    assert stored["ingest_status"]["state"] == "validated"


def test_outcome_fits_gc_fraction_and_refuses_a_thin_matrix() -> None:
    low = accept(_explicit("ATGAAATAA", 1.0, _IDS[0]))
    high = accept(_explicit("ATGGGGGGG", 4.0, _IDS[1]))
    model = fit_outcome([low, high])
    assert model["state"] == "fitted"
    assert model["feature"] == "gc_fraction"
    assert model["n"] == 2
    assert model["measures"] == "protein_abundance"
    slope = model["coefficients"]["gc_fraction"]
    intercept = model["coefficients"]["intercept"]
    for document, value in ((low, 1.0), (high, 4.0)):
        cds = document["subject"]["cds"]
        fraction = sum(base in "GC" for base in cds) / len(cds)
        assert intercept + slope * fraction == pytest.approx(value)
    assert model["residual_sum_of_squares"] == pytest.approx(0.0)

    alone = fit_outcome([low])
    assert alone["state"] == "not_fitted"
    assert "coefficients" not in alone

    same = accept(_explicit("ATGAAATAA", 3.0, _IDS[2]))
    tied = fit_outcome([low, same])
    assert tied["state"] == "not_fitted"
    assert "coefficients" not in tied
    assert "rank-deficient" in tied["reason"]

    censored = accept(
        _explicit(
            "ATGGGGGGG",
            9.0,
            _IDS[2],
            dynamic_range={"censoring": "right", "upper_limit": 8.0},
        )
    )
    held = fit_outcome([low, censored])
    assert held["state"] == "not_fitted"
    assert "coefficients" not in held


def _call(method: str, url: str, headers: dict[str, str], body: dict | None = None):
    import urllib.error
    import urllib.request

    data = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read())


def _headers(tenant: str = TENANT) -> dict[str, str]:
    return {"Authorization": f"Bearer {tenant}", "Content-Type": "application/json"}


@pytest.fixture()
def base(tmp_path: Path):
    server = serve(ArtifactStore(tmp_path / "artifacts"))
    yield f"http://127.0.0.1:{server.server_address[1]}"
    server.shutdown()


def test_post_and_get_observations_are_tenant_scoped(base: str) -> None:
    first = _explicit("ATGAAATAA", 1.0, _IDS[0])
    status, payload = _call("POST", f"{base}/v1/observations", _headers(), first)
    assert status == 202
    assert payload["outcome_model"]["state"] == "not_fitted"
    assert payload["recalibration"]["state"] == "not_scheduled"
    assert "coefficients" not in payload["outcome_model"]

    second = _explicit("ATGGGGGGG", 4.0, _IDS[1])
    status, fitted = _call("POST", f"{base}/v1/observations", _headers(), second)
    assert status == 202
    assert fitted["outcome_model"]["state"] == "fitted"
    assert fitted["outcome_model"]["feature"] == "gc_fraction"
    assert fitted["outcome_model"]["n"] == 2
    assert fitted["recalibration"]["state"] == "fitted"

    status, listed = _call("GET", f"{base}/v1/observations", _headers())
    assert status == 200
    assert len(listed["observations"]) == 2
    assert listed["outcome_model"]["state"] == "fitted"

    status, other = _call("GET", f"{base}/v1/observations", _headers(OTHER))
    assert status == 200
    assert other["observations"] == []
    assert other["outcome_model"]["state"] == "not_fitted"

    status, recal = _call("GET", f"{base}/v1/recalibrations", _headers())
    assert status == 200
    assert recal["recalibrations"]
    assert recal["recalibrations"][0]["digest"] == fitted["recalibration"]["digest"]
    assert recal["recalibrations"][0]["model"] == "phase5_ols"

    status, hosts = _call("GET", f"{base}/v1/hosts", _headers())
    assert status == 200
    assert hosts["hosts"]
    assert all(item["installed"] for item in hosts["hosts"])
