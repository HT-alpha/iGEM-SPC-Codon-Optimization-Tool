"""PrimerForge routes: manual analysis and asynchronous search."""

from __future__ import annotations

import json
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

from codonforge.api import serve
from codonforge.modelreg import ArtifactStore
from codonforge.primers import service

TENANT = "t_localdev1"


def _call(method: str, url: str, body: dict | None = None):
    data = None if body is None else json.dumps(body).encode()
    headers = {"Authorization": f"Bearer {TENANT}", "Content-Type": "application/json"}
    request = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read())


@pytest.fixture()
def base(tmp_path: Path):
    server = serve(ArtifactStore(tmp_path / "artifacts"))
    url = f"http://127.0.0.1:{server.server_address[1]}"
    yield url
    server.shutdown()


def test_a_stable_gc_dimer_is_refused(base: str) -> None:
    oligo = "GCGCGCGCGCGCGCGCGCGC"
    status, payload = _call(
        "POST",
        f"{base}/v1/primers/analyze",
        {"forward": oligo, "reverse": oligo},
    )
    assert status == 200
    assert payload["refused"] is True
    assert payload["hetero_dimer_dG"] < -9
    assert any(flag["severity"] == "red" for flag in payload["flags"])
    assert payload["blast"] == "Specificity against genome not checked."
    assert "SantaLucia" in payload["parameter_citation"]


def test_an_amber_dg_is_flagged_and_not_refused() -> None:
    from codonforge import _core

    # The zone function is exercised in Rust. This pair is the published
    # CaPIP hetero-dimer, whose Table 2 ΔG is about −5.8, inside (−9, −4).
    payload = json.loads(_core.primer_analyze("CAGCCACTGATGCCAAGAGA", "TTGTGGAATGGCATGGCTCT"))
    assert payload["hetero_dimer_dG"] < -4
    assert payload["hetero_dimer_dG"] >= -9
    assert any(flag["severity"] == "amber" for flag in payload["flags"])
    assert payload["refused"] is False or any(
        flag["property"] != "hetero_dimer_dG" and flag["severity"] == "red"
        for flag in payload["flags"]
    )


def test_progress_is_visible_before_the_search_finishes(base: str) -> None:
    started = threading.Event()
    release = threading.Event()

    def gate() -> None:
        started.set()
        release.wait(3)

    service.PROGRESS_GATE = gate
    try:
        status, posted = _call(
            "POST",
            f"{base}/v1/primers/candidates",
            {
                "sequence": (
                    "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
                    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                    "GCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
                ),
                "region": {"start": 30, "end": 55},
                "config": {
                    "tm_target": 60.0,
                    "product_size_range": {"min": 20, "max": 90},
                    "flank_length": 40,
                    "primer_length_range": {"min": 18, "max": 22},
                },
            },
        )
        assert status == 202
        assert started.wait(2)
        code, polled = _call("GET", f"{base}/v1/primers/{posted['job_id']}")
        assert code == 200
        assert polled["state"] == "running"
        assert polled["progress"]["tiers"][0]["tier"] == "primer_generate"
        release.set()
        for _ in range(40):
            code, polled = _call("GET", f"{base}/v1/primers/{posted['job_id']}")
            if polled["state"] in {"succeeded", "failed"}:
                break
            time.sleep(0.1)
        assert polled["state"] in {"succeeded", "failed"}
        if polled["state"] == "failed":
            assert polled["problem"]["code"] == "CF-PRIM-001"
        else:
            assert polled["primer_result"]["pairs"]
            pair_id = polled["primer_result"]["pairs"][0]["pair_id"]
            code, pair = _call("GET", f"{base}/v1/primers/{posted['job_id']}/pair/{pair_id}")
            assert code == 200
            assert pair["forward_sequence"]
    finally:
        service.PROGRESS_GATE = None
        release.set()
