"""SBOL 3 export. `TST-SBOL-001`, `DEP-SBOL-001`.

The file is produced by pySBOL3 and checked with ``Document.validate()``.
The assertions below call that validator themselves. A flag on the export
function would not be evidence.

Two selections are easy to get wrong and are pinned with documents that
would not notice a shortcut:

* The recommended candidate is not ``candidates[0]``. A writer that exports
  the first row would emit the other CDS.
* The protein string in the sequence document is not the translation of the
  CDS. A writer that translates would not copy ``QQ``.
"""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import pytest

from dna_design_suite.export.sbol import (
    EXPORT_SCOPE_TEXT,
    NAMESPACE,
    SbolExportError,
    to_sbol,
)

_PREDICATE = NAMESPACE + "#"
CDS = "ATGGCTTAA"
OTHER = "ATGAAATAA"
FLANK5 = "AAA"
FLANK3 = "TTT"


def _need_sbol() -> None:
    if importlib.util.find_spec("sbol3") is not None:
        return
    if os.environ.get("CF_REQUIRE_SBOL") == "1":
        pytest.fail("sbol3 is required (CF_REQUIRE_SBOL=1) and is not installed")
    pytest.skip("sbol3 is not installed")


def _front() -> dict:
    def candidate(identifier: str, cds: str, note: dict | None) -> dict:
        return {
            "candidate_id": identifier,
            "sequences": [{"cds": cds}],
            "rationale": {"annotations": [note]} if note else {"annotations": []},
        }

    return {
        "job_id": "job_sbol_fixture",
        "recommended": {"candidate_id": "cand_recommended"},
        "candidates": [
            candidate("cand_other", OTHER, None),
            candidate(
                "cand_recommended",
                CDS,
                {
                    "tag": "baseline_deviation",
                    "span": {"start": 1, "end": 2},
                    "text": "second codon",
                },
            ),
        ],
    }


def _sequence_document() -> dict:
    return {
        "cds": CDS,
        "five_prime_flank": FLANK5,
        "three_prime_flank": FLANK3,
        "translation_table": 11,
        # Not the translation of ATGGCTTAA. The exporter must copy this field.
        "protein": "QQ",
        "locked_codons": [
            {"span": {"start": 2, "end": 3}, "reason": "stop_codon"},
            {"span": {"start": 0, "end": 1}, "reason": "start_codon"},
            {
                "span": {"start": 1, "end": 2},
                "reason": "user_locked",
                "note": "affinity tag kept verbatim",
            },
            {"span": {"start": 9, "end": 10}, "reason": "other", "note": "outside the CDS"},
        ],
    }


def _host() -> dict:
    return {
        "host_id": "escherichia_coli:511145/MG1655",
        "display_name": "E. coli K-12 MG1655",
        "taxid": 511145,
        "strain": "MG1655",
        "codon_usage": {"subset": "highly_expressed", "table": {"ATG": {"aa": "M", "fraction": 1}}},
    }


def _constraints() -> dict:
    return {
        "constraint_set_id": "cst_sbol",
        "constraints": [
            {
                "id": "c_bsai",
                "class": "blocking",
                "rule": {"kind": "restriction_site", "enzyme": "BsaI", "site": "GGTCTC"},
            }
        ],
    }


def _export() -> str:
    return to_sbol(
        _front(),
        _sequence_document(),
        _host(),
        _constraints(),
        design_host="escherichia_coli:511145/BL21",
    )


def _read(text: str):
    import sbol3

    document = sbol3.Document()
    document.read_string(text, sbol3.SORTED_NTRIPLES)
    return document


def _sequence(reference):
    """pySBOL3 returns a URI after ``read_string``; the sequence is a lookup."""
    sequence = reference.lookup()
    assert sequence is not None
    return sequence


def _literal(document, subject: str, predicate: str) -> str:
    from rdflib import URIRef

    value = document.graph().value(URIRef(subject), URIRef(predicate))
    assert value is not None, predicate
    return str(value)


def test_a_missing_sbol3_install_is_cf_sbol_001(monkeypatch: pytest.MonkeyPatch) -> None:
    """The refusal has to fire while the package is installed, or it is untested."""
    real = importlib.util.find_spec

    def hidden(name: str, *args, **kwargs):
        if name == "sbol3":
            return None
        return real(name, *args, **kwargs)

    monkeypatch.setattr(importlib.util, "find_spec", hidden)
    with pytest.raises(SbolExportError) as caught:
        _export()
    assert caught.value.code == "CF-SBOL-001"
    assert "sbol3" in caught.value.detail
    assert caught.value.exit_code == 2


def test_sbol_export_validates_and_follows_the_documents() -> None:
    _need_sbol()
    import sbol3

    text = _export()
    assert text == _export()
    document = _read(text)
    report = document.validate()
    assert list(report.errors) == [], report
    assert list(report.warnings) == [], report
    assert not document.find_all(lambda obj: isinstance(obj, sbol3.Constraint))

    base = f"{NAMESPACE}/job_sbol_fixture"
    dna = document.find(f"{base}/recommended_cds")
    assert sbol3.SBO_DNA in dna.types
    dna_sequence = _sequence(dna.sequences[0])
    assert dna_sequence.elements == FLANK5 + CDS + FLANK3
    assert dna_sequence.encoding == sbol3.IUPAC_DNA_ENCODING
    # The other candidate's CDS must not be the sequence. It is candidates[0].
    assert OTHER not in dna_sequence.elements

    cds_features = [feature for feature in dna.features if sbol3.SO_CDS in feature.roles]
    assert len(cds_features) == 1
    cds_range = cds_features[0].locations[0]
    assert (cds_range.start, cds_range.end) == (len(FLANK5) + 1, len(FLANK5) + len(CDS))

    by_name = {feature.name: feature for feature in dna.features}
    assert by_name["five_prime_flank"].locations[0].start == 1
    assert by_name["five_prime_flank"].locations[0].end == len(FLANK5)
    assert not by_name["five_prime_flank"].roles
    assert "Not annotated as a UTR" in by_name["five_prime_flank"].description
    assert by_name["three_prime_flank"].locations[0].end == len(FLANK5 + CDS + FLANK3)
    start = by_name["locked_start_codon_4_6"]
    assert _SO_in(start, "https://identifiers.org/SO:0000318")
    assert (start.locations[0].start, start.locations[0].end) == (4, 6)
    stop = by_name["locked_stop_codon_10_12"]
    assert _SO_in(stop, "https://identifiers.org/SO:0000319")
    assert (stop.locations[0].start, stop.locations[0].end) == (10, 12)
    locked = by_name["locked_user_locked_7_9"]
    assert not locked.roles
    assert "affinity tag kept verbatim" in locked.description
    note = by_name["rationale_baseline_deviation_7_9"]
    assert "second codon" in note.description
    assert not note.roles
    assert "locked_other" not in by_name

    protein = document.find(f"{base}/protein")
    assert sbol3.SBO_PROTEIN in protein.types
    protein_sequence = _sequence(protein.sequences[0])
    assert protein_sequence.elements == "QQ"
    assert protein_sequence.encoding == sbol3.IUPAC_PROTEIN_ENCODING

    collection = f"{base}/design"
    assert _literal(document, collection, _PREDICATE + "export_scope") == EXPORT_SCOPE_TEXT
    assert "not a round trip" in EXPORT_SCOPE_TEXT
    assert "writes no sbol:Constraint" in EXPORT_SCOPE_TEXT
    assert _literal(document, collection, _PREDICATE + "sbol3_version") == sbol3.__version__
    assert _literal(document, collection, _PREDICATE + "design_host") == (
        "escherichia_coli:511145/BL21"
    )
    host = json.loads(_literal(document, collection, _PREDICATE + "host_profile_json"))
    constraints = json.loads(_literal(document, collection, _PREDICATE + "constraint_set_json"))
    assert host == _host()
    assert constraints == _constraints()
    assert "c_bsai" in _literal(document, collection, _PREDICATE + "constraint_summary")
    assert "outside the CDS" in _literal(document, collection, _PREDICATE + "omitted_feature")


def _SO_in(feature, role: str) -> bool:
    return role in list(feature.roles)


def test_a_sequence_document_for_a_different_cds_is_refused() -> None:
    _need_sbol()
    document = _sequence_document()
    document["cds"] = OTHER
    with pytest.raises(SbolExportError) as caught:
        to_sbol(_front(), document, _host(), _constraints())
    assert caught.value.code == "CF-SBOL-002"
    assert "does not match" in caught.value.detail


def test_a_non_acgt_base_is_not_rewritten() -> None:
    _need_sbol()
    document = _sequence_document()
    document["five_prime_flank"] = "AAN"
    with pytest.raises(SbolExportError) as caught:
        to_sbol(_front(), document, _host(), _constraints())
    assert caught.value.code == "CF-SBOL-002"
    assert "N" in caught.value.detail


def test_design_cli_and_result_write_the_same_validated_file(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    _need_sbol()
    import sbol3

    from dna_design_suite.cli import main as cli_main
    from dna_design_suite.codon.orchestrator import design
    from dna_design_suite.modelreg import ArtifactStore

    sys_path = str(Path(__file__).parent)
    import sys

    sys.path.insert(0, sys_path)
    from test_cli import request_for

    store = ArtifactStore(tmp_path / "artifacts")
    request = request_for(objectives=["usage.cai"], generations=0, budget="fast")
    result = design(request, store)
    assert result.sequence_document()["five_prime_flank"]
    assert result.sequence_document()["three_prime_flank"]
    text = result.to_sbol()
    document = _read(text)
    report = document.validate()
    assert list(report.errors) == [], report
    assert list(report.warnings) == [], report
    assert not document.find_all(lambda obj: isinstance(obj, sbol3.Constraint))

    recommended = result.front["recommended"]["candidate_id"]
    chosen = next(
        item for item in result.front["candidates"] if item["candidate_id"] == recommended
    )
    sequence_document = result.sequence_document()
    construct = (
        sequence_document["five_prime_flank"]
        + chosen["sequences"][0]["cds"]
        + sequence_document["three_prime_flank"]
    )
    dna = next(
        obj
        for obj in document.objects
        if isinstance(obj, sbol3.Component) and sbol3.SBO_DNA in obj.types
    )
    assert _sequence(dna.sequences[0]).elements == construct
    cds_feature = next(feature for feature in dna.features if sbol3.SO_CDS in feature.roles)
    flank = len(sequence_document["five_prime_flank"])
    cds = chosen["sequences"][0]["cds"]
    assert (cds_feature.locations[0].start, cds_feature.locations[0].end) == (
        flank + 1,
        flank + len(cds),
    )
    locks = {item["reason"] for item in sequence_document["locked_codons"]}
    assert "start_codon" in locks and "stop_codon" in locks
    roles = {role for feature in dna.features for role in feature.roles}
    assert "https://identifiers.org/SO:0000318" in roles
    assert "https://identifiers.org/SO:0000319" in roles
    collection = next(obj for obj in document.objects if isinstance(obj, sbol3.Collection))
    host = json.loads(_literal(document, collection.identity, _PREDICATE + "host_profile_json"))
    constraints = json.loads(
        _literal(document, collection.identity, _PREDICATE + "constraint_set_json")
    )
    assert host == result.job.inputs.host_profile
    assert constraints == result.job.inputs.constraint_set
    assert _literal(document, collection.identity, _PREDICATE + "design_host") == request.host

    path = tmp_path / "request.json"
    path.write_text(request.to_json(), encoding="utf-8")
    assert (
        cli_main(["design", "--json", str(path), "--output", "sbol", "--store", str(store.root)])
        == 0
    )
    assert capsys.readouterr().out == text
    refused = cli_main(["terms", "--output", "sbol"])
    assert refused == 2
    assert "CF-SBOL-002" in capsys.readouterr().err


def test_rest_format_sbol_uses_the_stored_design_request(tmp_path: Path) -> None:
    _need_sbol()
    import sbol3

    from dna_design_suite.api import serve
    from dna_design_suite.modelreg import ArtifactStore

    sys_path = str(Path(__file__).parent)
    import sys
    import time

    sys.path.insert(0, sys_path)
    from test_api import _body, _headers

    server = serve(ArtifactStore(tmp_path / "artifacts"))
    base = f"http://127.0.0.1:{server.server_address[1]}"
    try:
        body = _body(generations=0, budget="fast")
        status, payload = _json("POST", f"{base}/v1/jobs", _headers(), body)
        assert status == 202
        job_id = payload["job_id"]
        deadline = time.time() + 60
        finished = payload
        while finished["state"] in {"queued", "running"} and time.time() < deadline:
            time.sleep(0.05)
            _status, finished = _json("GET", f"{base}/v1/jobs/{job_id}", _headers())
        assert finished["state"] == "succeeded"
        code, kind, raw = _raw("GET", f"{base}/v1/jobs/{job_id}/front?format=sbol", _headers())
        assert code == 200, raw
        assert kind.startswith("application/n-triples")
        document = _read(raw.decode())
        report = document.validate()
        assert list(report.errors) == [], report
        assert list(report.warnings) == [], report
        dna = next(
            obj
            for obj in document.objects
            if isinstance(obj, sbol3.Component) and sbol3.SBO_DNA in obj.types
        )
        front_status, front = _json("GET", f"{base}/v1/jobs/{job_id}/front", _headers())
        assert front_status == 200
        recommended = front["recommended"]["candidate_id"]
        chosen = next(item for item in front["candidates"] if item["candidate_id"] == recommended)
        assert _sequence(dna.sequences[0]).elements == chosen["sequences"][0]["cds"]
        collection = next(obj for obj in document.objects if isinstance(obj, sbol3.Collection))
        host = json.loads(_literal(document, collection.identity, _PREDICATE + "host_profile_json"))
        constraints = json.loads(
            _literal(document, collection.identity, _PREDICATE + "constraint_set_json")
        )
        assert host == body["host_profile"]
        assert constraints == body["constraints"]
        assert not document.find_all(lambda obj: isinstance(obj, sbol3.Constraint))
        hidden, _kind, hidden_body = _raw(
            "GET",
            f"{base}/v1/jobs/{job_id}/front?format=sbol",
            _headers("t_otherdev1"),
        )
        assert hidden == 403
        assert json.loads(hidden_body)["code"] == "CF-AUTH-003"
    finally:
        server.shutdown()


def _json(method: str, url: str, headers: dict, body: dict | None = None):
    code, kind, raw = _raw(method, url, headers, body)
    return code, json.loads(raw) if raw else None


def _raw(method: str, url: str, headers: dict, body: dict | None = None):
    import urllib.error
    import urllib.request

    data = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request) as response:
            return response.status, response.headers.get("Content-Type"), response.read()
    except urllib.error.HTTPError as exc:
        return exc.code, exc.headers.get("Content-Type"), exc.read()


def test_the_page_fetches_sbol_instead_of_writing_rdf() -> None:
    page = Path(__file__).resolve().parents[1] / "dna_design_suite" / "api" / "static" / "app.js"
    text = page.read_text(encoding="utf-8")
    assert 'data-export="sbol"' in text
    assert "format=sbol" in text
    assert "sbols.org" not in text
    assert "edam:format_1207" not in text
