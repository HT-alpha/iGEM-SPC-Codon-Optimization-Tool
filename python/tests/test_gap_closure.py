"""The tool-gap closures: seeds, usage import, measurements, feasibility, Nussinov.

Each assertion fails if the behaviour it names is removed. Dwell stays refused.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from dna_design_suite import _core
from dna_design_suite.active.measurements import build_observation, codon_model_statement
from dna_design_suite.active.observations import accept
from dna_design_suite.api.server import (
    analyse_feasibility,
    import_host_body,
    import_seeds_body,
    submit_measurement,
)
from dna_design_suite.codon.orchestrator.design import DesignError, DesignRequest, design
from dna_design_suite.codon.predictors.nussinov import NussinovEngine, paired_fraction_penalty
from dna_design_suite.codon.predictors.structure import select_engine
from dna_design_suite.ingest.usage_table import parse_usage_table
from dna_design_suite.modelreg import ArtifactStore

ROOT = Path(__file__).resolve().parents[2]
TENANT = "t_localdemo01"


class _Api:
    def __init__(self, store: ArtifactStore) -> None:
        self.store = store

    def submit_observation(self, tenant: str, body: dict) -> tuple[int, dict]:
        from dna_design_suite.api.server import Api

        return Api(self.store).submit_observation(tenant, body)


def test_seed_import_accepts_a_synonym_and_names_a_protein_diff() -> None:
    native = "ATGAAATAA"
    status, document = import_seeds_body(
        {"tool": "manual", "payload": "ATGAAGTAA", "target_cds": native, "translation_table": 11}
    )
    assert status == 200
    assert document["seeds"][0]["acceptance"]["accepted"] is True
    status, rejected = import_seeds_body(
        {"tool": "manual", "payload": "ATGTTTTAA", "target_cds": native, "translation_table": 11}
    )
    assert status == 200
    acceptance = rejected["seeds"][0]["acceptance"]
    assert acceptance["accepted"] is False
    assert acceptance["reason"] == "protein_mismatch"
    assert acceptance["protein_diff"][0]["position"] == 1
    empty, problem = import_seeds_body({})
    assert empty == 400
    assert problem["code"] == "CF-CLI-001"
    unknown, problem = import_seeds_body(
        {"tool": "not_a_tool", "payload": native, "target_cds": native}
    )
    assert unknown == 400
    assert problem["code"] == "CF-CLI-001"


def test_starting_cds_replaces_the_search_input_only_when_the_protein_matches(
    tmp_path: Path,
) -> None:
    profile = json.loads(
        (ROOT / "examples" / "kazusa_e_coli_316407.json").read_text(encoding="utf-8")
    )
    request = DesignRequest(
        host="escherichia_coli:511145/BL21",
        cds="ATGAAATAA",
        starting_cds="ATGAAGTAA",
        host_profile=profile,
        objectives=["usage.cai"],
        constraints={"schema_version": "1.0.0", "constraint_set_id": "cst_none", "constraints": []},
        budget="fast",
        generations=0,
        seed=1,
    )
    bare = DesignRequest(
        host=request.host,
        cds=request.cds,
        host_profile=profile,
        objectives=request.objectives,
        constraints=request.constraints,
        budget="fast",
        generations=0,
        seed=1,
    )
    assert "starting_cds" not in bare.identity()
    result = design(request, ArtifactStore(tmp_path / "artifacts"))
    assert result.job.state == "succeeded"
    designed = result.front["candidates"][0]["sequences"][0]["cds"]
    assert _core.translate(designed, 11) == _core.translate("ATGAAATAA", 11)
    from dataclasses import replace

    try:
        design(replace(request, starting_cds="ATGTTTTAA"), ArtifactStore(tmp_path / "bad"))
    except DesignError as exc:
        assert exc.code == "CF-SEQ-004"
        assert "1:K->F" in exc.detail
    else:
        raise AssertionError("a different protein must be CF-SEQ-004")


def test_usage_import_is_user_supplied_and_genome_wide_is_refused() -> None:
    text = "UUU 22.0( 10) UUC 16.0( 5)\n"
    profile = parse_usage_table(text, organism="Escherichia coli")
    assert profile["codon_usage"]["subset"] == "user_supplied"
    assert profile["codon_usage"]["table"]["TTT"]["relative_adaptiveness"] == 1.0
    assert profile["codon_usage"]["table"]["TTC"]["relative_adaptiveness"] == 0.5
    assert "F" not in profile["unobserved_families"]
    assert profile["cai_usable"] is False
    status, refused = import_host_body({"text": text, "label": "genome_wide"})
    assert status == 422
    assert refused["code"] == "CF-HOST-003"
    cds = parse_usage_table(">a\nATGAAATAA\n", source_format="cds")
    assert cds["codon_usage"]["subset"] == "user_supplied"
    assert cds["codon_usage"]["table"]["AAA"]["relative_adaptiveness"] == 1.0


def test_a_measurement_is_stored_and_does_not_claim_a_codon_model(tmp_path: Path) -> None:
    store = ArtifactStore(tmp_path / "artifacts")
    status, payload = submit_measurement(
        _Api(store),
        TENANT,
        {
            "cds": "ATGAAATAA",
            "host_id": "escherichia_coli:511145/BL21",
            "measures": "protein_abundance",
            "value": 12.5,
            "unit": "RFU",
            "assay": "fluorescence",
            "condition_id": "cond_lb_37",
            "control_role": "negative",
            "censoring": "none",
            "measured_on": "2026-10-18T00:00:00Z",
        },
    )
    assert status == 202, payload
    assert payload["codon_model"] == "not_trained"
    assert "does not change codon choices" in payload["statement"]
    assert payload["statement"] == codon_model_statement()
    stored = payload["observation"]
    assert stored["subject"]["cds"] == "ATGAAATAA"
    assert stored["readouts"][0]["value"] == 12.5
    assert stored["ingest_status"]["state"] == "validated"
    built = build_observation(
        TENANT,
        {
            "cds": "ATGAAATAA",
            "host_id": "escherichia_coli:511145/BL21",
            "measures": "protein_abundance",
            "value": 1,
            "unit": "RFU",
        },
    )
    accept(built)
    assert "trained" not in json.dumps(stored)


def test_nussinov_is_opt_in_and_is_not_a_free_energy(monkeypatch) -> None:
    import dna_design_suite.codon.predictors.structure as seam

    monkeypatch.setattr(seam, "available_engines", lambda: [])
    try:
        seam.select_engine()
    except Exception as exc:
        assert exc.code == "CF-STRU-003"
    else:
        raise AssertionError("a missing ViennaRNA install must still be CF-STRU-003")
    engine = select_engine("nussinov")
    assert isinstance(engine, NussinovEngine)
    identity = engine.identity()
    assert identity.name == "nussinov"
    assert identity.parameter_set == "unit_pair_AU_GC_GU_min_hairpin_3"
    values = list(engine.local_accessibility("AAAGGGGCCCUUU", 4, 8, 37.0))
    assert len(values) == 10
    assert all(0.0 <= value <= 1.0 for value in values)
    penalty = paired_fraction_penalty(values)
    assert 0.0 <= penalty <= 1.0
    assert not hasattr(engine, "global_mfe")


def test_feasibility_witness_fields_are_on_the_route() -> None:
    status, report = analyse_feasibility(
        {
            "cds": "ATGAAACGTATTCTGTAA",
            "constraints": {
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_none",
                "constraints": [],
            },
        }
    )
    assert status == 200
    assert report["satisfiability"] in {"proven_satisfiable", "proven_unsatisfiable", "unknown"}
    assert "necessary_conditions_only" in report
    assert "witness" in report


def test_dwell_is_still_refused(tmp_path: Path) -> None:
    profile = json.loads(
        (ROOT / "examples" / "kazusa_e_coli_316407.json").read_text(encoding="utf-8")
    )
    try:
        design(
            DesignRequest(
                host="escherichia_coli:511145/BL21",
                cds="ATGAAATAA",
                host_profile=profile,
                objectives=["epistasis.dwell_deviation"],
                budget="fast",
                generations=0,
            ),
            ArtifactStore(tmp_path / "artifacts"),
        )
    except DesignError as exc:
        assert exc.code == "CF-LAND-001"
    else:
        raise AssertionError("dwell must stay CF-LAND-001")


def test_rpsl_assembly_cai_does_not_fall(tmp_path: Path) -> None:
    cds = "".join((ROOT / "examples" / "rpsL.cds").read_text(encoding="utf-8").split()).upper()
    profile = json.loads(
        (ROOT / "examples" / "mg1655_ribosomal_usage.json").read_text(encoding="utf-8")
    )
    constraints = json.loads(
        (ROOT / "benchmarks" / "validation" / "constraints_assembly.json").read_text(
            encoding="utf-8"
        )
    )
    result = design(
        DesignRequest(
            host="escherichia_coli:511145/BL21",
            cds=cds,
            host_profile=profile,
            objectives=["usage.cai"],
            constraints=constraints,
            budget="fast",
            seed=7,
            generations=0,
            host_digest="sha256:" + "ab" * 32,
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )
    assert result.job.state == "succeeded", result.status.get("problem")
    vector = result.front["candidates"][0]["objective_vector"]
    assert vector[0] == 0.9866866273575934
    designed = result.front["candidates"][0]["sequences"][0]["cds"]
    assert hashlib.sha256(designed.encode()).hexdigest() == (
        "0624d04f50821da0804a0cf91fa59b9a6bde2ba8ef6156bb79e9680ba09eedfc"
    )
    assert _core.translate(designed, 11) == _core.translate(cds, 11)
    assert result.front["candidates"][0]["constraint_report"]["feasible"] is True
    codons = [designed[i : i + 3] for i in range(0, len(designed), 3)]
    assert codons[74] == "CAA"
    assert codons[88] == "GAT"
    cset = json.dumps(constraints)

    def blocked(index: int, codon: str, rule: str) -> None:
        trial = codons.copy()
        trial[index] = codon
        sequence = "".join(trial)
        report = json.loads(
            _core.constraint_violations(cset, sequence, [(0, len(sequence))], "assembled")
        )
        hits = [item["constraint_id"] for item in report["violations"] if item["actionable"]]
        assert rule in hits, hits

    blocked(74, "CAG", "c_psti")
    blocked(88, "GAC", "c_kmer12")
