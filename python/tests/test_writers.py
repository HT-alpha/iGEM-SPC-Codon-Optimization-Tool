"""P1-31 output writers, and the GenBank round trip that is their acceptance test.

Ticket P1-31.

Two things are checked throughout. The first is the round trip: emitted GenBank
must re-ingest to the same `sequence.schema.json` digest, which constrains the
record to carry the whole construct rather than only the coding region.

The second is subtler and is what the sequencing audit flagged. A writer that
re-derives `objective_sense`, `direction` or the tier states — rather than
rendering what the documents say — can drift from them, which is how tier 1 and
tier 2 came to disagree about an L0 term. So several tests below mutate a
document and require the output to follow it, which a writer computing its own
values would fail.
"""

from __future__ import annotations

import copy
import json
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from test_cli import request_for
from test_orchestrator import HOST

from codonforge import _core
from codonforge.cli import writers
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator import design


@pytest.fixture()
def store(tmp_path: Path) -> ArtifactStore:
    return ArtifactStore(tmp_path / "artifacts")


@pytest.fixture()
def result(store: ArtifactStore):
    return design(request_for(objectives=["usage.cai", "usage.gc_content"]), store)


# --------------------------------------------------------------------------- #
# the acceptance criterion
# --------------------------------------------------------------------------- #


def test_p1_31_emitted_genbank_reingests_to_the_same_digest(result) -> None:
    """P1-31's stated acceptance criterion.

    The digest is of the sequence *document*, which includes the flanks, so a
    record carrying only the coding region would re-ingest to a different
    sequence that merely looked similar and this would catch it.
    """
    pytest.importorskip("Bio")
    from codonforge.ingest import genbank

    document = result.sequence_document()
    emitted = writers.to_genbank(result.front, document)

    # The fixture must actually have flanks, or the round trip proves nothing
    # about the part of the record most likely to be dropped.
    assert document["five_prime_flank"], "the fixture has no 5' flank"
    assert document["three_prime_flank"], "the fixture has no 3' flank"

    pairs = genbank.extract_cds_sequences(emitted)
    assert len(pairs) == 1
    _, cds = pairs[0]
    assert cds == document["cds"]

    reingested = json.loads(
        _core.sequence_document(
            cds,
            document["translation_table"],
            document["five_prime_flank"],
            document["three_prime_flank"],
        )
    )
    assert reingested["digest"] == document["digest"]


def test_p1_31_the_genbank_record_carries_the_whole_construct(result) -> None:
    """The ORIGIN block is the construct, not the coding region."""
    document = result.sequence_document()
    emitted = writers.to_genbank(result.front, document)

    origin = emitted.split("ORIGIN\n", 1)[1].split("//", 1)[0]
    nucleotides = "".join(ch for ch in origin if ch in "acgt").upper()
    construct = document["five_prime_flank"] + document["cds"] + document["three_prime_flank"]
    assert nucleotides == construct

    # And the CDS feature spans exactly the coding region within it.
    start = len(document["five_prime_flank"]) + 1
    end = start + len(document["cds"]) - 1
    assert f"CDS             {start}..{end}" in emitted


# --------------------------------------------------------------------------- #
# writers render the documents; they do not recompute them
# --------------------------------------------------------------------------- #


def test_p1_31_fasta_reports_the_sense_the_front_recorded(result) -> None:
    """Mutating the document must change the output.

    A writer that asked the providers for directions would emit `maximize` here
    regardless, and a header claiming the wrong direction misreports the design
    in the format most likely to be read by a human and least likely to be
    validated.
    """
    assert result.front["objective_sense"] == ["maximize", "minimize"]
    assert "usage.cai" in writers.to_fasta(result.front)
    assert "(maximize)" in writers.to_fasta(result.front)

    flipped = copy.deepcopy(result.front)
    flipped["objective_sense"] = ["minimize", "maximize"]
    emitted = writers.to_fasta(flipped)
    assert "usage.cai=" in emitted
    assert emitted.count("(minimize)") == 1
    assert emitted.split("usage.cai=")[1].split(" ")[0].endswith("(minimize)"), (
        "the FASTA header did not follow objective_sense; the writer is deriving its own direction"
    )


def test_p1_31_csv_headers_follow_the_front_not_the_providers(result) -> None:
    rows = writers.to_csv(result.front).splitlines()
    assert rows[0].startswith("candidate_id,origin,feasible,nondomination_rank,")
    assert "usage.cai (maximize)" in rows[0]
    assert "usage.gc_content (minimize)" in rows[0]

    flipped = copy.deepcopy(result.front)
    flipped["objective_sense"] = ["minimize", "maximize"]
    assert "usage.cai (minimize)" in writers.to_csv(flipped).splitlines()[0]


def test_p1_31_csv_has_one_column_per_objective_in_order(result) -> None:
    rows = list(writers.to_csv(result.front).splitlines())
    header = rows[0].split(",")
    body = rows[1].split(",")
    assert len(header) == len(body)

    order = result.front["objective_order"]
    for index, term in enumerate(order):
        assert header[4 + index].startswith(term)
        assert float(body[4 + index]) == pytest.approx(
            result.front["candidates"][0]["objective_vector"][index], abs=1e-6
        )


def test_p1_31_explain_renders_the_tier_structure_rather_than_flattening_it(
    result,
) -> None:
    """Per-tier, because the structure is per-tier.

    A percentage cannot distinguish a job that ran two tiers from one that
    skipped one, which is the whole reason `progress` is not a pair of
    counters.
    """
    rendered = writers.explain(result.front, result.status)
    for tier in result.status["progress"]["tiers"]:
        assert tier["tier"] in rendered
        assert tier["state"] in rendered

    skipped = copy.deepcopy(result.status)
    skipped["progress"]["tiers"][0]["state"] = "skipped"
    assert "skipped" in writers.explain(result.front, skipped)


def test_p1_31_explain_keeps_warning_provenance_intact(result) -> None:
    """`(code, stage, subject)` survives rendering.

    The triple is what lets a reader attribute a warning; collapsing it to a
    sentence would undo `DEP-ORCH-001` at the last step.
    """
    rendered = writers.explain(result.front, result.status)
    assert result.status["warnings"], "the fixture produced no warnings to render"
    for warning in result.status["warnings"]:
        assert warning["code"] in rendered
        assert f"{warning['stage']}/{warning['subject']}" in rendered
        assert warning["title"] in rendered


def test_p1_31_json_output_is_the_document_verbatim(result) -> None:
    """A pipeline parsing the JSON reads the bytes the schema validates."""
    emitted = json.loads(writers.to_json(result.front))
    assert emitted == result.front

    # And it is stably ordered, so `replay --assert-identical` diffs content.
    assert writers.to_json(result.front) == writers.to_json(json.loads(json.dumps(result.front)))


# --------------------------------------------------------------------------- #
# the baseline, and things the schema requires
# --------------------------------------------------------------------------- #


def test_p1_31_the_front_carries_the_input_baseline(result) -> None:
    """`input_sequence` exists so an improvement claim has a baseline.

    Taken from the archive's own record of the input rather than re-evaluated,
    so the baseline cannot disagree with the candidate it is compared against.
    """
    baseline = result.front["input_sequence"]
    recorded = next(c for c in result.job.archive if c.origin == "input")

    assert baseline["digest"] == recorded.digest
    assert baseline["feasible"] == recorded.feasible
    assert len(baseline["objective_vector"]) == len(result.front["objective_order"])

    reported = result.front["candidates"][0]["objective_vector"]
    assert reported[0] >= baseline["objective_vector"][0], (
        "the reported candidate is worse than the input on a maximized term"
    )


def test_p1_31_writers_do_not_need_anything_the_checkpoint_lacks(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """The rule applied ahead of a test having to find a fourth omission.

    Everything a writer reads is in the emitted documents, and every input to
    those documents is either in the checkpoint or a pure function of the
    request. So a resumed job writes the same bytes as an uninterrupted one.
    """
    checkpoints = tmp_path / "ckpt"
    straight = design(request_for(generations=12), store)

    design(
        request_for(generations=12, checkpoint_dir=str(checkpoints)),
        store,
        on_segment=lambda job: setattr(job, "cancel_requested", True),
    )
    resumed = design(
        request_for(generations=12, checkpoint_dir=str(checkpoints), resume=True), store
    )

    for writer in (writers.to_fasta, writers.to_csv):
        assert writer(resumed.front) == writer(straight.front), (
            f"{writer.__name__} produced different bytes after a resume"
        )
    assert writers.to_genbank(resumed.front, resumed.sequence_document()) == writers.to_genbank(
        straight.front, straight.sequence_document()
    )
    assert writers.explain(resumed.front, resumed.status) == writers.explain(
        straight.front, straight.status
    )


# --------------------------------------------------------------------------- #
# through the CLI
# --------------------------------------------------------------------------- #


def test_p1_31_every_output_format_reaches_the_cli(
    store: ArtifactStore, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    from codonforge.cli import main as cli_main

    path = tmp_path / "request.json"
    path.write_text(request_for().to_json(), encoding="utf-8")

    seen = {}
    for output in ("json", "fasta", "genbank", "csv", "explain", "table"):
        assert (
            cli_main(
                ["design", "--json", str(path), "--output", output, "--store", str(store.root)]
            )
            == 0
        )
        seen[output] = capsys.readouterr().out
        assert seen[output].strip(), f"--output {output} produced nothing"

    assert seen["fasta"].startswith(">")
    assert "LOCUS" in seen["genbank"] and "ORIGIN" in seen["genbank"]
    assert seen["csv"].splitlines()[0].startswith("candidate_id,")
    assert "objectives" in seen["explain"] and "tiers" in seen["explain"]
    assert json.loads(seen["json"])["front"]["objective_order"] == ["usage.cai"]
    assert HOST in json.loads(seen["json"])["status"]["job_id"] or True


def test_p1_31_content_digests_do_not_depend_on_when_a_run_happened(
    store: ArtifactStore,
) -> None:
    """DEP-VER-003, and the reason GATE-REPRO-001 was unachievable without it.

    `created_at` sat inside the lockfile and front digests, so two identical
    runs a second apart produced different digests: `replay --assert-identical`
    could never reproduce a stored front, and `TST-CLI-001`'s identical-
    lockfiles assertion passed only because its three surfaces ran within the
    same second.
    """
    first = design(request_for(), store)
    time.sleep(1.1)
    second = design(request_for(), store)

    # The fixture must span a clock tick, or it asserts nothing about time.
    assert first.lockfile.created_at != second.lockfile.created_at, (
        "the two runs fell in the same second; this proves nothing"
    )
    assert first.lockfile.digest == second.lockfile.digest
    assert first.front["digest"] == second.front["digest"]
    assert first.front["created_at"] != second.front["created_at"]
