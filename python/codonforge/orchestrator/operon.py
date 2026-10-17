"""Joint search of two or more cistrons on one demand ledger.

Each cistron is designed on its own. The coupling is the flux-weighted codon
demand, scored by the same `burden.growth_cost` a single sequence uses. A
later pass accepts a synonymous swap only when that shared cost falls and the
sequence stays feasible. Fluxes come from the caller. This module does not
invent expression rates, gene-copy numbers, or s-values.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

from codonforge import _core
from codonforge.modelreg import ArtifactStore

from .design import DesignError, DesignRequest, DesignResult, design

_EMPTY_CONSTRAINTS = {
    "schema_version": "1.0.0",
    "constraint_set_id": "cst_none",
    "constraints": [],
}
_BASE = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "U": 3,
    "a": 0,
    "c": 1,
    "g": 2,
    "t": 3,
    "u": 3,
}


@dataclass(frozen=True)
class Cistron:
    """One coding sequence and the expression flux the caller assigns it."""

    cds: str
    flux: float
    name: str = ""


@dataclass
class OperonResult:
    """Sequences after the joint passes, and the costs those sequences have."""

    sequences: list[str]
    proteins: list[str]
    fluxes: list[float]
    joint_growth_cost: float
    solo_growth_costs: list[float]
    independent_sequences: list[str]
    independent_joint_growth_cost: float
    swaps_per_pass: list[int]
    jobs: list[DesignResult] = field(default_factory=list)


def codon_index(codon: str) -> int:
    """Slot in the 64-vector the burden term uses. A=0, C=1, G=2, T=3."""
    index = 0
    for base in codon:
        try:
            value = _BASE[base]
        except KeyError as exc:
            raise DesignError("CF-SEQ-001", f"{codon} is not a codon", exit_code=2) from exc
        index = (index << 2) | value
    if len(codon) != 3:
        raise DesignError("CF-SEQ-001", f"{codon} is not a codon", exit_code=2)
    return index


def codon_counts(cds: str, translation_table: int) -> list[float]:
    """Sense-codon counts. Stops are omitted, matching `GrowthCostTerm`."""
    protein = _core.translate(cds, translation_table)
    if len(cds) != len(protein) * 3:
        raise DesignError(
            "CF-SEQ-001",
            f"translation length {len(protein)} does not match {len(cds)} nt",
            exit_code=2,
        )
    slots = [0.0] * 64
    for index, amino in enumerate(protein):
        if amino == "*":
            continue
        codon = cds[index * 3 : (index + 1) * 3].upper().replace("U", "T")
        slots[codon_index(codon)] += 1.0
    return slots


def _weighted_counts(
    sequences: list[str], fluxes: list[float], translation_table: int
) -> list[float]:
    """Flux-weighted codon demand across the operon."""
    total = [0.0] * 64
    for cds, flux in zip(sequences, fluxes, strict=True):
        for slot, count in enumerate(codon_counts(cds, translation_table)):
            total[slot] += flux * count
    return total


def _growth(host_profile: dict, translation_table: int, counts: list[float]) -> float:
    payload = json.dumps(
        {
            "translation_table": translation_table,
            "counts": counts,
            "host_profile": host_profile,
        }
    )
    return float(json.loads(_core.burden_of_demand(payload))["growth_cost"])


def _ledger_counts(
    sequences: list[str],
    fluxes: list[float],
    translation_table: int,
    previous: list[float] | None = None,
    swap: tuple[int, str, str] | None = None,
) -> list[float]:
    """Flux-weighted codon demand.

    A one-codon synonym edits two slots of the previous vector. That edit
    matches `_weighted_counts` on the rewritten sequences, which is the
    recount a joint pass would otherwise repeat at every synonym.
    """
    if previous is None or swap is None:
        return _weighted_counts(sequences, fluxes, translation_table)
    index, old, new = swap
    counts = list(previous)
    # Stops are omitted, the same rule `codon_counts` uses. A sense codon
    # and its synonym are both counted; the caller has already kept the protein.
    if _core.translate(old, translation_table) == "*":
        return counts
    flux = fluxes[index]
    counts[codon_index(old)] -= flux
    counts[codon_index(new)] += flux
    return counts


def _joint_cost(
    sequences: list[str],
    fluxes: list[float],
    translation_table: int,
    host_profile: dict,
) -> float:
    return _growth(
        host_profile,
        translation_table,
        _weighted_counts(sequences, fluxes, translation_table),
    )


def _feasible(
    cds: str,
    constraints: dict,
    flank5: str,
    flank3: str,
) -> bool:
    construct = flank5 + cds + flank3
    start = len(flank5)
    report = json.loads(
        _core.constraint_violations(
            json.dumps(constraints),
            construct,
            [(start, start + len(cds))],
            "assembled",
        )
    )
    if any(item["actionable"] for item in report["violations"]):
        return False
    from codonforge.orchestrator import local_structure

    if not local_structure.rules(constraints):
        return True
    return local_structure.blocking_distance(construct, constraints, {}) <= 0.0


def design_operon(
    cistrons: list[Cistron],
    store: ArtifactStore,
    *,
    host: str,
    host_profile: dict,
    constraints: dict | None = None,
    objectives: list[str] | None = None,
    translation_table: int = 11,
    budget: str = "fast",
    seed: int = 0,
    generations: int = 0,
    passes: int = 2,
    five_prime_flank: str = "",
    three_prime_flank: str = "",
    host_digest: str = "sha256:" + "0" * 64,
) -> OperonResult:
    """Design each cistron, then alternate synonymous swaps on the shared ledger.

    The per-cistron jobs use `objectives`. The joint passes then keep a swap
    only when `burden.growth_cost` of the flux-weighted counts falls. A swap
    that breaks a constraint or the protein is refused.
    """
    if len(cistrons) < 2:
        raise DesignError(
            "CF-CLI-001",
            "design_operon needs at least two cistrons; a single CDS is design()",
            exit_code=2,
        )
    if passes < 1:
        raise DesignError("CF-CLI-001", "passes must be at least 1", exit_code=2)
    rules = constraints if constraints is not None else dict(_EMPTY_CONSTRAINTS)
    rules = json.loads(json.dumps(rules))
    if "constraints" not in rules:
        rules["constraints"] = []
    requested = list(objectives) if objectives is not None else ["burden.growth_cost"]
    fluxes: list[float] = []
    for index, cistron in enumerate(cistrons):
        name = cistron.name or f"cistron {index}"
        if isinstance(cistron.flux, bool) or not isinstance(cistron.flux, (int, float)):
            raise DesignError(
                "CF-CLI-001",
                f"{name}: flux must be a positive finite number supplied by the caller",
                exit_code=2,
            )
        flux = float(cistron.flux)
        if flux <= 0.0 or flux != flux or flux == float("inf"):
            raise DesignError(
                "CF-CLI-001",
                f"{name}: flux must be a positive finite number supplied by the caller",
                exit_code=2,
            )
        fluxes.append(flux)

    jobs: list[DesignResult] = []
    sequences: list[str] = []
    proteins: list[str] = []
    for cistron in cistrons:
        protein = _core.translate(cistron.cds, translation_table)
        result = design(
            DesignRequest(
                host=host,
                cds=cistron.cds,
                host_profile=host_profile,
                objectives=requested,
                constraints=rules,
                translation_table=translation_table,
                budget=budget,
                seed=seed,
                generations=generations,
                five_prime_flank=five_prime_flank,
                three_prime_flank=three_prime_flank,
                host_digest=host_digest,
            ),
            store,
        )
        if result.front is None:
            raise DesignError("CF-CLI-002", "a cistron design produced no front", exit_code=1)
        designed = result.front["candidates"][0]["sequences"][0]["cds"]
        designed_protein = _core.translate(designed, translation_table)
        if designed_protein != protein:
            raise DesignError(
                "CF-SEQ-004",
                "per-cistron design changed the protein",
                exit_code=1,
            )
        jobs.append(result)
        sequences.append(designed)
        proteins.append(protein)

    independent = list(sequences)
    independent_cost = _joint_cost(independent, fluxes, translation_table, host_profile)
    ledger = _ledger_counts(sequences, fluxes, translation_table)
    swaps_per_pass: list[int] = []
    for _pass in range(passes):
        accepted = 0
        for index, cds in enumerate(sequences):
            for position in range(len(cds) // 3):
                codon = cds[position * 3 : (position + 1) * 3]
                options = sorted(_core.synonyms(codon, translation_table))
                best = cds
                best_ledger = ledger
                best_cost = _growth(host_profile, translation_table, ledger)
                for synonym in options:
                    if synonym == codon:
                        continue
                    trial = cds[: position * 3] + synonym + cds[(position + 1) * 3 :]
                    if _core.translate(trial, translation_table) != proteins[index]:
                        continue
                    if not _feasible(trial, rules, five_prime_flank, three_prime_flank):
                        continue
                    trial_sequences = list(sequences)
                    trial_sequences[index] = trial
                    trial_ledger = _ledger_counts(
                        trial_sequences,
                        fluxes,
                        translation_table,
                        ledger,
                        (index, codon, synonym),
                    )
                    trial_cost = _growth(host_profile, translation_table, trial_ledger)
                    if trial_cost < best_cost or (
                        trial_cost == best_cost and best != cds and trial < best
                    ):
                        best = trial
                        best_ledger = trial_ledger
                        best_cost = trial_cost
                if best != cds:
                    sequences[index] = best
                    cds = best
                    ledger = best_ledger
                    accepted += 1
        swaps_per_pass.append(accepted)

    solos = [
        _growth(
            host_profile,
            translation_table,
            _weighted_counts([sequence], [flux], translation_table),
        )
        for sequence, flux in zip(sequences, fluxes, strict=True)
    ]
    return OperonResult(
        sequences=sequences,
        proteins=proteins,
        fluxes=fluxes,
        joint_growth_cost=_joint_cost(sequences, fluxes, translation_table, host_profile),
        solo_growth_costs=solos,
        independent_sequences=independent,
        independent_joint_growth_cost=independent_cost,
        swaps_per_pass=swaps_per_pass,
        jobs=jobs,
    )
