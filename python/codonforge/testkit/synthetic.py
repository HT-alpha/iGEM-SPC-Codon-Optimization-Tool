"""The synthetic suite: instances whose right answer is known in advance.

Ticket P1-32. Specification: `planning/05-testing/03-synthetic-data-suite.md`.

Phase 1 ships four of the seven generators — `SYN-01`, `SYN-04`, `SYN-05` and
`SYN-07`. `SYN-02`, `SYN-03` and `SYN-06` need the Potts fitter, the shared
resource ledger and the burden estimators, which are Phase 2's.

## Why ties are a first-class requirement here

`DEP-SYN-001`. A generator that only produces instances with one obvious right
answer cannot test tie-breaking, and tie-breaking is where determinism bugs
live. That is not a hypothesis: the tier-2 determinism fixture used distinct
weights, so the optimum was unique, the beam's tie-break was never invoked, and
an implementation keying block RNG streams on iteration order passed it. A
suite of unambiguous instances would have that blind spot at the scale
everything is regressed against, so both `SYN-01` and `SYN-04` ship degenerate
configurations alongside their distinct ones.
"""

from __future__ import annotations

import functools
import itertools
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from typing import Any

from codonforge import _core

#: The reduced codon table: each designable residue restricted to 3 synonyms.
#:
#: "Restricted" is the specification's word and it is doing real work. Under
#: table 1 only isoleucine is exactly 3-fold degenerate, so a genuine `3⁸`
#: space cannot be built from unrestricted amino acids — the residues below all
#: have at least 3 synonyms and are capped at the first 3 in a fixed order.
#: `M` and the stop are fixed, so the *design space* is over the 8 interior
#: codons and has exactly `3⁸ = 6 561` members, which is what `TST-SYN-001`
#: states.
REDUCED_TABLE: dict[str, tuple[str, ...]] = {
    "M": ("ATG",),
    "A": ("GCA", "GCC", "GCG"),
    "G": ("GGA", "GGC", "GGG"),
    "I": ("ATA", "ATC", "ATT"),
    "L": ("CTA", "CTC", "CTG"),
    "P": ("CCA", "CCC", "CCG"),
    "R": ("CGA", "CGC", "CGG"),
    "S": ("TCA", "TCC", "TCG"),
    "T": ("ACA", "ACC", "ACG"),
    "V": ("GTA", "GTC", "GTG"),
    "*": ("TAA",),
}

#: The 8 designable residues. One fixed start, one fixed stop, 8 free codons.
TOY_PROTEIN = "MAGILPRSTV*"[:1] + "AGILPRST" + "*"


# --------------------------------------------------------------------------- #
# SYN-01 — toy exhaustive design space
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ToySpace:
    """A design space small enough to enumerate completely.

    The only place in the system where a search result can be compared against
    the *true* optimum rather than against another algorithm.
    """

    protein: str
    #: Per-codon-position weight for each synonym, keyed by codon.
    weights: dict[str, float]
    #: True when several distinct sequences attain the optimum exactly.
    degenerate: bool
    #: Solved at construction; see `__post_init__`.
    _optima: tuple[str, ...] = field(default=(), repr=False, compare=False)
    _optimum_value: float = field(default=0.0, repr=False, compare=False)

    @property
    def positions(self) -> list[tuple[str, ...]]:
        return [REDUCED_TABLE[residue] for residue in self.protein]

    @property
    def size(self) -> int:
        size = 1
        for synonyms in self.positions:
            size *= len(synonyms)
        return size

    def enumerate(self) -> Iterator[str]:
        """Every sequence in the space, in a fixed order."""
        for choice in itertools.product(*self.positions):
            yield "".join(choice)

    def score(self, cds: str) -> float:
        """The additive objective. Decomposable, so tier 1's DP is exact."""
        codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
        return sum(self.weights[codon] for codon in codons)

    def __post_init__(self) -> None:
        # Solved once, at construction. The specification's own claim for this
        # generator is that a 6 561-member space "is enumerable in
        # milliseconds, so this test can run on every commit" - and it is, but
        # only if the enumeration happens once. Recomputing it per call turned
        # a caller's `all(score(x) <= optimum_value() for x in enumerate())`
        # into 43 million scores and 56 seconds, which is the fixture failing
        # its own design goal rather than a slow test.
        best: float | None = None
        found: list[str] = []
        for cds in self.enumerate():
            value = self.score(cds)
            if best is None or value > best:
                best, found = value, [cds]
            elif value == best:
                found.append(cds)
        object.__setattr__(self, "_optima", tuple(found))
        object.__setattr__(self, "_optimum_value", best)

    def optima(self) -> list[str]:
        """Every sequence attaining the maximum, not just one of them.

        Returning the full set is what makes a tie testable: a caller can
        assert that a search returned *an* optimum and that it returned the
        same one every time, which are different claims and both matter.
        """
        return list(self._optima)

    def optimum_value(self) -> float:
        return self._optimum_value


@functools.cache
def syn01_distinct() -> ToySpace:
    """The unambiguous configuration: exactly one optimum."""
    weights = {codon: 0.0 for codon in _all_codons()}
    for index, residue in enumerate(sorted(set(TOY_PROTEIN))):
        for rank, codon in enumerate(REDUCED_TABLE[residue]):
            # Distinct across the whole table, so no two assignments tie.
            weights[codon] = 1.0 + index * 0.37 + rank * 0.11

    space = ToySpace(protein=TOY_PROTEIN, weights=weights, degenerate=False)
    assert space.size == 3**8, f"design space is {space.size}, not 3^8"
    assert len(space.optima()) == 1, "the distinct configuration must have one optimum"
    return space


@functools.cache
def syn01_degenerate(tied_residues: str = "AG") -> ToySpace:
    """The degenerate configuration: several sequences attain the optimum.

    `DEP-SYN-001`. Every synonym of the named residues scores identically, so
    the number of optima is `∏ |Syn(a)|` over the tied residues — known in
    closed form from the construction rather than discovered by enumeration,
    which is what makes it a *ground truth* rather than an observation.
    """
    weights = {codon: 0.0 for codon in _all_codons()}
    for index, residue in enumerate(sorted(set(TOY_PROTEIN))):
        for rank, codon in enumerate(REDUCED_TABLE[residue]):
            if residue in tied_residues:
                # Flat: every synonym equally good, so the choice is a tie.
                weights[codon] = 5.0
            else:
                # Strictly preferred, so the space keeps structure and the tie
                # is localized rather than the whole space being flat.
                weights[codon] = 1.0 + index * 0.37 - rank * 0.11

    space = ToySpace(protein=TOY_PROTEIN, weights=weights, degenerate=True)
    predicted = 1
    for residue in TOY_PROTEIN:
        if residue in tied_residues:
            predicted *= len(REDUCED_TABLE[residue])
    found = len(space.optima())
    assert found == predicted, (
        f"the degenerate configuration has {found} optima, not the {predicted} its "
        "construction predicts; the tie is not where it was designed to be"
    )
    return space


# --------------------------------------------------------------------------- #
# SYN-04 — planted constraint violations
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PlantedViolation:
    """A sequence with a known violation and a known set of minimal fixes."""

    cds: str
    constraint_set: dict[str, Any]
    #: Every single-codon substitution that removes the violation, as
    #: `(codon_index, codon)`. More than one entry is a deliberate tie.
    minimal_fixes: tuple[tuple[int, str], ...]
    #: Edit distance of the minimal fix, in codons.
    minimal_distance: int
    description: str

    @property
    def fix_count(self) -> int:
        return len(self.minimal_fixes)


def _motif_set(set_id: str, motif: str, scope: str = "cds") -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "constraint_set_id": set_id,
        "constraints": [
            {
                "id": f"no_{motif.lower()}",
                "class": "blocking",
                "scope": scope,
                "rule": {
                    "kind": "forbidden_motif",
                    "motif": motif,
                    "strand": "sense",
                    "max_occurrences": 0,
                },
            }
        ],
    }


#: The scaffold every SYN-04 case is planted in.
#:
#: `TGG` is tryptophan, which has one codon, so a motif overlapping it can only
#: be broken elsewhere — which is how a *unique* fix is constructible at all.
SYN04_SCAFFOLD = "ATGTGGAAATTTGGATCCGTTTAA"


@functools.cache
def syn04_unique_fix() -> PlantedViolation:
    """One violation with exactly one minimal single-codon fix.

    `TGGAAA` spans the fixed tryptophan and the lysine after it, so only the
    lysine can move and only one of its synonyms escapes. The fix set is
    computed against the real checker rather than asserted — a fixture stating
    its own ground truth would only be testing that two hand-written facts
    agree — and the count is then checked against the construction.
    """
    constraint_set = _motif_set("cst_syn04_unique", "TGGAAA")
    fixes = _single_codon_fixes(SYN04_SCAFFOLD, constraint_set)
    assert len(fixes) == 1, f"expected a unique fix, found {len(fixes)}: {fixes}"
    return PlantedViolation(
        cds=SYN04_SCAFFOLD,
        constraint_set=constraint_set,
        minimal_fixes=fixes,
        minimal_distance=1,
        description="motif overlapping the single-codon tryptophan: one escape",
    )


@functools.cache
def syn04_multiple_equal_fixes() -> PlantedViolation:
    """One violation with several equally-minimal fixes.

    `DEP-SYN-001` and `TST-CONS-002`'s diversity property. The multi-fix
    construction is what makes randomized tie-breaking *testable* rather than
    merely asserted: with one fix available, a repair that always returns the
    same answer is indistinguishable from one that chooses among equals.
    """
    constraint_set = _motif_set("cst_syn04_multi", "TTTGGA")
    fixes = _single_codon_fixes(SYN04_SCAFFOLD, constraint_set)
    assert len(fixes) >= 3, f"expected several equal fixes, found {len(fixes)}"
    return PlantedViolation(
        cds=SYN04_SCAFFOLD,
        constraint_set=constraint_set,
        minimal_fixes=fixes,
        minimal_distance=1,
        description="motif with several synonymous escapes, all equally minimal",
    )


@functools.cache
def syn04_no_single_codon_fix() -> PlantedViolation:
    """A violation no single substitution can remove.

    The distance-2 case, and the reason `minimal_distance` is a field rather
    than a constant. `TGGA` spans tryptophan and the codon after it in a
    position where every synonym of the mutable codon re-forms the motif, so
    repair must either take two edits or give up — and "gave up" is a different
    outcome from "fixed it", which a suite of one-edit cases cannot distinguish.
    """
    constraint_set = _motif_set("cst_syn04_none", "TGGA")
    fixes = _single_codon_fixes(SYN04_SCAFFOLD, constraint_set)
    assert fixes == (), f"expected no single-codon fix, found {fixes}"
    return PlantedViolation(
        cds=SYN04_SCAFFOLD,
        constraint_set=constraint_set,
        minimal_fixes=(),
        minimal_distance=2,
        description="motif no single synonymous substitution can remove",
    )


def _single_codon_fixes(cds: str, constraint_set: dict[str, Any]) -> tuple[tuple[int, str], ...]:
    """Every single-codon substitution that makes the sequence feasible.

    Computed by enumeration against the real checker, not asserted. A fixture
    that stated its own ground truth would be testing that two hand-written
    facts agree.
    """
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    found: list[tuple[int, str]] = []
    for index, codon in enumerate(codons):
        for synonym in _core.synonyms(codon, 1):
            if synonym == codon:
                continue
            candidate = codons.copy()
            candidate[index] = synonym
            if _is_feasible("".join(candidate), constraint_set):
                found.append((index, synonym))
    return tuple(sorted(found))


def _is_feasible(cds: str, constraint_set: dict[str, Any]) -> bool:
    import json

    report = json.loads(
        _core.constraint_violations(json.dumps(constraint_set), cds, [(0, len(cds))], "assembled")
    )
    return not any(v["actionable"] for v in report["violations"])


# --------------------------------------------------------------------------- #
# SYN-05 — infeasible constraint sets
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class InfeasibleCase:
    """An input that must be refused, with which check should fire."""

    cds: str
    constraint_set: dict[str, Any]
    five_prime_flank: str
    three_prime_flank: str
    #: The `CF-*` code the analyser is expected to report.
    expected_code: str
    description: str


@functools.cache
def syn05_cases() -> list[InfeasibleCase]:
    """One purpose-built input per feasibility check.

    Phase 1 ships the two whose checks exist: a blocking violation confined to
    immutable backbone, and a position with no admissible synonym. The
    pigeonhole k-mer bound, the DP-continuation closure and the unattainable GC
    band are `P2`'s, and are absent rather than stubbed — a case whose expected
    outcome no implemented check can produce would assert nothing.
    """
    return [
        InfeasibleCase(
            cds="ATGAAATTTTAA",
            # Construct scope, deliberately. Under `cds` scope the checker
            # never looks at the flank, so nothing is detected and the case
            # would pass for the wrong reason - it would be testing that an
            # unexamined region has no violations.
            constraint_set=_motif_set("cst_syn05_backbone", "GGATCC", scope="construct"),
            five_prime_flank="GGATCCTTGACAATTAATCATCGGCTCGTATAATGTGTGGA",
            three_prime_flank="",
            expected_code="",
            description=(
                "violation confined to the 5' flank: non-actionable by "
                "DEP-CONS-002, so it must NOT be reported as infeasible"
            ),
        ),
        InfeasibleCase(
            cds="ATGTGGAAATAA",
            constraint_set=_motif_set("cst_syn05_no_synonym", "TGG"),
            five_prime_flank="",
            three_prime_flank="",
            expected_code="CF-FEAS-001",
            description=(
                "tryptophan has exactly one codon, so forbidding TGG leaves the "
                "position with no admissible synonym"
            ),
        ),
    ]


# --------------------------------------------------------------------------- #
# SYN-07 — property-based strategies
# --------------------------------------------------------------------------- #


@dataclass
class Strategies:
    """Deterministic generators standing in for Hypothesis strategies.

    Seeded and reproducible rather than adaptive: a property failure that
    cannot be replayed is a report nobody can act on, and the suite is meant to
    run on every commit. Hypothesis itself is a Phase 2 dependency decision;
    these cover the same three shapes the specification names.
    """

    seed: int = 20260828
    _state: int = field(default=0, init=False)

    def __post_init__(self) -> None:
        self._state = self.seed or 1

    def _next(self) -> int:
        # splitmix64, matching the discipline the search RNG uses: a bijection,
        # so adjacent seeds do not collide.
        self._state = (self._state + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
        z = self._state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
        return z ^ (z >> 31)

    def _choice(self, options):
        return options[self._next() % len(options)]

    def arbitrary_protein(self, length: int = 8) -> str:
        residues = [r for r in REDUCED_TABLE if r not in ("M", "*")]
        return "M" + "".join(self._choice(residues) for _ in range(length - 2)) + "*"

    def arbitrary_cds(self, length: int = 8) -> str:
        protein = self.arbitrary_protein(length)
        return "".join(self._choice(REDUCED_TABLE[r]) for r in protein)

    def arbitrary_substitution(self, cds: str) -> tuple[int, str] | None:
        """A synonymous substitution at a position that has one."""
        codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
        options = [
            (index, synonym)
            for index, codon in enumerate(codons)
            for synonym in _core.synonyms(codon, 1)
            if synonym != codon
        ]
        return self._choice(options) if options else None


def _all_codons() -> list[str]:
    return ["".join(c) for c in itertools.product("TCAG", repeat=3)]


#: Every generator, so a caller can iterate the suite rather than name it.
GENERATORS: dict[str, Callable[[], Any]] = {
    "SYN-01/distinct": syn01_distinct,
    "SYN-01/degenerate": syn01_degenerate,
    "SYN-04/unique-fix": syn04_unique_fix,
    "SYN-04/multiple-equal-fixes": syn04_multiple_equal_fixes,
    "SYN-04/no-single-codon-fix": syn04_no_single_codon_fix,
    "SYN-05": syn05_cases,
}
