"""Interface stubs for the CodonForge predictor plugin seam.

NOT IMPLEMENTATION CODE. This is a typed interface listing, provided so that the
landscape contract can be reviewed as types rather than only as JSON Schema.
Every structure here has a normative JSON Schema counterpart in
`planning/02-schemas/`; where the two disagree, the JSON Schema wins.

The Rust side of the same contract (`cf-landscape-api`) mirrors this file trait
for trait; the `contract_*` tests in each module assert that both sides accept
the same golden documents.

THE CAPABILITY MODEL IS THE POINT OF THIS FILE.  A predictor declares what it
can do, and the search engine plans around it:

    L0  Scalar   score a full sequence                       -- every predictor
    L1  Delta    exact change for a single codon substitution -- decomposable ones
    L2  Gradient d(term)/d(relaxed codon logits)              -- differentiable ones

Only L0 is mandatory. Requiring L2 universally would force every predictor into
one autograd runtime and destroy the swappable-layer goal; see ADR-0002.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Literal, Mapping, Protocol, Sequence, runtime_checkable

# --------------------------------------------------------------------------- #
# value objects (mirrors of common.schema.json)
# --------------------------------------------------------------------------- #

TermId = str          # namespaced, e.g. "epistasis.dwell_deviation"
HostId = str          # e.g. "escherichia_coli:562/K-12_MG1655"
ArtifactId = str      # e.g. "cf:epistasis_potts:escherichia_coli:562:1.4.0"
Digest = str          # "sha256:<64 hex>"
CodonIndex = int


class Capability(Enum):
    L0_SCALAR = 0
    L1_DELTA = 1
    L2_GRADIENT = 2


class Direction(Enum):
    MINIMIZE = "minimize"
    MAXIMIZE = "maximize"


class BudgetTier(Enum):
    FAST = "fast"
    STANDARD = "standard"
    EXACT = "exact"


@dataclass(frozen=True)
class Uncertainty:
    method: str
    calibrated: bool
    sigma: float | None = None
    interval: tuple[float, float] | None = None
    coverage: float = 0.9
    n_members: int | None = None


@dataclass(frozen=True)
class Provenance:
    resolved_by: Literal["pin", "tenant_channel", "host_channel",
                         "nearest_host", "base", "builtin"]
    artifact_id: ArtifactId | None = None
    digest: Digest | None = None
    engine: str | None = None
    engine_version: str | None = None
    fallback: bool = False


@dataclass(frozen=True)
class TermValue:
    value: float
    unit: str
    direction: Direction
    provenance: Provenance
    reference_range: tuple[float, float] | None = None
    uncertainty: Uncertainty | None = None
    capability: Capability = Capability.L0_SCALAR


# --------------------------------------------------------------------------- #
# the sequence type (mirror of sequence.schema.json)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CodonSequence:
    """Frame safety is structural: `codons` holds codon INDICES (0..63), and the
    only mutating operation exposed anywhere is a synonymous substitution at a
    codon index. Untranslated context lives in the flank fields and is never
    edited, so no operation can shift the reading frame.
    """

    codons: Sequence[int]
    translation_table: int
    protein: str
    five_prime_flank: str
    three_prime_flank: str
    locked: frozenset[CodonIndex]
    digest: Digest

    def substitute(self, index: CodonIndex, codon: int) -> "CodonSequence":
        """Returns a new sequence. Raises if `codon` is not synonymous with the
        codon currently at `index`, or if `index` is locked. Enforcing the
        synonymy precondition here, rather than by convention, is what makes
        protein invariance a type-level property instead of a test-level one.
        """
        ...


@dataclass(frozen=True)
class HostContext:
    host_id: HostId
    host_profile_digest: Digest
    condition_id: str | None


@dataclass(frozen=True)
class ConstructRef:
    construct_digest: Digest
    cistron_id: str | None


# --------------------------------------------------------------------------- #
# the delta table (the reason tier 2 is affordable)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class DeltaTable:
    """Exact change in each term for every synonymous alternative at every listed
    position.

    `exact` distinguishes closed-form deltas from finite-difference or surrogate
    estimates. The repair operator requires `exact=True` for its minimality
    property to hold, and asserts it rather than assuming it.
    """

    base_digest: Digest
    exact: bool
    positions: Sequence[CodonIndex]
    # entries[position_slot][alternative_slot][term] -> delta
    entries: Sequence[Sequence[Mapping[TermId, float]]]
    alternatives: Sequence[Sequence[int]]


# --------------------------------------------------------------------------- #
# the predictor protocols
# --------------------------------------------------------------------------- #


@runtime_checkable
class ScalarPredictor(Protocol):
    """L0. The only mandatory interface."""

    @property
    def provides(self) -> Sequence[TermId]: ...

    @property
    def artifact(self) -> Provenance: ...

    def capability(self, term: TermId) -> Capability: ...

    def score(
        self,
        sequence: CodonSequence,
        host: HostContext,
        construct: ConstructRef,
        terms: Sequence[TermId],
        tier: BudgetTier,
        want_uncertainty: bool = False,
    ) -> Mapping[TermId, TermValue]: ...

    def score_batch(
        self,
        sequences: Sequence[CodonSequence],
        host: HostContext,
        construct: ConstructRef,
        terms: Sequence[TermId],
        tier: BudgetTier,
        want_uncertainty: bool = False,
    ) -> Sequence[Mapping[TermId, TermValue]]:
        """Batch entry point. This is the ONLY method the out-of-core worker
        exposes across the process boundary, and it is called at most once per
        generation (ADR-0001). A predictor with no batching advantage may
        implement it as a loop over `score`.
        """
        ...


@runtime_checkable
class DeltaPredictor(ScalarPredictor, Protocol):
    """L1. Implemented by every decomposable term: CAI, GC, k-mer, the Potts
    context model, the tRNA/burden ledger, windowed structure accessibility.
    """

    def delta(
        self,
        base: CodonSequence,
        index: CodonIndex,
        new_codon: int,
        terms: Sequence[TermId],
    ) -> Mapping[TermId, float]:
        """Exact change from a single synonymous substitution. Must satisfy
        delta == score(after) - score(before) to within 1e-9; TST-EPI-002 and
        its siblings enforce this, because a wrong delta silently corrupts every
        tier-2 result while leaving tier-3 results looking plausible.
        """
        ...

    def delta_table(
        self,
        base: CodonSequence,
        positions: Sequence[CodonIndex],
        terms: Sequence[TermId],
    ) -> DeltaTable: ...

    def rebase(self, base: CodonSequence) -> None:
        """Reset the incremental caches to a new base sequence. Called once per
        block-CA restart, not per move.
        """
        ...


@runtime_checkable
class GradientPredictor(ScalarPredictor, Protocol):
    """L2. Rare and optional.

    The design space is discrete, so a gradient only exists with respect to a
    continuous relaxation. That relaxation is fixed by contract: a per-position
    simplex over the synonymous set. Any other reading of "gradient" is rejected
    at review.
    """

    def gradient(
        self,
        base: CodonSequence,
        terms: Sequence[TermId],
        temperature: float = 1.0,
    ) -> Mapping[TermId, Sequence[Sequence[float | None]]]:
        """shape = [codon_count, max_synonyms]; entries beyond a position's
        synonym count are None.
        """
        ...


# --------------------------------------------------------------------------- #
# swap seams
# --------------------------------------------------------------------------- #


@runtime_checkable
class StructureEngine(Protocol):
    """A folding backend. Stateless and pure.

    AMD-013 / DEP-STRU-002. DNA Design Suite ships no folding engine: ViennaRNA and
    LinearFold are user-installed and detected at runtime, never vendored. Since
    this suite redistributes neither, their licences' no-fee-redistribution
    condition is satisfied by construction.

    Three things this protocol deliberately lacks, each removed for a reason:

      * `local_accessibility_incremental` - it took a changed range and no
        sequence, which only means anything if the engine retains state between
        calls. That assumed an in-repo kernel, which AMD-013 removed.
        A stateless C library cannot implement it, and requiring it would put a
        cache inside every future engine instead of once in P2-10.
      * `global_mfe` - AMD-001 keeps global MFE off the per-move path, and the
        cheapest enforcement is to leave the search loop nothing to call. An
        exact-tier re-ranking path needs a separate interface. `DEP-STRU-007`
        is that separate report: ViennaRNA `RNA.fold` when it is installed,
        otherwise a Nussinov Watson-Crick pair count up to 400 nt, labelled
        as a pair count. It is not `DEP-STRU-006` (the opt-in paired-fraction
        objective) and not a method on this protocol.
      * `supports_*` capability flags - for methods that no longer exist.

    Engines are called from the Python side as an L0 term through P2-05's batch
    protocol, so the Rust hot path never crosses the boundary per move.
    """

    def identity(self) -> EngineIdentity: ...

    def local_accessibility(
        self, seq: str, window_nt: int, span_nt: int, temperature_c: float
    ) -> Sequence[float]:
        """One value per window *end* position, not per sequence position.

        A window shorter than its own length cannot end near the 5' edge, and
        engines report exactly 0.0 there. Returning that verbatim would claim
        those positions are maximally inaccessible when the question was never
        asked, so undefined positions have no entry at all.
        """
        ...


@dataclass(frozen=True)
class EngineIdentity:
    """What the lockfile records about the engine that produced the energies.

    DEP-STRU-001's reproducibility half, which survives the AMD-013 licence
    correction unchanged: which engine, version and parameter table were in
    force decides every number a run reported.
    """

    name: str
    version: str
    parameter_set: str


@runtime_checkable
class CapacityEstimator(Protocol):
    """Burden capacity estimation, tiered by available host data.

    tier 1: fit to measured burden/expression curves
    tier 2: derived from ribosome count x elongation rate
    tier 3: relative-to-native demand profile (genome only)

    The tier is not an implementation detail. A tier-3 estimate is a
    demand-profile mismatch score, not a growth-rate prediction, and the unit it
    reports differs accordingly. TST-BURD-005 asserts this.
    """

    tier: Literal[1, 2, 3]

    def estimate(self, host_profile: object) -> tuple[float, str, Uncertainty]: ...


@runtime_checkable
class TermRegistry(Protocol):
    """Runtime plugin registry. This is the cycle-breaker: `cf-landscape`
    depends on this protocol, never on a concrete predictor crate, which is what
    makes the layers independently swappable (see
    planning/01-modules/02-module-dependency-graph.md).
    """

    def register(self, predictor: ScalarPredictor) -> None: ...

    def resolve(self, term: TermId) -> ScalarPredictor: ...

    def capabilities(self) -> Mapping[TermId, Capability]: ...
