//! Capacity estimators — the three tiers of `θ`.
//!
//! Ticket P2-13. Spec: `planning/03-algorithms/05-burden-penalty-and-capacity.md` §2.
//! `ADR-0004`, `DEP-BURD-003`.
//!
//! `θ` is the threshold of the saturating burden penalty, and how it is
//! estimated depends entirely on what data the host has. The tier is selected
//! from `HostProfile.data_completeness`, recorded in the profile, propagated
//! into term provenance, and surfaced to the user — because a tier-3 burden
//! number means something categorically different from a tier-1 one.
//!
//! ## The honesty requirement is load-bearing, not decorative
//!
//! Tier 3 is where most hosts land. Its output is a **demand-profile mismatch
//! score**, not a predicted growth-rate cost, and `DEP-BURD-003` says reporting
//! the same unit at all three tiers would be the single most misleading thing
//! DNA Design Suite could do. So the unit, the calibration flag and the rationale tag
//! are part of the return type rather than something a caller attaches
//! afterwards and might forget to: [`CapacityEstimate`] cannot be constructed
//! without them, and tier 3 cannot be constructed claiming to be calibrated.
//!
//! ## What is not fabricated here
//!
//! Tier 3 assumes a host's tRNA pools are on average matched to its own native
//! codon demand. That is the same empirical regularity underlying tAI, so it is
//! a documented approximation rather than an invention — but it is only sound
//! if the native demand actually comes from the host's genome. An estimator
//! that returned a plausible constant would satisfy every semantic check in
//! `TST-BURD-005` while being entirely disconnected from the host, which is why
//! the tests for this module assert that the output *moves with* the supplied
//! `trna.gene_copy_number` and `tai.s_values`.

use std::collections::BTreeMap;

use cf_types::{code, CfError, Result};

/// What the burden number means at a given tier.
///
/// Not a formatting detail. At tier 3 `u = D/S` is a mismatch score in `[0,1)`;
/// at tiers 1 and 2 it is a utilization ratio against an absolute capacity.
/// Presenting the second as the first is the failure `DEP-BURD-003` names.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThetaUnit {
    /// Tier 3. Dimensionless mismatch against the host's own native profile.
    FractionOfNativeDemand,
    /// Tier 1. Fit to measured burden/expression curves.
    AaPerCellPerSecond,
    /// Tier 2. Ribosome count x elongation rate.
    RibosomeOccupancyFraction,
}

impl ThetaUnit {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::FractionOfNativeDemand => "fraction_of_native_demand",
            Self::AaPerCellPerSecond => "aa_per_cell_per_s",
            Self::RibosomeOccupancyFraction => "ribosome_occupancy_fraction",
        }
    }

    /// The unit the API reports for `burden.growth_cost` at this tier.
    ///
    /// Tier 3 reports `dimensionless` because the number is not a growth-rate
    /// prediction and must not be dressed as one.
    pub fn reported_unit(&self) -> &'static str {
        match self {
            Self::FractionOfNativeDemand => "dimensionless",
            _ => "relative_growth_rate",
        }
    }
}

/// A capacity estimate, carrying the semantics it must be reported with.
///
/// The honesty fields are not optional and not settable after the fact. A
/// tier-3 estimate that could be constructed and *then* marked calibrated would
/// erode exactly the way `ADR-0004` warns about — one caller at a time, each
/// individually reasonable.
#[derive(Clone, Debug, PartialEq)]
pub struct CapacityEstimate {
    tier: u8,
    theta: f64,
    unit: ThetaUnit,
    /// Per-anticodon supply shares, summing to 1 at tier 3.
    supply: BTreeMap<[u8; 3], f64>,
    g_max: f64,
    calibrated: bool,
    rationale_tags: Vec<&'static str>,
}

impl CapacityEstimate {
    pub fn tier(&self) -> u8 {
        self.tier
    }

    pub fn theta(&self) -> f64 {
        self.theta
    }

    pub fn unit(&self) -> ThetaUnit {
        self.unit
    }

    pub fn supply(&self) -> &BTreeMap<[u8; 3], f64> {
        &self.supply
    }

    pub fn g_max(&self) -> f64 {
        self.g_max
    }

    /// Whether the estimate rests on measurements in this host.
    ///
    /// False at tier 3 by construction; there is no setter.
    pub fn calibrated(&self) -> bool {
        self.calibrated
    }

    pub fn rationale_tags(&self) -> &[&'static str] {
        &self.rationale_tags
    }

    /// The unit `burden.growth_cost` is reported in.
    pub fn reported_unit(&self) -> &'static str {
        self.unit.reported_unit()
    }
}

/// What the host has, as read from `data_completeness`.
#[derive(Clone, Copy, Debug, Default)]
pub struct DataCompleteness {
    pub has_burden_panel: bool,
    pub has_ribosome_count: bool,
    pub has_trna_gene_copies: bool,
    pub has_codon_usage: bool,
}

impl DataCompleteness {
    /// Select the tier from what the host actually has.
    ///
    /// Highest tier the data supports, never higher. A host with a burden panel
    /// gets tier 1; without it but with ribosome counts, tier 2; otherwise
    /// tier 3, which is where most hosts land and which is the case
    /// `ADR-0004` says getting right is the difference between DNA Design Suite being
    /// usable outside model organisms and not.
    pub fn tier(&self) -> u8 {
        if self.has_burden_panel {
            1
        } else if self.has_ribosome_count {
            2
        } else {
            3
        }
    }
}

/// Tier 3 — relative-to-native normalization, the no-data case.
///
/// ```text
/// weight_g <- CAI(gene, highly_expressed_reference)^rho     # rho = 2
/// D_native_a <- sum_g weight_g . sum_x n_g(x) . w(x, a)
/// S_a <- D_native_a / sum_b D_native_b
/// theta_a <- 1
/// ```
///
/// `native_demand` is the host's own expression-weighted codon demand per
/// anticodon, computed from its annotated genome. Passing it in rather than
/// deriving it here keeps `cf-burden` free of a genome-parsing dependency, and
/// it is the value the tests vary to check the estimate actually tracks the
/// host rather than returning a plausible constant.
pub fn fit_capacity_tier3(native_demand: &BTreeMap<[u8; 3], f64>) -> Result<CapacityEstimate> {
    if native_demand.is_empty() {
        return Err(CfError::new(
            code::BURD_NO_CAPACITY_DATA,
            "tier 3 needs the host's native codon demand; with no genome-derived \
             demand there is nothing to normalize against, and a default supply \
             profile would be a fabricated host",
        ));
    }
    let total: f64 = native_demand.values().sum();
    // Written out rather than as `!(total > 0.0)`: NaN must be refused too, and
    // the negated form hides that it is doing two jobs.
    if !total.is_finite() || total <= 0.0 {
        return Err(CfError::new(
            code::BURD_NO_CAPACITY_DATA,
            "the host's native demand sums to zero; it cannot be normalized",
        ));
    }

    let supply: BTreeMap<[u8; 3], f64> = native_demand
        .iter()
        .map(|(anticodon, demand)| (*anticodon, demand / total))
        .collect();

    Ok(CapacityEstimate {
        tier: 3,
        // theta = 1 by construction: at this tier u = D/S is a ratio against
        // the host's own profile, so the threshold is "matches native".
        theta: 1.0,
        unit: ThetaUnit::FractionOfNativeDemand,
        supply,
        // g_max = 1 so growth_cost lands in [0,1) as a mismatch score rather
        // than being scaled by an asymptotic growth cost this tier cannot know.
        g_max: 1.0,
        calibrated: false,
        rationale_tags: vec!["low_confidence_host"],
    })
}

/// Tier 2 — ribosome count x elongation rate.
pub fn fit_capacity_tier2(
    ribosomes_per_cell: f64,
    elongation_rate_aa_per_s: f64,
    native_demand: &BTreeMap<[u8; 3], f64>,
) -> Result<CapacityEstimate> {
    if !ribosomes_per_cell.is_finite()
        || ribosomes_per_cell <= 0.0
        || !elongation_rate_aa_per_s.is_finite()
        || elongation_rate_aa_per_s <= 0.0
    {
        return Err(CfError::new(
            code::BURD_NO_CAPACITY_DATA,
            "tier 2 needs a positive ribosome count and elongation rate",
        ));
    }
    let mut estimate = fit_capacity_tier3(native_demand)?;
    estimate.tier = 2;
    estimate.theta = ribosomes_per_cell * elongation_rate_aa_per_s;
    estimate.unit = ThetaUnit::RibosomeOccupancyFraction;
    estimate.g_max = 0.6;
    // Still not calibrated *in this host against measured burden* — tier 2
    // derives capacity from physiology, not from a burden/expression curve.
    // Calling it calibrated would claim a measurement nobody took.
    estimate.calibrated = false;
    estimate.rationale_tags = vec![];
    Ok(estimate)
}

/// Tier 1 — fit to measured burden/expression data.
///
/// Requires a panel of >= 20 constructs spanning >= 1 order of magnitude in
/// expression, measured in the same host with both an expression and a
/// growth-rate readout. The threshold is a real requirement rather than a
/// suggestion: a curve fit to 5 points over half a decade of expression will
/// produce a `θ` with a confident-looking value and no information in it.
pub fn fit_capacity_tier1(
    panel: &[(f64, f64)],
    native_demand: &BTreeMap<[u8; 3], f64>,
) -> Result<CapacityEstimate> {
    if panel.len() < 20 {
        return Err(CfError::new(
            code::BURD_NO_CAPACITY_DATA,
            format!(
                "tier 1 needs at least 20 measured constructs, got {}; a curve fit to \
                 fewer produces a confident-looking theta with no information in it",
                panel.len()
            ),
        ));
    }
    let expressions: Vec<f64> = panel.iter().map(|(e, _)| *e).filter(|e| *e > 0.0).collect();
    if expressions.len() != panel.len() {
        return Err(CfError::new(
            code::BURD_NO_CAPACITY_DATA,
            "tier 1 panel contains a non-positive expression readout",
        ));
    }
    let low = expressions.iter().cloned().fold(f64::INFINITY, f64::min);
    let high = expressions
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    if high / low < 10.0 {
        return Err(CfError::new(
            code::BURD_NO_CAPACITY_DATA,
            "tier 1 panel spans less than one order of magnitude in expression; the \
             saturating region is not observed, so theta is unidentifiable",
        ));
    }

    let mut estimate = fit_capacity_tier3(native_demand)?;
    estimate.tier = 1;
    // theta at the half-maximal growth cost, located on the measured curve.
    estimate.theta = half_maximal_expression(panel);
    estimate.unit = ThetaUnit::AaPerCellPerSecond;
    estimate.g_max = panel
        .iter()
        .map(|(_, cost)| *cost)
        .fold(f64::NEG_INFINITY, f64::max)
        .clamp(0.0, 1.0);
    estimate.calibrated = true;
    estimate.rationale_tags = vec![];
    Ok(estimate)
}

/// Expression at which the measured growth cost reaches half its maximum.
fn half_maximal_expression(panel: &[(f64, f64)]) -> f64 {
    let max_cost = panel
        .iter()
        .map(|(_, cost)| *cost)
        .fold(f64::NEG_INFINITY, f64::max);
    let target = max_cost / 2.0;
    let mut sorted: Vec<(f64, f64)> = panel.to_vec();
    sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    for window in sorted.windows(2) {
        let (x0, y0) = window[0];
        let (x1, y1) = window[1];
        if (y0 - target) * (y1 - target) <= 0.0 && (y1 - y0).abs() > f64::EPSILON {
            return x0 + (target - y0) * (x1 - x0) / (y1 - y0);
        }
    }
    sorted.last().map(|(x, _)| *x).unwrap_or(f64::NAN)
}
