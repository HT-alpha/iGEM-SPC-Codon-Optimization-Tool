//! The saturating burden penalty, as a ranked objective.
//!
//! Ticket P2-14. Specification: `planning/03-algorithms/05-burden-penalty-and-capacity.md` §1.
//! `DEP-BURD-002`, `DEP-BURD-003`, `DEP-BURD-005`.
//!
//! The reported value is `g`, bounded in `[0, g_max]`. The search hinge `H` is
//! not this term: it has no unit a user can read, and the specification keeps
//! it off the objective vector. Tiers 1 and 2 see `g` directly. `g` is not
//! flat on the designs this objective exists to separate — a codon assignment
//! piled onto one synonym versus one spread across the host's own mix — so
//! the hinge is not required to give those tiers a direction.
//!
//! ## What the number is
//!
//! At the tier this build can actually fill from a codon-usage table or a
//! gene-copy map, `g` is a demand-profile mismatch in `[0, 1)`, not a
//! predicted growth-rate cost. `calibrated` stays false. The unit is
//! `dimensionless`. A tier-1 panel fit is a different function
//! ([`crate::fit_capacity_tier1`]) and is not implied by a number this term
//! returns.

use std::collections::BTreeMap;

use cf_landscape_api::{
    Capability, Edit, EvalContext, InteractionRadius, PositionDeltas, Provenance, Rebased,
    TermProvider, TermRole,
};
use cf_types::{
    code, CfError, Codon, CodonSequence, Direction, GeneticCode, Result, TermId, Unit, CODON_COUNT,
};
use serde_json::Value;

use crate::primary_decoder;

/// Soft-max order from §1.1. Specified, not fitted here.
pub const AGGREGATION_P: f64 = 4.0;

/// Saturation exponent from §1.1. Specified, not fitted here.
pub const SATURATION_N: f64 = 2.0;

/// `g_max` at tier 3. §2.3 sets it to 1 so the score stays in `[0, 1)`.
const TIER3_G_MAX: f64 = 1.0;

/// Power mean of order `p`.
///
/// §1.1 writes an L_p norm and, in the same paragraph, says `p = 4` sits
/// between a plain mean and a hard max. An L_p norm is at least the max and
/// grows with the number of pools, so it is not between those two. The power
/// mean is the aggregate that sentence describes. `DEP-BURD-005`.
pub fn power_mean(values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mean: f64 = values.iter().map(|value| value.powf(p)).sum::<f64>() / values.len() as f64;
    mean.powf(1.0 / p)
}

/// `g(u) = g_max · u^n / (1 + u^n)`.
pub fn growth_cost(u_eff: f64, n: f64, g_max: f64) -> f64 {
    if !u_eff.is_finite() || u_eff <= 0.0 {
        return 0.0;
    }
    let powered = u_eff.powf(n);
    g_max * powered / (1.0 + powered)
}

/// How supply was obtained. Recorded in provenance so a mismatch score is not
/// read as a growth-rate prediction.
#[derive(Clone, Debug, PartialEq)]
enum Supply {
    /// Synonymous shares proportional to `relative_adaptiveness`.
    UsagePrior { weights: Box<[f64; CODON_COUNT]> },
    /// Shares proportional to `trna.gene_copy_number`.
    GeneCopy { share: BTreeMap<String, f64> },
}

/// `burden.growth_cost`. Minimize. L1.
#[derive(Clone, Debug)]
pub struct GrowthCostTerm {
    term: TermId,
    supply: Supply,
    /// Engine name, so the lockfile can say which prior produced the number.
    engine: &'static str,
}

impl GrowthCostTerm {
    /// Build the term from a host profile.
    ///
    /// Gene-copy counts win when the profile has a non-empty
    /// `trna.gene_copy_number`: that is a supply prior. Otherwise the codon
    /// usage table is the native demand prior (tier 3). A profile with
    /// neither is refused. No copy numbers are invented for hosts that did
    /// not bring them.
    pub fn from_host_profile(profile: &Value) -> Result<Self> {
        if let Some(supply) = gene_copy_supply(profile)? {
            return Ok(Self {
                term: TermId::new("burden.growth_cost")?,
                supply: Supply::GeneCopy { share: supply },
                engine: "burden_tier3_gene_copy",
            });
        }
        let weights = usage_weights(profile)?;
        Ok(Self {
            term: TermId::new("burden.growth_cost")?,
            supply: Supply::UsagePrior {
                weights: Box::new(weights),
            },
            engine: "burden_tier3_usage_prior",
        })
    }

    fn score_sequence(&self, sequence: &CodonSequence) -> Result<f64> {
        self.score_counts(&Self::counts(sequence), sequence.table())
    }

    fn score_counts(&self, counts: &[u64; CODON_COUNT], code: &GeneticCode) -> Result<f64> {
        let utilizations = match &self.supply {
            Supply::UsagePrior { weights } => usage_utilizations(counts, weights.as_ref(), code)?,
            Supply::GeneCopy { share } => gene_copy_utilizations(counts, share, code)?,
        };
        Ok(growth_cost(
            power_mean(&utilizations, AGGREGATION_P),
            SATURATION_N,
            TIER3_G_MAX,
        ))
    }

    fn counts(sequence: &CodonSequence) -> [u64; CODON_COUNT] {
        let mut counts = [0u64; CODON_COUNT];
        for codon in sequence.codons() {
            if sequence.table().is_stop(*codon) {
                continue;
            }
            counts[codon.index_usize()] += 1;
        }
        counts
    }
}

fn gene_copy_supply(profile: &Value) -> Result<Option<BTreeMap<String, f64>>> {
    let Some(table) = profile
        .pointer("/trna/gene_copy_number")
        .and_then(Value::as_object)
    else {
        return Ok(None);
    };
    if table.is_empty() {
        return Ok(None);
    }
    let mut copies: BTreeMap<String, f64> = BTreeMap::new();
    for (anticodon, value) in table {
        let count = value.as_f64().ok_or_else(|| {
            CfError::new(
                code::BURD_NO_CAPACITY_DATA,
                format!("trna.gene_copy_number.{anticodon} is not a number"),
            )
        })?;
        if count < 0.0 || !count.is_finite() {
            return Err(CfError::new(
                code::BURD_NO_CAPACITY_DATA,
                format!("trna.gene_copy_number.{anticodon} is {count}"),
            ));
        }
        let key = anticodon.to_ascii_uppercase();
        if key.len() != 3
            || !key
                .bytes()
                .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'U'))
        {
            return Err(CfError::new(
                code::BURD_NO_CAPACITY_DATA,
                format!("trna.gene_copy_number key {anticodon} is not an anticodon"),
            ));
        }
        copies.insert(key.replace('U', "T"), count);
    }
    let total: f64 = copies.values().sum();
    if total <= 0.0 {
        return Err(CfError::new(
            code::BURD_NO_CAPACITY_DATA,
            "trna.gene_copy_number sums to zero, so it cannot be a supply prior",
        ));
    }
    Ok(Some(
        copies
            .into_iter()
            .map(|(anticodon, count)| (anticodon, count / total))
            .collect(),
    ))
}

fn usage_weights(profile: &Value) -> Result<[f64; CODON_COUNT]> {
    let table = profile
        .pointer("/codon_usage/table")
        .and_then(Value::as_object)
        .ok_or_else(|| {
            CfError::new(
                code::BURD_NO_CAPACITY_DATA,
                "burden.growth_cost needs trna.gene_copy_number or codon_usage.table; \
                 this host profile has neither, and a default tRNA pool would be a \
                 fabricated host",
            )
        })?;
    let mut weights = [0.0f64; CODON_COUNT];
    for (codon, entry) in table {
        let codon = Codon::parse(codon)?;
        let weight = entry
            .get("relative_adaptiveness")
            .and_then(Value::as_f64)
            .ok_or_else(|| {
                CfError::new(
                    code::HOST_PROFILE_INVALID,
                    format!("codon_usage.table.{codon} has no relative_adaptiveness"),
                )
            })?;
        if !(0.0..=1.0).contains(&weight) {
            return Err(CfError::new(
                code::HOST_PROFILE_INVALID,
                format!("relative_adaptiveness for {codon} is {weight}, outside [0, 1]"),
            ));
        }
        weights[codon.index_usize()] = weight;
    }
    Ok(weights)
}

/// Per-codon utilizations `q / p` inside each multi-codon family.
///
/// `p` is the host's synonymous share, proportional to relative adaptiveness.
/// `q` is the construct's synonymous share. Amino-acid composition is fixed
/// by the protein, so the mismatch that synonymous choice can move is inside
/// the family (`DEP-BURD-001`). A codon the host never uses (`p = 0`) and the
/// construct does is charged at `q / min p` over the family's positive
/// shares. That floor is read from the host table.
fn usage_utilizations(
    counts: &[u64; CODON_COUNT],
    weights: &[f64; CODON_COUNT],
    code: &GeneticCode,
) -> Result<Vec<f64>> {
    let mut utilizations = Vec::new();
    let mut seen = [false; 256];
    for codon in Codon::all() {
        if counts[codon.index_usize()] == 0 || code.is_stop(codon) {
            continue;
        }
        let amino = code.amino_acid(codon) as usize;
        if seen[amino] {
            continue;
        }
        seen[amino] = true;
        let family = code.synonyms(codon);
        if family.len() <= 1 {
            continue;
        }
        let weight_sum: f64 = family.iter().map(|item| weights[item.index_usize()]).sum();
        let count_sum: u64 = family.iter().map(|item| counts[item.index_usize()]).sum();
        if count_sum == 0 {
            continue;
        }
        if weight_sum <= 0.0 {
            return Err(CfError::new(
                code::BURD_NO_CAPACITY_DATA,
                format!(
                    "amino acid {} has demand in this sequence and no positive \
                     relative_adaptiveness in the host table, so there is no native \
                     mix to compare against",
                    code.amino_acid(codon) as char
                ),
            ));
        }
        let mut positive = Vec::new();
        for item in family {
            let share = weights[item.index_usize()] / weight_sum;
            if share > 0.0 {
                positive.push(share);
            }
        }
        let floor = positive.iter().copied().fold(f64::INFINITY, f64::min);
        for item in family {
            let share = weights[item.index_usize()] / weight_sum;
            let observed = counts[item.index_usize()] as f64 / count_sum as f64;
            if share > 0.0 {
                utilizations.push(observed / share);
            } else if observed > 0.0 {
                utilizations.push(observed / floor);
            }
        }
    }
    Ok(utilizations)
}

fn gene_copy_utilizations(
    counts: &[u64; CODON_COUNT],
    share: &BTreeMap<String, f64>,
    code: &GeneticCode,
) -> Result<Vec<f64>> {
    let mut demand: BTreeMap<String, f64> = BTreeMap::new();
    let mut total = 0.0;
    for codon in Codon::all() {
        let count = counts[codon.index_usize()];
        if count == 0 || code.is_stop(codon) {
            continue;
        }
        let Some(decoder) = primary_decoder(codon) else {
            continue;
        };
        *demand.entry(decoder.as_str()).or_insert(0.0) += count as f64;
        total += count as f64;
    }
    if total <= 0.0 {
        return Ok(Vec::new());
    }
    let floor = share
        .values()
        .copied()
        .filter(|value| *value > 0.0)
        .fold(f64::INFINITY, f64::min);
    if !floor.is_finite() {
        return Err(CfError::new(
            code::BURD_NO_CAPACITY_DATA,
            "trna.gene_copy_number has no positive pool",
        ));
    }
    let mut pools: Vec<String> = share.keys().cloned().collect();
    for key in demand.keys() {
        if !share.contains_key(key) {
            pools.push(key.clone());
        }
    }
    pools.sort();
    pools.dedup();
    let mut utilizations = Vec::new();
    for pool in pools {
        let supply = share.get(&pool).copied().unwrap_or(0.0);
        let observed = demand.get(&pool).copied().unwrap_or(0.0) / total;
        if supply > 0.0 {
            utilizations.push(observed / supply);
        } else if observed > 0.0 {
            utilizations.push(observed / floor);
        }
    }
    Ok(utilizations)
}

/// Counts pinned by `rebase`, so a delta does not rescan the sequence.
#[derive(Clone, Copy, Debug)]
struct BurdenState {
    counts: [u64; CODON_COUNT],
}

impl TermProvider for GrowthCostTerm {
    fn term_id(&self) -> &TermId {
        &self.term
    }

    fn unit(&self) -> Unit {
        // Tier 3. Not `relative_growth_rate`: that unit would claim a
        // prediction this prior does not make (`DEP-BURD-003`).
        Unit::new(Unit::DIMENSIONLESS)
    }

    fn direction(&self) -> Direction {
        Direction::Minimize
    }

    fn capability(&self) -> Capability {
        Capability::L1
    }

    fn max_role(&self) -> TermRole {
        TermRole::ParetoAxis
    }

    /// One codon changes only its own count. The penalty is global, and the
    /// delta is exact from those counts, so the interaction radius is one
    /// position.
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }

    fn provenance(&self) -> Provenance {
        Provenance::builtin(self.engine, env!("CARGO_PKG_VERSION"))
    }

    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        self.score_sequence(ctx.sequence)
    }

    fn rebase(&self, ctx: &EvalContext<'_>) -> Result<Rebased> {
        let state = BurdenState {
            counts: GrowthCostTerm::counts(ctx.sequence),
        };
        Ok(Rebased::with_state(
            self.score_sequence(ctx.sequence)?,
            state,
        ))
    }

    fn delta(&self, ctx: &EvalContext<'_>, base: &Rebased, edits: &[Edit]) -> Result<f64> {
        let code = ctx.sequence.table();
        let mut counts = base
            .state::<BurdenState>()
            .map(|state| state.counts)
            .unwrap_or_else(|| GrowthCostTerm::counts(ctx.sequence));
        let before = self.score_counts(&counts, code)?;
        for edit in edits {
            if edit.codon_index >= ctx.sequence.codon_count() {
                return Err(CfError::new(
                    code::SEQ_INDEX_OUT_OF_RANGE,
                    format!("codon {} is outside the sequence", edit.codon_index),
                ));
            }
            let old = ctx.sequence.codons()[edit.codon_index];
            if code.amino_acid(old) != code.amino_acid(edit.codon) {
                return Err(CfError::new(
                    code::BURD_NON_SYNONYMOUS,
                    format!("{old} and {} do not encode the same amino acid", edit.codon),
                ));
            }
            if !code.is_stop(old) {
                let slot = &mut counts[old.index_usize()];
                *slot = slot.saturating_sub(1);
            }
            if !code.is_stop(edit.codon) {
                counts[edit.codon.index_usize()] += 1;
            }
        }
        let after = self.score_counts(&counts, code)?;
        Ok(after - before)
    }

    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        positions
            .iter()
            .map(|&index| {
                if index >= ctx.sequence.codon_count() {
                    return Err(CfError::new(
                        code::SEQ_INDEX_OUT_OF_RANGE,
                        format!("codon {index} is outside the sequence"),
                    ));
                }
                let current = ctx.sequence.codons()[index];
                let alternatives = ctx
                    .sequence
                    .synonyms_at(index)
                    .iter()
                    .copied()
                    .filter(|codon| *codon != current)
                    .map(|codon| {
                        self.delta(ctx, base, &[Edit::new(index, codon)])
                            .map(|change| (codon, change))
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(PositionDeltas {
                    codon_index: index,
                    alternatives,
                })
            })
            .collect()
    }
}
