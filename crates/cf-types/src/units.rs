//! Units, objective direction, and the objective vector.
//!
//! Ticket P1-08. Contract: `planning/02-schemas/common.schema.json` (`Unit`,
//! `ObjectiveDirection`) and `landscape_response.schema.json`
//! (`objective_vector`).
//!
//! Two rules from the planning package are enforced by these types rather than
//! by convention:
//!
//! 1. **Every numeric term carries a unit.** A unitless objective vector is not
//!    comparable across model versions, so `TermValue` has no constructor that
//!    omits one.
//! 2. **`UNKNOWN` is a first-class value and is never imputed.** A provider that
//!    missed its latency budget yields `None`, which is semantically distinct
//!    from any number. Imputing a mean or a worst case would bias the front in
//!    a way that is invisible in the output, so the type makes the distinction
//!    unrepresentable-away: the vector holds `Option<f64>`.

use crate::ids::TermId;

/// Direction of improvement for a term.
///
/// Declared per term in every response so that no module hard-codes a sign
/// convention. All internal optimization is minimization; `Maximize` terms are
/// negated exactly once, at ingress to the metrics layer.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum Direction {
    Minimize,
    Maximize,
}

impl Direction {
    pub fn as_str(self) -> &'static str {
        match self {
            Direction::Minimize => "minimize",
            Direction::Maximize => "maximize",
        }
    }

    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "minimize" => Some(Direction::Minimize),
            "maximize" => Some(Direction::Maximize),
            _ => None,
        }
    }

    /// Convert a raw term value into minimization sense.
    #[inline]
    pub fn to_minimization(self, value: f64) -> f64 {
        match self {
            Direction::Minimize => value,
            Direction::Maximize => -value,
        }
    }
}

/// A unit drawn from the documented vocabulary, or a free-form extension.
///
/// Kept as a string rather than a closed enum because the vocabulary grows with
/// each new term and a closed enum would make adding a term a breaking change —
/// the same reasoning that keeps `requested_terms` free of a wildcard.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Unit(String);

impl Unit {
    pub const DIMENSIONLESS: &'static str = "dimensionless";
    pub const Z_SCORE: &'static str = "z_score";
    pub const KCAL_PER_MOL: &'static str = "kcal/mol";
    pub const RELATIVE_GROWTH_RATE: &'static str = "relative_growth_rate";
    pub const NATS_PER_CODON: &'static str = "nats/codon";
    pub const CODONS_PER_SECOND: &'static str = "codons/s";
    pub const MOLECULES_PER_CELL: &'static str = "molecules/cell";
    pub const LOG2_FOLD_CHANGE: &'static str = "log2_fold_change";

    /// The vocabulary the schema documents. Extensions are permitted but are
    /// reported by `is_documented` so a reviewer can see them.
    pub const DOCUMENTED: &'static [&'static str] = &[
        Self::KCAL_PER_MOL,
        Self::DIMENSIONLESS,
        Self::Z_SCORE,
        Self::CODONS_PER_SECOND,
        Self::MOLECULES_PER_CELL,
        Self::RELATIVE_GROWTH_RATE,
        Self::NATS_PER_CODON,
        Self::LOG2_FOLD_CHANGE,
    ];

    pub fn new(unit: impl Into<String>) -> Self {
        Self(unit.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn is_documented(&self) -> bool {
        Self::DOCUMENTED.contains(&self.0.as_str())
    }
}

impl std::fmt::Display for Unit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

/// The objective vector handed to the search engine.
///
/// The ordering is carried explicitly. Nothing may infer it from map iteration:
/// two runs whose objective order differed would produce incomparable fronts
/// while looking identical in the output.
#[derive(Clone, Debug, PartialEq)]
pub struct ObjectiveVector {
    order: Vec<TermId>,
    sense: Vec<Direction>,
    values: Vec<Option<f64>>,
}

impl ObjectiveVector {
    pub fn new(
        order: Vec<TermId>,
        sense: Vec<Direction>,
        values: Vec<Option<f64>>,
    ) -> Result<Self, &'static str> {
        if order.len() != sense.len() || order.len() != values.len() {
            return Err("order, sense and values must have the same length");
        }
        if order.is_empty() {
            return Err("an objective vector must have at least one term");
        }
        let mut seen = std::collections::HashSet::new();
        for term in &order {
            if !seen.insert(term.as_str()) {
                return Err("duplicate term in the objective order");
            }
        }
        Ok(Self {
            order,
            sense,
            values,
        })
    }

    pub fn len(&self) -> usize {
        self.order.len()
    }

    pub fn is_empty(&self) -> bool {
        self.order.is_empty()
    }

    pub fn order(&self) -> &[TermId] {
        &self.order
    }

    pub fn sense(&self) -> &[Direction] {
        &self.sense
    }

    pub fn values(&self) -> &[Option<f64>] {
        &self.values
    }

    /// True when any term is `UNKNOWN`.
    ///
    /// Such a candidate is placed at the back of its feasibility class and
    /// re-queued for exact evaluation; it is never compared on the missing
    /// term and never imputed.
    pub fn has_unknown(&self) -> bool {
        self.values.iter().any(Option::is_none)
    }

    pub fn get(&self, term: &TermId) -> Option<Option<f64>> {
        self.order
            .iter()
            .position(|t| t == term)
            .map(|i| self.values[i])
    }

    /// The vector in minimization sense, with `UNKNOWN` preserved.
    pub fn minimization_values(&self) -> Vec<Option<f64>> {
        self.values
            .iter()
            .zip(&self.sense)
            .map(|(v, s)| v.map(|x| s.to_minimization(x)))
            .collect()
    }

    /// Whether two vectors are comparable at all: same order, same sense, and
    /// neither carrying an `UNKNOWN`.
    pub fn comparable_with(&self, other: &Self) -> bool {
        self.order == other.order
            && self.sense == other.sense
            && !self.has_unknown()
            && !other.has_unknown()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn terms(names: &[&str]) -> Vec<TermId> {
        names.iter().map(|n| TermId::new(*n).unwrap()).collect()
    }

    #[test]
    fn unit_direction_round_trip_and_sign() {
        for d in [Direction::Minimize, Direction::Maximize] {
            assert_eq!(Direction::parse(d.as_str()), Some(d));
        }
        assert_eq!(Direction::Minimize.to_minimization(1.5), 1.5);
        assert_eq!(Direction::Maximize.to_minimization(1.5), -1.5);
        assert_eq!(Direction::parse("descending"), None);
    }

    #[test]
    fn unit_objective_vector_shape_is_checked() {
        let order = terms(&["usage.cai", "burden.growth_cost"]);
        assert!(ObjectiveVector::new(
            order.clone(),
            vec![Direction::Maximize, Direction::Minimize],
            vec![Some(0.7), Some(0.03)]
        )
        .is_ok());
        assert!(ObjectiveVector::new(
            order.clone(),
            vec![Direction::Maximize],
            vec![Some(0.7), Some(0.03)]
        )
        .is_err());
        assert!(ObjectiveVector::new(vec![], vec![], vec![]).is_err());
        assert!(ObjectiveVector::new(
            terms(&["usage.cai", "usage.cai"]),
            vec![Direction::Maximize, Direction::Maximize],
            vec![Some(1.0), Some(1.0)]
        )
        .is_err());
    }

    /// UNKNOWN is preserved through every transformation; nothing imputes it.
    #[test]
    fn unit_unknown_is_never_imputed() {
        let v = ObjectiveVector::new(
            terms(&["usage.cai", "lm.naturalness"]),
            vec![Direction::Maximize, Direction::Maximize],
            vec![Some(0.7), None],
        )
        .unwrap();
        assert!(v.has_unknown());
        assert_eq!(v.minimization_values(), vec![Some(-0.7), None]);

        let w = ObjectiveVector::new(
            terms(&["usage.cai", "lm.naturalness"]),
            vec![Direction::Maximize, Direction::Maximize],
            vec![Some(0.8), Some(-3.0)],
        )
        .unwrap();
        assert!(
            !w.comparable_with(&v),
            "a vector with UNKNOWN is not comparable"
        );
        assert!(w.comparable_with(&w));
    }

    #[test]
    fn unit_documented_units() {
        assert!(Unit::new(Unit::Z_SCORE).is_documented());
        assert!(!Unit::new("furlongs_per_fortnight").is_documented());
    }
}
