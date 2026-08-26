//! Term roles, and the type that keeps a diagnostic term out of a ranking.
//!
//! Ticket P1-20, departure `DEP-LAND-002`.
//!
//! `AMD-002` says folding kinetics ships as a diagnostic and becomes a Pareto
//! axis only if `EXP-06` justifies it. Stated as policy, that survives exactly
//! as long as nobody appends a term id to an objective-vector order by hand.
//! [`ObjectiveAxis`] makes it survive longer: the objective vector is built
//! from `ObjectiveAxis` values, and the only constructor consults the role. A
//! diagnostic term has no path into a ranking that does not go through a
//! compile error or an explicit `CF-LAND-004`.
//!
//! This is the same device as `ActionableViolation` in `cf-constraints`: where
//! a rule has to hold across code that has not been written yet, encode it in a
//! type rather than in a comment.

use cf_types::{code, CfError, Result, TermId};
use serde::{Deserialize, Serialize};

/// What a term is *for* in a given job.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum TermRole {
    /// A dimension of the Pareto front. Participates in dominance.
    #[default]
    ParetoAxis,
    /// Demoted from an axis to a bound — the remedy for many-objective front
    /// explosion. Evaluated and enforced, but not an axis.
    EpsilonConstraint,
    /// Computed and reported, never ranked on.
    ///
    /// The reason this variant exists at all is `AMD-002`: a term can be worth
    /// showing a user long before it is worth optimizing against.
    Diagnostic,
}

impl TermRole {
    pub fn as_str(self) -> &'static str {
        match self {
            TermRole::ParetoAxis => "pareto_axis",
            TermRole::EpsilonConstraint => "epsilon_constraint",
            TermRole::Diagnostic => "diagnostic",
        }
    }

    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "pareto_axis" => Some(TermRole::ParetoAxis),
            "epsilon_constraint" => Some(TermRole::EpsilonConstraint),
            "diagnostic" => Some(TermRole::Diagnostic),
            _ => None,
        }
    }

    /// Whether a term in this role may occupy a slot in the objective vector.
    #[inline]
    pub fn is_axis(self) -> bool {
        matches!(self, TermRole::ParetoAxis)
    }
}

impl std::fmt::Display for TermRole {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A term that has been checked to be a Pareto axis.
///
/// The objective vector is assembled from these and from nothing else, so
/// "a diagnostic term cannot influence a ranking" is a property of the type
/// system rather than a property of reviewer attention.
#[derive(Clone, PartialEq, Eq, Debug, Hash)]
pub struct ObjectiveAxis(TermId);

impl ObjectiveAxis {
    /// Promote a term to an axis, or refuse.
    ///
    /// Refuses `epsilon_constraint` as well as `diagnostic`. An epsilon
    /// constraint is enforced, not ranked; silently treating it as an axis
    /// would change the dimensionality of the front, which is the failure this
    /// type exists to prevent, just from the other direction.
    pub fn new(term: TermId, role: TermRole) -> Result<Self> {
        if !role.is_axis() {
            return Err(CfError::new(
                code::LAND_ROLE_NOT_AXIS,
                format!("term {term} has role {role} and cannot be a Pareto axis"),
            ));
        }
        Ok(Self(term))
    }

    pub fn term(&self) -> &TermId {
        &self.0
    }

    pub fn into_term(self) -> TermId {
        self.0
    }
}

impl std::fmt::Display for ObjectiveAxis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn term(name: &str) -> TermId {
        TermId::new(name).unwrap()
    }

    #[test]
    fn unit_role_round_trips() {
        for r in [
            TermRole::ParetoAxis,
            TermRole::EpsilonConstraint,
            TermRole::Diagnostic,
        ] {
            assert_eq!(TermRole::parse(r.as_str()), Some(r));
        }
        assert_eq!(TermRole::parse("tiebreak"), None);
        assert_eq!(TermRole::default(), TermRole::ParetoAxis);
    }

    /// The AMD-002 gate, at the level of the single type that enforces it.
    #[test]
    fn unit_only_a_pareto_axis_becomes_an_objective_axis() {
        assert!(ObjectiveAxis::new(term("usage.cai"), TermRole::ParetoAxis).is_ok());

        let refused = ObjectiveAxis::new(term("kinetics.pause_concordance"), TermRole::Diagnostic)
            .unwrap_err();
        assert_eq!(refused.code, "CF-LAND-004");

        assert_eq!(
            ObjectiveAxis::new(term("usage.gc_content"), TermRole::EpsilonConstraint)
                .unwrap_err()
                .code,
            "CF-LAND-004"
        );
    }
}
