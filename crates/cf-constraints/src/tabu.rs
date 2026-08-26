//! The tabu ledger — the concrete anti-thrashing mechanism.
//!
//! Ticket P1-18a. Contract: `planning/03-algorithms/06-constraint-repair.md` §4.
//!
//! SPLIT NOTE. `P1-18` originally bundled the ledger with its wiring into the
//! search loop. The ledger itself must land with `P1-17`, because repair's
//! signature takes it; the wiring (`tick` driven by the generation clock, and
//! the leak consulted by the mutation operator) needs `P1-23`, and so does
//! `TST-CONS-005`. `P1-18b` carries those.
//!
//! THE FAILURE THIS PREVENTS. Mutation proposes edit `E`, repair reverts it,
//! the next generation proposes `E` again — forever. Recording the codon repair
//! *reverted to* discourages the mutation operator from re-proposing the edit
//! that caused the reversion.
//!
//! THE 10 % LEAK IS DELIBERATE. A hard tabu can make a genuinely good region
//! unreachable if the constraint that triggered it was only violated in a
//! context that no longer applies. The leak lets the population re-test a
//! tabooed edit occasionally, at a rate low enough not to re-establish the
//! cycle.

use std::collections::HashMap;

use cf_types::Codon;

/// Default probability that a tabooed edit is allowed through anyway.
pub const DEFAULT_LEAK: f64 = 0.10;

/// Per-run, shared across the population, decaying on a generation clock.
#[derive(Clone, Debug)]
pub struct TabuLedger {
    /// `(region, codon index, codon) -> generations remaining`.
    entries: HashMap<(usize, usize, Codon), u32>,
    leak: f64,
}

impl Default for TabuLedger {
    fn default() -> Self {
        Self::new(DEFAULT_LEAK)
    }
}

impl TabuLedger {
    pub fn new(leak: f64) -> Self {
        assert!(
            (0.0..1.0).contains(&leak),
            "leak must be a probability below 1; a leak of 1 disables the ledger"
        );
        Self {
            entries: HashMap::new(),
            leak,
        }
    }

    pub fn leak(&self) -> f64 {
        self.leak
    }

    /// Record that repair reverted `reverted` at this position.
    ///
    /// The codon recorded is the one repair moved *away from*, because that is
    /// what the mutation operator would otherwise re-propose.
    pub fn record(&mut self, region: usize, codon_index: usize, reverted: Codon, tenure: u32) {
        if tenure == 0 {
            return;
        }
        self.entries.insert((region, codon_index, reverted), tenure);
    }

    /// One generation of decay. Entries reaching zero are removed.
    pub fn tick(&mut self) {
        self.entries.retain(|_, remaining| {
            *remaining = remaining.saturating_sub(1);
            *remaining > 0
        });
    }

    pub fn is_tabu(&self, region: usize, codon_index: usize, codon: Codon) -> bool {
        self.entries.contains_key(&(region, codon_index, codon))
    }

    pub fn remaining(&self, region: usize, codon_index: usize, codon: Codon) -> u32 {
        *self
            .entries
            .get(&(region, codon_index, codon))
            .unwrap_or(&0)
    }

    /// Whether to suppress a proposed edit.
    ///
    /// `roll` is supplied by the caller rather than drawn here, so the decision
    /// is deterministic given the run's RNG plan and reproducible in a test.
    pub fn should_suppress(
        &self,
        region: usize,
        codon_index: usize,
        codon: Codon,
        roll: f64,
    ) -> bool {
        self.is_tabu(region, codon_index, codon) && roll >= self.leak
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn clear(&mut self) {
        self.entries.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn codon(text: &str) -> Codon {
        Codon::parse(text).unwrap()
    }

    #[test]
    fn unit_record_and_decay() {
        let mut ledger = TabuLedger::new(DEFAULT_LEAK);
        ledger.record(0, 4, codon("CTG"), 3);
        assert!(ledger.is_tabu(0, 4, codon("CTG")));
        assert_eq!(ledger.remaining(0, 4, codon("CTG")), 3);
        // Position and region are part of the key.
        assert!(!ledger.is_tabu(0, 5, codon("CTG")));
        assert!(!ledger.is_tabu(1, 4, codon("CTG")));

        for expected in [2u32, 1] {
            ledger.tick();
            assert_eq!(ledger.remaining(0, 4, codon("CTG")), expected);
        }
        ledger.tick();
        assert!(!ledger.is_tabu(0, 4, codon("CTG")));
        assert!(
            ledger.is_empty(),
            "expired entries must be removed, not kept at zero"
        );
    }

    #[test]
    fn unit_zero_tenure_records_nothing() {
        let mut ledger = TabuLedger::default();
        ledger.record(0, 1, codon("CTG"), 0);
        assert!(ledger.is_empty());
    }

    /// The leak is a real probability, not a rounding artefact.
    #[test]
    fn unit_leak_lets_roughly_one_in_ten_through() {
        let mut ledger = TabuLedger::new(0.10);
        ledger.record(0, 0, codon("CTG"), 5);

        // A roll below the leak passes; at or above it, the edit is suppressed.
        assert!(!ledger.should_suppress(0, 0, codon("CTG"), 0.05));
        assert!(ledger.should_suppress(0, 0, codon("CTG"), 0.10));
        assert!(ledger.should_suppress(0, 0, codon("CTG"), 0.99));

        // Over a uniform sweep, close to 10% get through.
        let passed = (0..1000)
            .filter(|i| !ledger.should_suppress(0, 0, codon("CTG"), *i as f64 / 1000.0))
            .count();
        assert_eq!(passed, 100, "the leak must be exactly the stated rate");

        // A codon that was never recorded is never suppressed.
        assert!(!ledger.should_suppress(0, 0, codon("TTG"), 0.99));
    }

    #[test]
    #[should_panic(expected = "leak must be a probability")]
    fn unit_leak_of_one_is_refused() {
        // A leak of 1 would silently disable the ledger, which is worse than
        // an error because the thrashing it prevents is invisible.
        TabuLedger::new(1.0);
    }
}
