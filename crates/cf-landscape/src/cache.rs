//! The evaluation cache.
//!
//! Ticket P1-22, `TST-LAND-001`.
//!
//! ## Why the namespace is mandatory
//!
//! The key carries the landscape lockfile digest. Two jobs that resolved
//! different artifact versions for the same term ask the same question and are
//! entitled to different answers, and a cache that could not tell them apart
//! would serve one job's answer to the other — silently, and only when the
//! cache happened to be warm, which is the worst possible reproducibility bug
//! because it does not reproduce.
//!
//! The digest arrives as a *value*, supplied by the orchestrator. `DEP-LAND-001`:
//! `cf-landscape` is a `rust-core` crate and `cf-modelreg` is a Python package,
//! so the dependency could not have gone the other way even if it had been
//! wanted.
//!
//! ## Why a delta query is keyed by `(base, edits)`
//!
//! The obvious key for any query is the digest of the sequence being scored.
//! For a delta query that would mean materializing the mutated sequence and
//! hashing its canonical document — `O(L)` work, and a SHA-256 at that, to
//! serve a lookup that exists so the caller can avoid `O(L)` work. Keying on
//! the base digest plus the edit list is exact, costs the length of the edit
//! list, and is sound: two routes to the same sequence occupy two entries and
//! both are correct.

use std::collections::HashMap;

use cf_landscape_api::Edit;
use cf_types::{Digest, TermId};

/// Identity of one cached term evaluation.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct CacheKey {
    namespace: String,
    base: String,
    edits: Vec<(usize, u8)>,
    term: String,
}

impl CacheKey {
    pub fn new(namespace: &str, base: &Digest, edits: &[Edit], term: &TermId) -> Self {
        let mut edits: Vec<(usize, u8)> = edits
            .iter()
            .map(|e| (e.codon_index, e.codon.index()))
            .collect();
        // The same set of substitutions is the same sequence however the caller
        // happened to order them.
        edits.sort_unstable();
        Self {
            namespace: namespace.to_owned(),
            base: base.as_str().to_owned(),
            edits,
            term: term.as_str().to_owned(),
        }
    }
}

/// A bounded, insertion-ordered evaluation cache.
///
/// Eviction is first-in-first-out rather than least-recently-used. The access
/// pattern is a sweep over a population, so recency carries little information
/// and an LRU's bookkeeping would cost more than it saved.
#[derive(Debug)]
pub struct EvaluationCache {
    namespace: String,
    entries: HashMap<CacheKey, f64>,
    order: std::collections::VecDeque<CacheKey>,
    capacity: usize,
    hits: u64,
    misses: u64,
}

impl EvaluationCache {
    /// Default capacity, in term evaluations.
    ///
    /// Sized to hold a generation of a few hundred candidates across a handful
    /// of terms with room to spare. Not derived from a measurement, which is
    /// why it is a named constant a caller can override rather than a literal
    /// buried in a constructor.
    pub const DEFAULT_CAPACITY: usize = 65_536;

    pub fn new(namespace: impl Into<String>) -> Self {
        Self::with_capacity(namespace, Self::DEFAULT_CAPACITY)
    }

    pub fn with_capacity(namespace: impl Into<String>, capacity: usize) -> Self {
        Self {
            namespace: namespace.into(),
            entries: HashMap::new(),
            order: std::collections::VecDeque::new(),
            capacity: capacity.max(1),
            hits: 0,
            misses: 0,
        }
    }

    pub fn namespace(&self) -> &str {
        &self.namespace
    }

    pub fn key(&self, base: &Digest, edits: &[Edit], term: &TermId) -> CacheKey {
        CacheKey::new(&self.namespace, base, edits, term)
    }

    pub fn get(&mut self, key: &CacheKey) -> Option<f64> {
        match self.entries.get(key) {
            Some(value) => {
                self.hits += 1;
                Some(*value)
            }
            None => {
                self.misses += 1;
                None
            }
        }
    }

    pub fn put(&mut self, key: CacheKey, value: f64) {
        if self.entries.contains_key(&key) {
            return;
        }
        while self.entries.len() >= self.capacity {
            match self.order.pop_front() {
                Some(oldest) => {
                    self.entries.remove(&oldest);
                }
                None => break,
            }
        }
        self.order.push_back(key.clone());
        self.entries.insert(key, value);
    }

    pub fn hits(&self) -> u64 {
        self.hits
    }

    pub fn misses(&self) -> u64 {
        self.misses
    }

    pub fn reset_counters(&mut self) {
        self.hits = 0;
        self.misses = 0;
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
