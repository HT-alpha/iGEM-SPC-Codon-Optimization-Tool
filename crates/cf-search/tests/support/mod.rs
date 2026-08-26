//! Fixtures shared by the `cf-search` integration tests.
//!
//! Each test binary compiles this module separately, so a fixture used by one
//! and not another is dead code from that binary's point of view.
#![allow(dead_code)]

use cf_landscape_api::HostContext;
use cf_search::CodonWeights;
use cf_types::{Digest, HostId};

pub fn host() -> HostContext {
    HostContext::new(
        HostId::new("escherichia_coli:511145/BL21").unwrap(),
        Digest::of_bytes(b"host-profile-fixture"),
    )
}

/// A deterministic, non-uniform weight table.
///
/// Non-uniform on purpose: a uniform table would make a weighted-sampling
/// assertion pass for the wrong reason.
pub fn weights() -> CodonWeights {
    let mut table = [0.0f64; 64];
    for (index, slot) in table.iter_mut().enumerate() {
        *slot = ((index as f64 * 17.0) % 23.0 + 1.0) / 23.0;
    }
    let weights = CodonWeights::from_relative_adaptiveness(table).unwrap();
    assert!(!weights.is_uniform());
    weights
}
