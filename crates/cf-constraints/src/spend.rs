//! Wall-clock and resident-memory caps.
//!
//! `DEP-BUDG-002`. The numbers are the published ADR-0008 walls and the
//! laptop profile's 8 GB, computed by the caller. This module only refuses
//! when a cap that was actually supplied has already been passed. Absent
//! fields are unlimited, so a caller that does not set them is unchanged.

use std::cell::Cell;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use cf_types::{code, CfError, Result};

/// How long a resident-set reading may be reused.
///
/// `enforce` sits on the seed-repair inner loop. Parsing `/proc/self/status`
/// on every call is a syscall the deadline check does not need. A cap still
/// trips within this interval, and a deadline is still checked on every call.
const RSS_SAMPLE_INTERVAL: Duration = Duration::from_millis(200);

thread_local! {
    static RSS_SAMPLE: Cell<Option<(Instant, Option<u64>)>> = const { Cell::new(None) };
}

/// A cap the caller measured. `None` means that axis is not limited.
#[derive(Clone, Copy, Debug, Default)]
pub struct SpendLimit {
    pub deadline: Option<Instant>,
    pub max_rss_bytes: Option<u64>,
}

impl SpendLimit {
    /// `deadline_unix_ms` is milliseconds since the Unix epoch.
    ///
    /// A deadline that is already past becomes an instant slightly before
    /// now, so the next [`Self::enforce`] refuses instead of racing the clock.
    pub fn from_caps(deadline_unix_ms: Option<u64>, max_rss_bytes: Option<u64>) -> Self {
        let deadline = deadline_unix_ms.map(|deadline_ms| {
            let now_ms = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .map(|elapsed| u64::try_from(elapsed.as_millis()).unwrap_or(u64::MAX))
                .unwrap_or(0);
            if deadline_ms <= now_ms {
                Instant::now()
                    .checked_sub(Duration::from_millis(1))
                    .unwrap_or_else(Instant::now)
            } else {
                Instant::now() + Duration::from_millis(deadline_ms - now_ms)
            }
        });
        Self {
            deadline,
            max_rss_bytes,
        }
    }

    pub fn enforce(self) -> Result<()> {
        if let Some(deadline) = self.deadline {
            if Instant::now() >= deadline {
                return Err(CfError::new(
                    code::BUDGET_EXCEEDED,
                    "the phase exceeded its wall-clock budget",
                ));
            }
        }
        if let Some(cap) = self.max_rss_bytes {
            if let Some(rss) = resident_bytes_cached() {
                if rss > cap {
                    return Err(CfError::new(
                        code::BUDGET_EXCEEDED,
                        format!(
                            "resident set is {rss} bytes, which is above the memory cap of {cap} bytes"
                        ),
                    ));
                }
            }
        }
        Ok(())
    }
}

fn resident_bytes_cached() -> Option<u64> {
    let now = Instant::now();
    if let Some((at, bytes)) = RSS_SAMPLE.get() {
        if now.saturating_duration_since(at) < RSS_SAMPLE_INTERVAL {
            return bytes;
        }
    }
    let bytes = resident_bytes();
    RSS_SAMPLE.set(Some((now, bytes)));
    bytes
}

/// `VmRSS` from `/proc/self/status`, in bytes.
///
/// `None` when the file is absent, so a platform without it does not invent
/// a reading and does not refuse.
fn resident_bytes() -> Option<u64> {
    let text = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("VmRSS:") else {
            continue;
        };
        let kb = rest.split_whitespace().next()?.parse::<u64>().ok()?;
        return Some(kb.saturating_mul(1024));
    }
    None
}
