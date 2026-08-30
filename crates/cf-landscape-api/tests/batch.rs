//! TST-LM-001 and the batch-yield protocol.
//!
//! Ticket P2-05.
//!
//! `TST-LM-001` — "`score_batch` equals a loop over `score` for every batch
//! size including 1 and `E_max`" — is applied to a stub, which is what the
//! ticket asks for and is the only honest thing available: the L0 terms it will
//! eventually guard arrive at `P2-11` and `P3-05`. The stub is written to be
//! *batch-sensitive on purpose*, so a batched implementation that quietly
//! differed from its per-item one would be caught rather than agreeing by
//! construction.

use cf_landscape_api::batch::{BatchProtocol, BatchResults, BatchScorer, Yield};
use cf_types::{CodonSequence, Result, TermId};

const E_MAX: usize = 64;

fn sequence(cds: &str) -> CodonSequence {
    CodonSequence::from_cds(cds, 1).unwrap()
}

fn corpus(count: usize) -> Vec<CodonSequence> {
    // Distinct sequences, so a scorer that returned the same value for every
    // item would be caught by the comparisons below rather than agreeing with
    // itself.
    let synonyms = ["CTG", "CTA", "CTC", "CTT", "TTA", "TTG"];
    (0..count)
        .map(|index| {
            let mut cds = String::from("ATG");
            for position in 0..4 {
                cds.push_str(synonyms[(index + position) % synonyms.len()]);
            }
            cds.push_str("TAA");
            sequence(&cds)
        })
        .collect()
}

/// A stub L0 term.
///
/// Its score depends only on the item, but it is computed through a genuinely
/// batched path — a running accumulator over the batch — so a bug that let
/// batch position leak into the value would show up. A stub that returned a
/// constant, or that looped internally, would make `TST-LM-001` vacuous: the
/// batched and per-item paths would be the same code.
#[derive(Debug)]
struct StubL0 {
    term: TermId,
}

impl StubL0 {
    fn new() -> Self {
        Self {
            term: TermId::new("lm.naturalness").unwrap(),
        }
    }
}

impl BatchScorer for StubL0 {
    fn term(&self) -> &TermId {
        &self.term
    }

    fn score_batch(&self, items: &[CodonSequence]) -> Result<Vec<Option<f64>>> {
        // A shared accumulator across the batch, reset per item. Written this
        // way because it is the shape a real batched model has, and it is where
        // a cross-item leak would live — `LeakyL0` below is the same code with
        // the reset removed, which is what makes `TST-LM-001` able to fail.
        #[allow(unused_assignments)]
        let mut scratch = 0.0f64;
        let mut out = Vec::with_capacity(items.len());
        for item in items {
            scratch = 0.0;
            for (position, codon) in item.codons().iter().enumerate() {
                scratch += (codon.index_usize() as f64) / (position as f64 + 1.0);
            }
            out.push(Some(scratch));
        }
        Ok(out)
    }
}

/// The same stub with a deliberate cross-item leak, for the negative control.
#[derive(Debug)]
struct LeakyL0 {
    term: TermId,
}

impl BatchScorer for LeakyL0 {
    fn term(&self) -> &TermId {
        &self.term
    }
    fn score_batch(&self, items: &[CodonSequence]) -> Result<Vec<Option<f64>>> {
        let mut scratch = 0.0f64;
        let mut out = Vec::with_capacity(items.len());
        for item in items {
            for (position, codon) in item.codons().iter().enumerate() {
                scratch += (codon.index_usize() as f64) / (position as f64 + 1.0);
            }
            out.push(Some(scratch));
        }
        Ok(out)
    }
}

// --------------------------------------------------------------------------- //
// TST-LM-001
// --------------------------------------------------------------------------- //

#[test]
fn tst_lm_001_score_batch_equals_a_loop_over_score_at_every_batch_size() {
    let scorer = StubL0::new();
    let items = corpus(E_MAX);

    let one_at_a_time: Vec<Option<f64>> = items
        .iter()
        .map(|item| scorer.score_one(item).unwrap())
        .collect();

    // Every size from 1 to E_MAX, not a sample: the sizes where batching goes
    // wrong are the boundaries, and a sample is exactly what misses them.
    for size in 1..=E_MAX {
        let mut batched: Vec<Option<f64>> = Vec::new();
        for chunk in items.chunks(size) {
            batched.extend(scorer.score_batch(chunk).unwrap());
        }
        assert_eq!(batched.len(), one_at_a_time.len());
        for (index, (a, b)) in batched.iter().zip(&one_at_a_time).enumerate() {
            match (a, b) {
                (Some(a), Some(b)) => assert!(
                    (a - b).abs() < 1e-12,
                    "batch size {size}, item {index}: {a} vs {b}"
                ),
                (None, None) => {}
                _ => panic!("batch size {size}, item {index}: one path returned UNKNOWN"),
            }
        }
    }
}

#[test]
fn tst_lm_001_the_fixture_can_actually_detect_a_batch_dependent_scorer() {
    // Without this, the test above passes for a stub whose batched and
    // per-item paths are the same code, and asserts nothing about batching.
    let leaky = LeakyL0 {
        term: TermId::new("lm.naturalness").unwrap(),
    };
    let items = corpus(8);

    let one_at_a_time: Vec<Option<f64>> = items
        .iter()
        .map(|item| leaky.score_one(item).unwrap())
        .collect();
    let batched = leaky.score_batch(&items).unwrap();

    assert_ne!(
        batched, one_at_a_time,
        "the leaky scorer agreed with itself, so TST-LM-001 above would pass for a \
         genuinely batch-dependent implementation too"
    );
}

// --------------------------------------------------------------------------- //
// one crossing per generation
// --------------------------------------------------------------------------- //

#[test]
fn p2_05_a_generation_yields_exactly_one_batch() {
    // ADR-0001's cap, asserted rather than assumed. The protocol refuses a
    // second `step` while a batch is outstanding, which is what makes "one
    // crossing per generation" a property of the type rather than a discipline
    // the caller has to remember.
    let mut protocol = BatchProtocol::new();
    let candidates: Vec<usize> = (0..10).collect();

    let (first, chosen) = protocol.step(&candidates, E_MAX).unwrap();
    assert_eq!(first, Yield::Score(0, 10));
    assert_eq!(chosen.len(), 10);

    let error = protocol
        .step(&candidates, E_MAX)
        .expect_err("a second step in one generation is a second crossing");
    assert_eq!(error.code, "CF-LAND-011");

    protocol
        .absorb(&BatchResults {
            generation: 0,
            values: vec![Some(1.0); 10],
        })
        .unwrap();
    protocol.advance().unwrap();

    // The next generation may step again, and does so from generation 1.
    let (second, _) = protocol.step(&candidates, E_MAX).unwrap();
    assert_eq!(second, Yield::Score(1, 10));
}

#[test]
fn p2_05_advancing_with_a_batch_outstanding_is_refused() {
    let mut protocol = BatchProtocol::new();
    protocol.step(&[0, 1, 2], E_MAX).unwrap();
    let error = protocol
        .advance()
        .expect_err("advancing would absorb the results into the wrong generation");
    assert_eq!(error.code, "CF-LAND-011");
}

// --------------------------------------------------------------------------- //
// rationing: truncation defers rather than discards
// --------------------------------------------------------------------------- //

#[test]
fn p2_05_truncation_defers_the_remainder_to_the_next_generation() {
    // §6 component (d). A bounded budget must delay an evaluation, not drop
    // it: an individual dropped for budget is `UNKNOWN` forever, and nothing
    // in the output says it was ever a candidate.
    let mut protocol = BatchProtocol::new();
    let candidates: Vec<usize> = (0..10).collect();

    let (_, chosen) = protocol.step(&candidates, 4).unwrap();
    assert_eq!(chosen, vec![0, 1, 2, 3]);
    assert_eq!(protocol.deferred(), &[4, 5, 6, 7, 8, 9]);

    protocol
        .absorb(&BatchResults {
            generation: 0,
            values: vec![Some(0.0); 4],
        })
        .unwrap();
    protocol.advance().unwrap();

    // The deferred items go to the head of the next queue, ahead of new ones.
    let (_, next) = protocol.step(&[100, 101], 4).unwrap();
    assert_eq!(next, vec![4, 5, 6, 7]);
    assert_eq!(protocol.deferred(), &[8, 9, 100, 101]);
}

#[test]
fn p2_05_an_item_deferred_and_reproposed_is_not_scored_twice() {
    let mut protocol = BatchProtocol::new();
    let (_, _) = protocol.step(&[0, 1, 2, 3], 2).unwrap();
    protocol
        .absorb(&BatchResults {
            generation: 0,
            values: vec![Some(0.0); 2],
        })
        .unwrap();
    protocol.advance().unwrap();

    // 2 and 3 are deferred; proposing them again must not duplicate them.
    let (_, next) = protocol.step(&[2, 3, 4], 10).unwrap();
    assert_eq!(next, vec![2, 3, 4]);
    assert!(protocol.deferred().is_empty());
}

#[test]
fn p2_05_an_empty_queue_yields_ready_rather_than_an_empty_batch() {
    let mut protocol = BatchProtocol::new();
    let (outcome, chosen) = protocol.step(&[], E_MAX).unwrap();
    assert_eq!(outcome, Yield::Ready);
    assert!(chosen.is_empty());
    assert!(
        !protocol.is_waiting(),
        "an empty batch left the protocol waiting"
    );
}

// --------------------------------------------------------------------------- //
// absorb refuses what it cannot match
// --------------------------------------------------------------------------- //

#[test]
fn p2_05_results_of_the_wrong_length_are_refused_not_aligned() {
    // Aligning by position would attribute one sequence's score to another,
    // and the scores would be plausible - just wrong, with nothing downstream
    // able to tell.
    let mut protocol = BatchProtocol::new();
    protocol.step(&[0, 1, 2], E_MAX).unwrap();

    let error = protocol
        .absorb(&BatchResults {
            generation: 0,
            values: vec![Some(1.0), Some(2.0)],
        })
        .expect_err("two results for a three-item batch must be refused");
    assert_eq!(error.code, "CF-LAND-012");
    assert!(
        protocol.is_waiting(),
        "a refused absorb cleared the batch anyway"
    );
}

#[test]
fn p2_05_results_from_another_generation_are_refused() {
    let mut protocol = BatchProtocol::new();
    protocol.step(&[0, 1], E_MAX).unwrap();
    let error = protocol
        .absorb(&BatchResults {
            generation: 7,
            values: vec![Some(1.0), Some(2.0)],
        })
        .expect_err("late results from another generation must not be applied");
    assert_eq!(error.code, "CF-LAND-012");
}

#[test]
fn p2_05_absorbing_twice_is_refused() {
    let mut protocol = BatchProtocol::new();
    protocol.step(&[0, 1], E_MAX).unwrap();
    let results = BatchResults {
        generation: 0,
        values: vec![Some(1.0), Some(2.0)],
    };
    protocol.absorb(&results).unwrap();
    let error = protocol
        .absorb(&results)
        .expect_err("a duplicated absorb would double-count the evaluations");
    assert_eq!(error.code, "CF-LAND-012");
    assert_eq!(protocol.evaluations(), 2);
}

#[test]
fn p2_05_unknown_results_are_not_counted_as_evaluations() {
    // An item the scorer could not evaluate has not consumed the budget the
    // rationing exists to protect, and counting it would ration against work
    // that never happened.
    let mut protocol = BatchProtocol::new();
    protocol.step(&[0, 1, 2], E_MAX).unwrap();
    protocol
        .absorb(&BatchResults {
            generation: 0,
            values: vec![Some(1.0), None, Some(3.0)],
        })
        .unwrap();
    assert_eq!(protocol.evaluations(), 2);
    assert_eq!(protocol.absorbed(), &[Some(1.0), None, Some(3.0)]);
}

// --------------------------------------------------------------------------- //
// snapshot / resume, including mid-batch
// --------------------------------------------------------------------------- //

#[test]
fn p2_05_a_protocol_snapshotted_mid_batch_resumes_still_waiting() {
    // This is the state a checkpoint is *most* likely to be taken in: the
    // search is idle precisely while the scorer runs. A resume that forgot the
    // outstanding batch would step again and issue a second crossing for the
    // same generation.
    let mut protocol = BatchProtocol::new();
    protocol.step(&[0, 1, 2, 3, 4], 3).unwrap();

    let resumed = BatchProtocol::resume(&protocol.snapshot()).unwrap();
    assert_eq!(resumed, protocol);
    assert!(resumed.is_waiting());
    assert_eq!(resumed.deferred(), &[3, 4]);

    let mut resumed = resumed;
    let error = resumed
        .step(&[9], 3)
        .expect_err("a resumed protocol must still know a batch is outstanding");
    assert_eq!(error.code, "CF-LAND-011");
}

#[test]
fn p2_05_the_snapshot_round_trips_exactly() {
    let mut protocol = BatchProtocol::new();
    protocol.step(&[0, 1, 2], 2).unwrap();
    protocol
        .absorb(&BatchResults {
            generation: 0,
            values: vec![Some(0.1), None],
        })
        .unwrap();
    protocol.advance().unwrap();
    protocol.step(&[5, 6], 8).unwrap();

    let restored = BatchProtocol::resume(&protocol.snapshot()).unwrap();
    assert_eq!(restored, protocol);
    assert_eq!(restored.snapshot(), protocol.snapshot());
}

/// Adversarial floats through the snapshot.
///
/// The `serde_json` float-roundtrip defect was latent until one organically
/// generated objective value happened to expose it, and by then every
/// checkpoint in the system had been affected since the first float entered
/// one. Waiting for the next boundary case to arrive on its own is the same
/// bet a second time, so this fixture brings them deliberately: denormals, the
/// smallest and largest finite values, repeating binary fractions, and values
/// needing all 17 significant digits.
#[test]
fn p2_05_adversarial_floats_survive_the_snapshot_exactly() {
    let adversarial: Vec<Option<f64>> = vec![
        Some(0.0),
        Some(-0.0),
        Some(f64::MIN_POSITIVE),
        Some(f64::MIN_POSITIVE / 2.0), // subnormal
        Some(5e-324),                  // smallest subnormal
        Some(f64::MAX),
        Some(f64::MIN),
        Some(0.1),
        Some(0.2),
        Some(0.1 + 0.2),
        Some(1.0 / 3.0),
        Some(17.0 / 42.0),
        Some(0.404_761_904_761_904_77),
        Some(2.225_073_858_507_201_4e-308),
        Some(9_007_199_254_740_993.0), // 2^53 + 1, unrepresentable as written
        Some(1.797_693_134_862_315_7e308),
        None,
        Some(-1.0 / 3.0),
    ];

    let mut protocol = BatchProtocol::new();
    let indices: Vec<usize> = (0..adversarial.len()).collect();
    protocol.step(&indices, adversarial.len()).unwrap();
    protocol
        .absorb(&BatchResults {
            generation: 0,
            values: adversarial.clone(),
        })
        .unwrap();

    // Through *text*, not through a `serde_json::Value`. A snapshot that never
    // leaves Rust holds an `f64` directly and no parser is involved, so a
    // fixture that round-tripped the value alone would not exercise the code
    // that was broken: the defect was in parsing decimal text back to binary.
    // A checkpoint on disk is text, so text is what has to be tested.
    let serialized = serde_json::to_string(&protocol.snapshot()).unwrap();
    let reparsed: serde_json::Value = serde_json::from_str(&serialized).unwrap();
    let restored = BatchProtocol::resume(&reparsed).unwrap();

    for (index, (a, b)) in restored.absorbed().iter().zip(&adversarial).enumerate() {
        match (a, b) {
            (Some(a), Some(b)) => assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "value {index} changed across the snapshot: {a:.17e} vs {b:.17e}"
            ),
            (None, None) => {}
            _ => panic!("value {index}: UNKNOWN was gained or lost"),
        }
    }
    assert_eq!(restored, protocol);

    // And the fixture is genuinely adversarial: these values are not all
    // exactly representable in short decimal form, so a shortest-repr
    // round-trip that was not correctly rounded would differ.
    assert!(adversarial.iter().flatten().any(|v| format!("{v}")
        .parse::<f64>()
        .map(|p| p.to_bits())
        != Ok(v.to_bits())
        || format!("{v:.17e}").len() > 20));
}

#[test]
fn p2_05_a_malformed_snapshot_is_refused() {
    let mut protocol = BatchProtocol::new();
    protocol.step(&[0], 1).unwrap();
    let good = protocol.snapshot();

    for missing in [
        "generation",
        "outstanding",
        "absorbed",
        "evaluations",
        "deferred",
    ] {
        let mut broken = good.clone();
        broken.as_object_mut().unwrap().remove(missing);
        let error = BatchProtocol::resume(&broken)
            .expect_err("a snapshot missing a required field must not load");
        assert_eq!(error.code, "CF-LAND-012");
        assert!(error.detail.contains(missing), "{}", error.detail);
    }
}
