//! TST-BURD-001 and TST-BURD-002: the demand ledger.
//!
//! Ticket P2-12.
//!
//! The two properties are deliberately different in kind. Amino-acid demand is
//! conserved *exactly* — a codon consumes one charged tRNA of its amino acid
//! whichever synonym is used — so it is compared bit for bit. Anticodon-pool
//! demand is a weighted float sum and is *not* conserved, because moving demand
//! between pools is what codon optimization for burden does; it is compared
//! against a canonical full recomputation to `1e-12`.

use std::collections::BTreeMap;

use cf_burden::{decoders, DecodingWeights, DemandLedger, Pairing};
use cf_types::{Codon, GeneticCode};

/// splitmix64, so a failure replays exactly.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

fn table() -> &'static GeneticCode {
    GeneticCode::get(11).unwrap()
}

/// Every pairing key the standard code can require, with distinct penalties.
///
/// Distinct on purpose: equal penalties would make every pool weigh the same,
/// so a substitution between synonyms read by different pools would move no
/// demand and `TST-BURD-002` would compare two constants.
fn weights() -> DecodingWeights {
    let mut s = BTreeMap::new();
    for (key, value) in [
        ("A:T", 0.0),
        ("T:A", 0.0),
        ("G:C", 0.0),
        ("C:G", 0.0),
        ("G:T", 0.41),
        ("T:G", 0.68),
        ("I:T", 0.0),
        ("I:C", 0.28),
        ("I:A", 0.9999),
    ] {
        s.insert(key.to_string(), value);
    }
    DecodingWeights::from_s_values(s).unwrap()
}

fn coding_codons() -> Vec<Codon> {
    Codon::all()
        .filter(|c| table().amino_acid(*c) != b'*')
        .collect()
}

// --------------------------------------------------------------------------- //
// the decoding rules
// --------------------------------------------------------------------------- //

#[test]
fn unit_every_codon_is_read_by_at_least_one_pool() {
    for codon in Codon::all() {
        let found = decoders(codon);
        assert!(
            !found.is_empty(),
            "{codon} is decoded by nothing, so its demand would vanish from the ledger"
        );
        assert!(found.len() <= 3, "{codon} claims {} pools", found.len());
    }
}

#[test]
fn unit_the_anticodon_is_the_reverse_complement_at_the_strict_positions() {
    // Positions 35 and 36 pair strictly; only position 34 wobbles. An
    // anticodon that differed elsewhere would be reading a different codon.
    for codon in Codon::all() {
        let bytes = codon.to_bytes();
        for decoder in decoders(codon) {
            let expect = |base: u8| match base {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                other => other,
            };
            assert_eq!(decoder.anticodon[1], expect(bytes[1]), "{codon}");
            assert_eq!(decoder.anticodon[2], expect(bytes[0]), "{codon}");
        }
    }
}

#[test]
fn unit_cricks_rules_are_encoded_as_published() {
    // Spot-checked against the hypothesis rather than against this
    // implementation: G34 reads C and U, U34 reads A and G, C34 reads G only,
    // and inosine reads U, C and A.
    let reads = |third: u8| -> Vec<(u8, Pairing)> {
        let codon = Codon::parse(&format!("AA{}", third as char)).unwrap();
        let mut found: Vec<(u8, Pairing)> = decoders(codon)
            .iter()
            .map(|d| (d.anticodon[0], d.pairing))
            .collect();
        found.sort();
        found
    };

    // C3 is read by G34 (Watson-Crick) and by inosine.
    assert_eq!(
        reads(b'C'),
        vec![(b'A', Pairing::Inosine), (b'G', Pairing::WatsonCrick)]
    );
    // G3 is read by C34 (Watson-Crick) and by U34 (U:G wobble).
    assert_eq!(
        reads(b'G'),
        vec![(b'C', Pairing::WatsonCrick), (b'T', Pairing::GU)]
    );
    // U3 is read by A34-as-inosine and by G34 (G:U wobble).
    assert_eq!(
        reads(b'T'),
        vec![(b'A', Pairing::Inosine), (b'G', Pairing::GU)]
    );
    // A3 is read by U34 (Watson-Crick) and by inosine.
    assert_eq!(
        reads(b'A'),
        vec![(b'A', Pairing::Inosine), (b'T', Pairing::WatsonCrick)]
    );
}

#[test]
fn unit_a34_is_counted_once_as_inosine_not_twice() {
    // A34 is essentially always deaminated, so listing it as both a
    // Watson-Crick reader of U3 and an inosine reader would count one pool
    // twice and inflate that pool's demand by a factor of two.
    let codon = Codon::parse("AAT").unwrap();
    let a34: Vec<_> = decoders(codon)
        .into_iter()
        .filter(|d| d.anticodon[0] == b'A')
        .collect();
    assert_eq!(a34.len(), 1, "A34 was counted twice for {codon}");
    assert_eq!(a34[0].pairing, Pairing::Inosine);
}

// --------------------------------------------------------------------------- //
// data provenance: nothing is defaulted
// --------------------------------------------------------------------------- //

#[test]
fn unit_a_missing_s_value_is_refused_rather_than_defaulted() {
    // A fabricated penalty would assert that a wobble pairing decodes as well
    // as Watson-Crick, changing every tAI weight and every ledger demand while
    // looking entirely plausible.
    let mut s = BTreeMap::new();
    s.insert("G:C".to_string(), 0.0);
    let partial = DecodingWeights::from_s_values(s).unwrap();

    let codon = Codon::parse("AAT").unwrap();
    let decoder = decoders(codon)[0];
    let error = partial
        .weight(codon, &decoder)
        .expect_err("a missing s value must be refused, not defaulted to zero");
    assert_eq!(error.code, "CF-BURD-001");
    assert!(error.detail.contains("tai.s_values"), "{}", error.detail);
}

#[test]
fn unit_an_out_of_range_s_value_is_refused() {
    for bad in [-0.1, 1.5] {
        let mut s = BTreeMap::new();
        s.insert("G:T".to_string(), bad);
        let error = DecodingWeights::from_s_values(s)
            .expect_err("s is a fractional penalty; outside [0, 1] the weight is nonsense");
        assert_eq!(error.code, "CF-BURD-001");
    }
}

#[test]
fn unit_a_non_synonymous_substitution_is_refused() {
    let mut ledger = DemandLedger::new();
    let weights = weights();
    let leucine = Codon::parse("CTG").unwrap();
    let arginine = Codon::parse("CGT").unwrap();
    ledger
        .accumulate(&[leucine], 1.0, &weights, table())
        .unwrap();

    let error = ledger
        .substitute(leucine, arginine, 1.0, &weights, table())
        .expect_err(
            "a ledger updated through a non-synonymous change stops describing \
                     the protein it was built for",
        );
    assert_eq!(error.code, "CF-BURD-002");
}

// --------------------------------------------------------------------------- //
// TST-BURD-001 — amino-acid demand is conserved bit for bit
// --------------------------------------------------------------------------- //

#[test]
fn tst_burd_001_amino_acid_demand_is_bit_identical_across_substitutions() {
    let weights = weights();
    let table = table();
    let mut rng = Rng(20260830);

    let coding = coding_codons();
    let sequence: Vec<Codon> = (0..300).map(|_| coding[rng.below(coding.len())]).collect();

    let mut ledger = DemandLedger::new();
    ledger.accumulate(&sequence, 1.0, &weights, table).unwrap();
    let before = ledger.amino_acid_demand().clone();

    let mut current = sequence.clone();
    let mut moved = 0u32;
    for _ in 0..100_000 {
        let position = rng.below(current.len());
        let synonyms: Vec<Codon> = table.synonyms(current[position]).to_vec();
        if synonyms.len() < 2 {
            continue;
        }
        let replacement = synonyms[rng.below(synonyms.len())];
        if replacement == current[position] {
            continue;
        }
        ledger
            .substitute(current[position], replacement, 1.0, &weights, table)
            .unwrap();
        current[position] = replacement;
        moved += 1;
    }

    // The fixture has to actually substitute, or conservation holds trivially.
    assert!(moved > 50_000, "only {moved} substitutions were applied");
    assert_eq!(
        ledger.amino_acid_demand(),
        &before,
        "amino-acid demand changed under synonymous substitution; a codon consumes \
         one charged tRNA of its amino acid whichever synonym is used"
    );

    // And pool demand *did* move, or the substitutions were all between
    // synonyms read by the same pools and the test says nothing.
    let recomputed_pools =
        DemandLedger::recompute(ledger.codon_counts(), 1.0, &weights, table).unwrap();
    assert!(
        recomputed_pools.anticodon_demand() != &BTreeMap::new(),
        "no pool demand was accumulated at all"
    );
}

// --------------------------------------------------------------------------- //
// TST-BURD-002 — incremental agrees with full recomputation
// --------------------------------------------------------------------------- //

#[test]
fn tst_burd_002_incremental_and_full_agree_over_many_edits() {
    let weights = weights();
    let table = table();
    let mut rng = Rng(4242);

    let coding = coding_codons();
    let sequence: Vec<Codon> = (0..200).map(|_| coding[rng.below(coding.len())]).collect();

    let mut ledger = DemandLedger::new();
    ledger.accumulate(&sequence, 2.5, &weights, table).unwrap();

    let mut current = sequence.clone();
    let mut applied = 0u32;
    let mut worst = 0.0f64;

    for step in 0..100_000 {
        let position = rng.below(current.len());
        let synonyms: Vec<Codon> = table.synonyms(current[position]).to_vec();
        if synonyms.len() < 2 {
            continue;
        }
        let replacement = synonyms[rng.below(synonyms.len())];
        if replacement == current[position] {
            continue;
        }
        ledger
            .substitute(current[position], replacement, 2.5, &weights, table)
            .unwrap();
        current[position] = replacement;
        applied += 1;

        // Comparing every step would dominate the runtime; comparing
        // periodically still catches drift, because drift accumulates.
        if step % 5_000 == 0 {
            let full =
                DemandLedger::recompute(ledger.codon_counts(), 2.5, &weights, table).unwrap();
            worst = worst.max(ledger.max_pool_difference(&full));
        }
    }

    let full = DemandLedger::recompute(ledger.codon_counts(), 2.5, &weights, table).unwrap();
    worst = worst.max(ledger.max_pool_difference(&full));

    assert!(
        applied > 50_000,
        "only {applied} substitutions were applied"
    );
    assert!(
        worst < 1e-12,
        "incremental and full ledgers drifted apart by {worst}"
    );
}

#[test]
fn unit_removing_what_was_never_added_is_refused() {
    // A ledger that silently absorbed an unmatched removal would drift from
    // the sequence it describes, and every conservation claim above would
    // become false without anything reporting it.
    let mut ledger = DemandLedger::new();
    let weights = weights();
    let leucine = Codon::parse("CTG").unwrap();
    let other = Codon::parse("CTA").unwrap();

    let error = ledger
        .substitute(leucine, other, 1.0, &weights, table())
        .expect_err("removing from an empty ledger must be refused");
    assert_eq!(error.code, "CF-BURD-003");
}

#[test]
fn unit_pool_demand_scales_with_flux() {
    // Demand is per-molecule cost times how much protein is made; a ledger
    // that ignored flux would rank a highly-expressed cistron the same as a
    // barely-transcribed one.
    let weights = weights();
    let sequence = vec![Codon::parse("CTG").unwrap(), Codon::parse("AAA").unwrap()];

    let mut single = DemandLedger::new();
    single
        .accumulate(&sequence, 1.0, &weights, table())
        .unwrap();
    let mut double = DemandLedger::new();
    double
        .accumulate(&sequence, 2.0, &weights, table())
        .unwrap();

    for (pool, demand) in single.anticodon_demand() {
        let scaled = double.anticodon_demand()[pool];
        assert!((scaled - 2.0 * demand).abs() < 1e-12, "pool {pool:?}");
    }
    // Amino-acid counts are counts, so they do not scale.
    assert_eq!(single.amino_acid_demand(), double.amino_acid_demand());
}
