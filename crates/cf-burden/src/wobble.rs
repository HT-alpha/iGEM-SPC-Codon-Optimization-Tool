//! Codon → anticodon decoding: which tRNA pools can read which codon.
//!
//! Ticket P2-12. Specification: `05-burden-penalty-and-capacity.md`, and
//! `host_profile.schema.json`'s `trna` and `tai` blocks.
//!
//! ## What is a rule here and what is data
//!
//! This module encodes *rules* and refuses to invent *data*, and the line is
//! worth stating because it is the whole reason the module is shaped this way.
//!
//! **Rules** — which anticodon can pair with which codon — are structural
//! biochemistry, not measurements: Crick's wobble hypothesis (Crick, *J. Mol.
//! Biol.* 19:548, 1966) plus the inosine extension. They are the same for every
//! organism, in the way the genetic code table already in `cf-types` is, and
//! they are encoded here.
//!
//! **Data** — how *efficiently* each pairing decodes, and how many tRNA genes an
//! organism has — is measured, organism-specific, and comes from the host
//! profile: `tai.s_values` for the pairing penalties (dos Reis et al., *NAR*
//! 32:5036, 2004) and `trna.gene_copy_number` for the pools. None of it is
//! defaulted. A host profile missing an `s_value` for a pairing its own codons
//! require is an error naming the pairing, not a silent zero — a fabricated
//! penalty would change every tAI weight and every ledger demand while looking
//! entirely plausible.

use std::collections::BTreeMap;

use cf_types::{code, CfError, Codon, Result};

/// The named decoding table. `host_profile.trna.wobble_rules_id`.
pub const CRICK_WOBBLE_V1: &str = "crick_wobble_v1";

/// How an anticodon's wobble base pairs with a codon's third base.
///
/// The key is `anticodon position 34` against `codon position 3`, which is the
/// only position where non-Watson–Crick pairing is permitted; positions 35 and
/// 36 pair strictly.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Pairing {
    /// Standard Watson–Crick.
    WatsonCrick,
    /// G34 reading U3, or U34 reading G3.
    GU,
    /// Inosine (deaminated A34) reading U3, C3 or A3.
    Inosine,
}

impl Pairing {
    /// The key used to look up a penalty in `tai.s_values`.
    ///
    /// Written as the pairing the *anticodon* makes with the *codon*, which is
    /// the convention `host_profile.schema.json` documents with its `'G:U'` and
    /// `'I:C'` examples.
    pub fn s_value_key(&self, anticodon_base: u8, codon_base: u8) -> String {
        match self {
            Self::WatsonCrick => format!("{}:{}", anticodon_base as char, codon_base as char),
            Self::GU => format!("{}:{}", anticodon_base as char, codon_base as char),
            Self::Inosine => format!("I:{}", codon_base as char),
        }
    }
}

/// The anticodon that reads a codon, and how it pairs at the wobble position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Decoder {
    /// Anticodon 5'→3', which is the reverse complement of the codon at
    /// positions 35 and 36 and the wobble base at 34.
    pub anticodon: [u8; 3],
    pub pairing: Pairing,
}

impl Decoder {
    pub fn as_str(&self) -> String {
        String::from_utf8_lossy(&self.anticodon).into_owned()
    }
}

/// Watson–Crick complement.
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        other => other,
    }
}

/// Which anticodon wobble bases can read a given codon third base.
///
/// Crick's rules, plus inosine. Written as a match rather than a table so each
/// arm can carry the reason it exists:
///
/// * **G34** reads C3 (Watson–Crick) and U3 (G·U wobble).
/// * **U34** reads A3 (Watson–Crick) and G3 (U·G wobble).
/// * **C34** reads G3 only.
/// * **A34** is essentially always deaminated to inosine in the tRNAs that
///   have it, so it is *not* listed as reading U3 by Watson–Crick pairing.
///   Listing it would double-count the same pool: the gene is annotated A34
///   and the molecule decodes as I34.
/// * **I34** reads U3, C3 and A3.
///
/// A codon is therefore read by one, two or three pools.
fn wobble_bases_reading(codon_third: u8) -> Vec<(u8, Pairing)> {
    match codon_third {
        b'T' => vec![
            (b'A', Pairing::WatsonCrick),
            (b'G', Pairing::GU),
            (b'A', Pairing::Inosine),
        ],
        b'C' => vec![(b'G', Pairing::WatsonCrick), (b'A', Pairing::Inosine)],
        b'A' => vec![(b'T', Pairing::WatsonCrick), (b'A', Pairing::Inosine)],
        b'G' => vec![(b'C', Pairing::WatsonCrick), (b'T', Pairing::GU)],
        _ => Vec::new(),
    }
}

/// The one pool a codon is charged to when the host profile has no `s_values`.
///
/// Watson–Crick if the wobble table still lists it, otherwise G·U, otherwise
/// inosine. Charging the codon to every permitted pool equally would assert
/// that an unmeasured wobble pairing decodes as well as Watson–Crick. That
/// assertion is a fabricated efficiency, so it is not made. A host that has
/// `tai.s_values` should use [`DecodingWeights`] instead of this priority.
pub fn primary_decoder(codon: Codon) -> Option<Decoder> {
    let mut found = decoders(codon);
    found.sort_by_key(|decoder| {
        (
            match decoder.pairing {
                Pairing::WatsonCrick => 0,
                Pairing::GU => 1,
                Pairing::Inosine => 2,
            },
            decoder.anticodon,
        )
    });
    found.into_iter().next()
}

/// The tRNA pools that can decode a codon.
///
/// Deduplicated and ordered, so two callers building the same ledger from the
/// same codons get the same accumulation order — the ledger is summed in
/// floating point, and a differently-ordered sum is a different number.
pub fn decoders(codon: Codon) -> Vec<Decoder> {
    let bytes = codon.to_bytes();
    // Anticodon 5'→3' is the reverse complement of the codon: position 34
    // pairs with codon position 3.
    let position_35 = complement(bytes[1]);
    let position_36 = complement(bytes[0]);

    let mut found: Vec<Decoder> = wobble_bases_reading(bytes[2])
        .into_iter()
        .map(|(wobble, pairing)| Decoder {
            anticodon: [wobble, position_35, position_36],
            pairing,
        })
        .collect();

    // A34 appears twice above for `T` — once as Watson–Crick, once as
    // inosine — and the inosine reading is the one that happens. Keeping both
    // would count one pool twice.
    found.retain(|d| !(d.anticodon[0] == b'A' && d.pairing == Pairing::WatsonCrick));
    found.sort_by_key(|d| (d.anticodon, d.pairing));
    found.dedup();
    found
}

/// A host's decoding weights: how well each pool reads each codon.
///
/// Built from the host profile's `s_values`. Every weight is `1 − s(pairing)`,
/// so a Watson–Crick pairing with `s = 0` weighs 1 and a penalised wobble
/// weighs less.
#[derive(Clone, Debug, Default)]
pub struct DecodingWeights {
    /// `s` by pairing key, exactly as the host profile supplied it.
    s_values: BTreeMap<String, f64>,
}

impl DecodingWeights {
    /// Read the `tai.s_values` block.
    ///
    /// Missing keys are *not* defaulted here — the error is raised at use, so
    /// it can name the codon and pairing that needed it rather than a key the
    /// caller has no context for.
    pub fn from_s_values(s_values: BTreeMap<String, f64>) -> Result<Self> {
        for (key, value) in &s_values {
            if !(0.0..=1.0).contains(value) {
                return Err(CfError::new(
                    code::BURD_WOBBLE_PARAMETER,
                    format!(
                        "s value for pairing {key} is {value}, outside [0, 1]; s is a \
                         fractional penalty and a value outside it makes the decoding \
                         weight negative or greater than one"
                    ),
                ));
            }
        }
        Ok(Self { s_values })
    }

    /// The weight with which `decoder` reads `codon`.
    ///
    /// Refuses rather than defaults. A missing `s` value is a gap in the host
    /// profile, and substituting zero would silently assert that a wobble
    /// pairing decodes as efficiently as Watson–Crick — a fabricated constant
    /// that changes every tAI weight and every ledger demand while looking
    /// entirely plausible.
    pub fn weight(&self, codon: Codon, decoder: &Decoder) -> Result<f64> {
        let key = decoder
            .pairing
            .s_value_key(decoder.anticodon[0], codon.to_bytes()[2]);
        match self.s_values.get(&key) {
            Some(s) => Ok(1.0 - s),
            None => Err(CfError::new(
                code::BURD_WOBBLE_PARAMETER,
                format!(
                    "the host profile has no s value for pairing {key}, which codon {codon} \
                     requires; supply it in tai.s_values rather than accepting a default, \
                     because a fabricated penalty changes every tAI weight and every \
                     ledger demand"
                ),
            )),
        }
    }

    pub fn contains(&self, key: &str) -> bool {
        self.s_values.contains_key(key)
    }
}
