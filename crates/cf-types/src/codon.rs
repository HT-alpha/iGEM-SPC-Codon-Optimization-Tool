//! Codon indexing.
//!
//! Ticket P1-05. Contract: `planning/02-schemas/sequence.schema.json`.
//!
//! The index convention is fixed once, here, and everything downstream depends
//! on it:
//!
//! ```text
//! index = 16*b0 + 4*b1 + b2,   with A=0, C=1, G=2, T=3
//! ```
//!
//! A `Codon` is a validated index in `0..64`, so an out-of-range codon is not
//! representable and table lookups need no bounds check on the hot path.

use crate::error::{code, CfError, Result};

/// Bases in index order. The order is load-bearing; changing it changes every
/// vendored genetic-code table.
pub const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Number of distinct codons.
pub const CODON_COUNT: usize = 64;

/// A codon, stored as its index in `0..64`.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Codon(u8);

impl Codon {
    /// Construct from an index, rejecting anything outside `0..64`.
    #[inline]
    pub fn from_index(index: u8) -> Result<Self> {
        if (index as usize) < CODON_COUNT {
            Ok(Self(index))
        } else {
            Err(CfError::new(
                code::SEQ_INDEX_OUT_OF_RANGE,
                format!("codon index {index} is not in 0..64"),
            ))
        }
    }

    /// Construct from an index already known to be in range.
    ///
    /// Used only where the index provably comes from `Codon::index` or from a
    /// bounded loop; it saturates rather than wrapping so a logic error becomes
    /// a wrong answer in one codon rather than an out-of-bounds read.
    #[inline]
    pub(crate) fn from_index_unchecked(index: u8) -> Self {
        debug_assert!((index as usize) < CODON_COUNT);
        Self(index & 0x3f)
    }

    /// Parse three ASCII bases.
    #[inline]
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        if bytes.len() != 3 {
            return Err(CfError::new(
                code::SEQ_NOT_IN_FRAME,
                format!("expected 3 bases, got {}", bytes.len()),
            ));
        }
        let mut index = 0u8;
        for &b in bytes {
            let value = base_index(b).ok_or_else(|| {
                CfError::new(
                    code::SEQ_INVALID_BASE,
                    format!("{:?} is not one of A, C, G, T", b as char),
                )
            })?;
            index = (index << 2) | value;
        }
        Ok(Self(index))
    }

    /// Parse three ASCII bases from a string.
    ///
    /// Named `parse` rather than `from_str` so that the inherent method does
    /// not shadow [`std::str::FromStr`], which is also implemented.
    pub fn parse(text: &str) -> Result<Self> {
        Self::from_bytes(text.as_bytes())
    }

    #[inline]
    pub fn index(self) -> u8 {
        self.0
    }

    #[inline]
    pub fn index_usize(self) -> usize {
        self.0 as usize
    }

    /// The three bases, in 5'->3' order.
    #[inline]
    pub fn to_bytes(self) -> [u8; 3] {
        [
            BASES[((self.0 >> 4) & 3) as usize],
            BASES[((self.0 >> 2) & 3) as usize],
            BASES[(self.0 & 3) as usize],
        ]
    }

    pub fn to_string_upper(self) -> String {
        String::from_utf8(self.to_bytes().to_vec()).expect("ASCII")
    }

    /// Every codon, in index order.
    pub fn all() -> impl Iterator<Item = Codon> {
        (0..CODON_COUNT as u8).map(Codon)
    }
}

/// Map an ASCII base to its index. Accepts lower case; `U`/`u` is transcoded to
/// `T` here so that no downstream module needs an RNA branch.
#[inline]
pub fn base_index(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

impl std::str::FromStr for Codon {
    type Err = CfError;
    fn from_str(text: &str) -> Result<Self> {
        Self::from_bytes(text.as_bytes())
    }
}

impl std::fmt::Display for Codon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let bytes = self.to_bytes();
        f.write_str(std::str::from_utf8(&bytes).expect("ASCII"))
    }
}

impl std::fmt::Debug for Codon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Codon({self}#{})", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// TST-CORE-002 (codec portion): the index codec is a bijection.
    #[test]
    fn unit_codon_index_bijection() {
        for index in 0..CODON_COUNT as u8 {
            let codon = Codon::from_index(index).unwrap();
            assert_eq!(codon.index(), index);
            let round = Codon::from_bytes(&codon.to_bytes()).unwrap();
            assert_eq!(round, codon, "round trip failed at index {index}");
        }
        let strings: std::collections::HashSet<String> =
            Codon::all().map(|c| c.to_string_upper()).collect();
        assert_eq!(strings.len(), CODON_COUNT, "codon strings are not distinct");
    }

    #[test]
    fn unit_index_convention_is_the_documented_one() {
        // index = 16*b0 + 4*b1 + b2 with A=0, C=1, G=2, T=3
        assert_eq!(Codon::parse("AAA").unwrap().index(), 0);
        assert_eq!(Codon::parse("AAC").unwrap().index(), 1);
        assert_eq!(Codon::parse("AAG").unwrap().index(), 2);
        assert_eq!(Codon::parse("AAT").unwrap().index(), 3);
        assert_eq!(Codon::parse("ACA").unwrap().index(), 4);
        assert_eq!(Codon::parse("CAA").unwrap().index(), 16);
        assert_eq!(Codon::parse("GAA").unwrap().index(), 32);
        assert_eq!(Codon::parse("TTT").unwrap().index(), 63);
        // ATG: A=0, T=3, G=2 -> 16*0 + 4*3 + 2
        assert_eq!(Codon::parse("ATG").unwrap().index(), 14);
    }

    #[test]
    fn unit_rna_is_transcoded_and_case_is_accepted() {
        assert_eq!(Codon::parse("AUG").unwrap(), Codon::parse("ATG").unwrap());
        assert_eq!(Codon::parse("aug").unwrap(), Codon::parse("ATG").unwrap());
    }

    #[test]
    fn unit_rejects_bad_input() {
        assert_eq!(Codon::parse("AT").unwrap_err().code, code::SEQ_NOT_IN_FRAME);
        assert_eq!(
            Codon::parse("ATGC").unwrap_err().code,
            code::SEQ_NOT_IN_FRAME
        );
        assert_eq!(
            Codon::parse("ATN").unwrap_err().code,
            code::SEQ_INVALID_BASE
        );
        assert_eq!(
            Codon::from_index(64).unwrap_err().code,
            code::SEQ_INDEX_OUT_OF_RANGE
        );
    }
}
