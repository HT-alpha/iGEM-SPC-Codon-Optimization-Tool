//! RFC 8785 canonical JSON (JCS) and content digests.
//!
//! Ticket P1-07. Contract: `planning/02-schemas/common.schema.json`
//! (`Sha256Digest`) and the digest fields of every document schema.
//!
//! WHY THIS IS NOT `serde_json::to_string`. A digest is only a stable identity
//! if two independent serializers of the same value produce the same bytes.
//! `serde_json` gives no such guarantee: map iteration order, float formatting
//! and escaping are all implementation choices. RFC 8785 pins all three:
//!
//!   * object keys are sorted by their UTF-16 code units, not by bytes;
//!   * numbers are serialized per ECMAScript `Number::toString`, which is
//!     *not* the same as Rust's `{}` once exponents get involved
//!     (`1e21` is `"1e+21"` in ECMAScript and `"1000000000000000000000"` in
//!     Rust's `Display`);
//!   * strings escape only what RFC 8259 requires, with the short forms.
//!
//! Getting any one of those wrong produces digests that agree on this machine
//! and disagree on the reviewer's, which would silently break every lockfile
//! and every cache key.

use sha2::{Digest as _, Sha256};

use crate::error::{code, CfError, Result};

/// A content digest, always carrying its algorithm prefix.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Digest(String);

impl Digest {
    /// Digest of raw bytes.
    pub fn of_bytes(bytes: &[u8]) -> Self {
        let mut hasher = Sha256::new();
        hasher.update(bytes);
        Self(format!("sha256:{:x}", hasher.finalize()))
    }

    /// Digest of a JSON value, over its RFC 8785 canonical form.
    pub fn of_json(value: &serde_json::Value) -> Result<Self> {
        Ok(Self::of_bytes(canonicalize(value)?.as_bytes()))
    }

    /// Digest of a document with one field removed.
    ///
    /// This is how every `digest` field in the schemas is defined: the digest
    /// of the document that carries it cannot include itself.
    pub fn of_document(value: &serde_json::Value, omit_field: &str) -> Result<Self> {
        match value {
            serde_json::Value::Object(map) => {
                let mut stripped = map.clone();
                stripped.remove(omit_field);
                Self::of_json(&serde_json::Value::Object(stripped))
            }
            _ => Self::of_json(value),
        }
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// Parse and validate a digest string.
    pub fn parse(text: &str) -> Result<Self> {
        let hex = text.strip_prefix("sha256:").ok_or_else(|| {
            CfError::new(
                code::SEQ_BAD_IDENTIFIER,
                format!("{text:?} lacks the sha256: prefix"),
            )
        })?;
        if hex.len() != 64
            || !hex
                .bytes()
                .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
        {
            return Err(CfError::new(
                code::SEQ_BAD_IDENTIFIER,
                format!("{text:?} is not sha256: followed by 64 lower-case hex digits"),
            ));
        }
        Ok(Self(text.to_string()))
    }
}

impl std::fmt::Display for Digest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::fmt::Debug for Digest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Abbreviated: a full digest in a panic message is unreadable.
        write!(f, "Digest({}…)", &self.0[..sha_prefix_len(&self.0)])
    }
}

fn sha_prefix_len(text: &str) -> usize {
    text.len().min(15)
}

/// Serialize a value to its RFC 8785 canonical form.
pub fn canonicalize(value: &serde_json::Value) -> Result<String> {
    let mut out = String::new();
    write_value(value, &mut out)?;
    Ok(out)
}

fn write_value(value: &serde_json::Value, out: &mut String) -> Result<()> {
    match value {
        serde_json::Value::Null => out.push_str("null"),
        serde_json::Value::Bool(true) => out.push_str("true"),
        serde_json::Value::Bool(false) => out.push_str("false"),
        serde_json::Value::Number(n) => out.push_str(&write_number(n)?),
        serde_json::Value::String(s) => write_string(s, out),
        serde_json::Value::Array(items) => {
            out.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write_value(item, out)?;
            }
            out.push(']');
        }
        serde_json::Value::Object(map) => {
            // RFC 8785: sort by the UTF-16 code units of the key. For ASCII
            // keys this coincides with byte order, but the spec is written in
            // UTF-16 and a key outside the BMP would order differently, so the
            // comparison is done properly rather than assumed away.
            let mut keys: Vec<(&String, Vec<u16>)> = map
                .keys()
                .map(|k| (k, k.encode_utf16().collect()))
                .collect();
            keys.sort_by(|a, b| a.1.cmp(&b.1));

            out.push('{');
            for (i, (key, _)) in keys.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write_string(key, out);
                out.push(':');
                write_value(&map[*key], out)?;
            }
            out.push('}');
        }
    }
    Ok(())
}

fn write_string(text: &str, out: &mut String) {
    out.push('"');
    for ch in text.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

/// Largest integer exactly representable as an IEEE-754 double.
const MAX_SAFE_INTEGER: i128 = 9_007_199_254_740_991;

fn write_number(n: &serde_json::Number) -> Result<String> {
    if let Some(v) = n.as_u64() {
        return integer_string(v as i128);
    }
    if let Some(v) = n.as_i64() {
        return integer_string(v as i128);
    }
    let v = n.as_f64().ok_or_else(|| {
        CfError::new(
            code::SEQ_NOT_CANONICAL,
            format!("number {n} is not representable"),
        )
    })?;
    es_number(v)
}

fn integer_string(v: i128) -> Result<String> {
    if v.abs() > MAX_SAFE_INTEGER {
        // JCS is defined over IEEE-754 doubles. An integer beyond 2^53 cannot
        // round-trip, so emitting one would produce a digest that a conforming
        // implementation in another language could not reproduce.
        return Err(CfError::new(
            code::SEQ_NOT_CANONICAL,
            format!("integer {v} exceeds the IEEE-754 safe integer range required by RFC 8785"),
        ));
    }
    Ok(v.to_string())
}

/// ECMAScript `Number::toString` (ECMA-262 6.1.6.1.20), as RFC 8785 requires.
fn es_number(v: f64) -> Result<String> {
    if v.is_nan() || v.is_infinite() {
        return Err(CfError::new(
            code::SEQ_NOT_CANONICAL,
            "NaN and infinity have no JSON representation",
        ));
    }
    if v == 0.0 {
        // Covers -0.0, which ECMAScript renders as "0".
        return Ok("0".to_string());
    }
    if v < 0.0 {
        return Ok(format!("-{}", es_number(-v)?));
    }

    // Rust's LowerExp gives the shortest round-tripping decimal, which is
    // exactly the (digits, exponent) pair the ECMAScript algorithm is defined
    // over.
    let sci = format!("{v:e}");
    let (mantissa, exponent) = sci
        .split_once('e')
        .expect("LowerExp always emits an exponent");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let exponent: i32 = exponent.parse().expect("LowerExp exponent is an integer");

    let k = digits.len() as i32;
    let n = exponent + 1; // value == 0.digits * 10^n

    let rendered = if k <= n && n <= 21 {
        // Integer with trailing zeros.
        let mut s = digits;
        for _ in 0..(n - k) {
            s.push('0');
        }
        s
    } else if 0 < n && n <= 21 {
        // Decimal point inside the digits.
        format!("{}.{}", &digits[..n as usize], &digits[n as usize..])
    } else if -6 < n && n <= 0 {
        // Leading "0." and (-n) zeros.
        let mut s = String::from("0.");
        for _ in 0..(-n) {
            s.push('0');
        }
        s.push_str(&digits);
        s
    } else {
        // Exponential form. ECMAScript writes a '+' for a positive exponent.
        let e = n - 1;
        let sign = if e >= 0 { '+' } else { '-' };
        if k == 1 {
            format!("{digits}e{sign}{}", e.abs())
        } else {
            format!("{}.{}e{sign}{}", &digits[..1], &digits[1..], e.abs())
        }
    };
    Ok(rendered)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// TST-CORE-003: canonical form is stable across key insertion order.
    #[test]
    fn unit_tst_core_003_key_order_does_not_matter() {
        let a: serde_json::Value =
            serde_json::from_str(r#"{"b":1,"a":2,"c":{"z":1,"y":2}}"#).unwrap();
        let b: serde_json::Value =
            serde_json::from_str(r#"{"c":{"y":2,"z":1},"a":2,"b":1}"#).unwrap();
        assert_eq!(canonicalize(&a).unwrap(), canonicalize(&b).unwrap());
        assert_eq!(
            canonicalize(&a).unwrap(),
            r#"{"a":2,"b":1,"c":{"y":2,"z":1}}"#
        );
        assert_eq!(Digest::of_json(&a).unwrap(), Digest::of_json(&b).unwrap());
    }

    /// TST-CORE-003: number formatting follows ECMAScript, not Rust's Display.
    ///
    /// The last case is written with more digits than an `f64` can hold on
    /// purpose: it is the classic RFC 8785 test vector for a value whose
    /// shortest round-tripping form has fewer digits than its literal.
    #[allow(clippy::excessive_precision)]
    #[test]
    fn unit_tst_core_003_es_number_formatting() {
        let cases: &[(f64, &str)] = &[
            (0.0, "0"),
            (-0.0, "0"),
            (1.0, "1"),
            (-1.0, "-1"),
            (0.1, "0.1"),
            (1.5, "1.5"),
            (100.0, "100"),
            (1e20, "100000000000000000000"),
            // The two boundary cases where Rust's Display and ECMAScript differ.
            (1e21, "1e+21"),
            (1e22, "1e+22"),
            (1e-6, "0.000001"),
            (1e-7, "1e-7"),
            (1.2345e-7, "1.2345e-7"),
            (f64::MAX, "1.7976931348623157e+308"),
            (5e-324, "5e-324"),
            (333333333.33333329, "333333333.3333333"),
        ];
        for (input, expected) in cases {
            assert_eq!(
                es_number(*input).unwrap(),
                *expected,
                "ECMAScript rendering of {input:?}"
            );
        }
        assert!(es_number(f64::NAN).is_err());
        assert!(es_number(f64::INFINITY).is_err());
    }

    #[test]
    fn unit_tst_core_003_float_and_integer_agree_where_they_should() {
        // 1 and 1.0 are the same JSON number and must digest identically.
        let as_int = json!({"v": 1});
        let as_float: serde_json::Value = serde_json::from_str(r#"{"v":1.0}"#).unwrap();
        assert_eq!(canonicalize(&as_int).unwrap(), r#"{"v":1}"#);
        assert_eq!(canonicalize(&as_float).unwrap(), r#"{"v":1}"#);
        assert_eq!(
            Digest::of_json(&as_int).unwrap(),
            Digest::of_json(&as_float).unwrap()
        );
    }

    #[test]
    fn unit_string_escaping_is_minimal_and_correct() {
        let value = json!({"s": "a\"b\\c\nd\te\u{0}f\u{1f}g\u{8}h\u{c}i\rj"});
        assert_eq!(
            canonicalize(&value).unwrap(),
            r#"{"s":"a\"b\\c\nd\te\u0000f\u001fg\bh\fi\rj"}"#
        );
        // Non-ASCII is emitted literally as UTF-8, not escaped.
        let unicode = json!({"s": "é☃𝄞"});
        assert_eq!(canonicalize(&unicode).unwrap(), "{\"s\":\"é☃𝄞\"}");
    }

    #[test]
    fn unit_keys_sort_by_utf16_code_units() {
        // U+FF3A (fullwidth Z) is a single UTF-16 unit 0xFF3A; U+1D11E (musical
        // G clef) is the surrogate pair 0xD834 0xDD1E, which sorts BEFORE it in
        // UTF-16 order but AFTER it in UTF-8 byte order. This is the case that
        // distinguishes a correct implementation from a plausible one.
        let value: serde_json::Value =
            serde_json::from_str("{\"\u{1D11E}\":1,\"\u{FF3A}\":2}").unwrap();
        let canonical = canonicalize(&value).unwrap();
        let clef = canonical.find('\u{1D11E}').unwrap();
        let fullwidth = canonical.find('\u{FF3A}').unwrap();
        assert!(
            clef < fullwidth,
            "UTF-16 order puts the surrogate pair first; got {canonical:?}"
        );
    }

    #[test]
    fn unit_digest_of_document_omits_the_digest_field() {
        let with = json!({"a": 1, "digest": "sha256:deadbeef"});
        let without = json!({"a": 1});
        assert_eq!(
            Digest::of_document(&with, "digest").unwrap(),
            Digest::of_json(&without).unwrap()
        );
    }

    #[test]
    fn unit_digest_parsing() {
        let d = Digest::of_bytes(b"");
        // Known SHA-256 of the empty string.
        assert_eq!(
            d.as_str(),
            "sha256:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert!(Digest::parse(d.as_str()).is_ok());
        assert!(Digest::parse("e3b0c442").is_err());
        assert!(Digest::parse("sha256:XYZ").is_err());
        assert!(
            Digest::parse(
                "sha256:E3B0C44298FC1C149AFBF4C8996FB92427AE41E4649B934CA495991B7852B855"
            )
            .is_err(),
            "upper-case hex must be rejected: it digests differently as a string"
        );
    }

    #[test]
    fn unit_oversized_integers_are_refused() {
        let value = json!({"v": 9_007_199_254_740_993u64});
        assert_eq!(
            canonicalize(&value).unwrap_err().code,
            code::SEQ_NOT_CANONICAL
        );
    }
}
