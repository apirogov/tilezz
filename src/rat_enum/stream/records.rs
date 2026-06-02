//! Length-prefixed, bias-encoded binary record format for the
//! streaming pipeline.
//!
//! Each record is `1 + L` bytes:
//!   * byte 0: length L (0..=255).
//!   * bytes 1..=L: bias-encoded angles, `byte = (angle as i16 + 128) as u8`.
//!
//! Why the bias? The DAFSA expects rats in `(length asc, integer-lex
//! asc)` order. Lex-comparing raw signed-byte representations would
//! interleave negative and positive angles wrong (`-1` as `0xFF`
//! sorts above `+1` as `0x01`). Adding 128 maps `i8::MIN` to `0` and
//! `i8::MAX` to `0xFF`, so byte-lex compares **agree** with integer
//! lex. Combined with the leading length byte, byte-lex order on
//! whole records is `(length, lex)` -- exactly the order Stage 3
//! consumes.

/// Offset added to each `i8` angle when encoding so byte-lex
/// comparison agrees with integer-lex comparison.
pub const ANGLE_BIAS: i16 = 128;

/// Append the bias-encoded record for `canonical` to `out`. Appends
/// `1 + canonical.len()` bytes total.
#[inline]
pub fn encode_record(canonical: &[i8], out: &mut Vec<u8>) {
    debug_assert!(
        canonical.len() <= u8::MAX as usize,
        "record longer than 255"
    );
    out.push(canonical.len() as u8);
    for &a in canonical {
        out.push((a as i16 + ANGLE_BIAS) as u8);
    }
}

/// Decode the first record at the start of `bytes`. Returns
/// `(remaining, decoded_angles)` or `None` if `bytes` is shorter
/// than the declared length.
#[inline]
pub fn decode_record(bytes: &[u8]) -> Option<(&[u8], Vec<i8>)> {
    if bytes.is_empty() {
        return None;
    }
    let len = bytes[0] as usize;
    if bytes.len() < 1 + len {
        return None;
    }
    let angles: Vec<i8> = bytes[1..=len]
        .iter()
        .map(|&b| (b as i16 - ANGLE_BIAS) as i8)
        .collect();
    Some((&bytes[1 + len..], angles))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Round-trip: encode then decode recovers the input.
    #[test]
    fn round_trip() {
        let cases: &[&[i8]] = &[&[], &[0], &[-5, 3, -2, 4], &[127, -128, 0, 1, -1]];
        for &seq in cases {
            let mut buf = Vec::new();
            encode_record(seq, &mut buf);
            let (rest, decoded) = decode_record(&buf).expect("decoded");
            assert!(rest.is_empty(), "trailing bytes");
            assert_eq!(decoded, seq);
        }
    }

    /// Byte-lex order of encoded records equals (length asc, integer
    /// lex asc) of the underlying sequences. The DAFSA depends on
    /// this property; it's the whole point of the bias.
    #[test]
    fn byte_lex_equals_length_then_int_lex() {
        let mut seqs: Vec<Vec<i8>> = vec![
            vec![1, 2, 3],
            vec![-1, 0],
            vec![1, 1],
            vec![],
            vec![1, 2],
            vec![-5, 0, 3],
            vec![-5, 0, 2],
            vec![1, -1, 0],
        ];
        // Reference: sort by (len, lex).
        seqs.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));

        // Test: encode all, sort by byte-lex, decode, check matches.
        let mut encoded: Vec<Vec<u8>> = seqs
            .iter()
            .map(|s| {
                let mut buf = Vec::new();
                encode_record(s, &mut buf);
                buf
            })
            .collect();
        encoded.sort();

        let decoded: Vec<Vec<i8>> = encoded
            .iter()
            .map(|e| decode_record(e).expect("decode").1)
            .collect();
        assert_eq!(decoded, seqs);
    }
}
