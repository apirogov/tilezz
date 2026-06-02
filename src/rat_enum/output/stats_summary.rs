//! Per-boundary-length statistics over a finished enumeration:
//! total, achiral, rotational-symmetry histogram. Driven by the
//! `--stats` CLI flag.

use crate::geom::rat::lex_min_rot;
use crate::stringmatch::repetition_factor;

/// True iff this canonical sequence equals the canonical form of its
/// mirror image. The mirror polygon's CCW angle sequence is just
/// `reverse(seq)` -- traversing the mirrored vertices in CCW order
/// (in original coordinates) hits the original vertices in reverse
/// order, and each turn keeps its sign because the mirror's
/// orientation flip and the CCW-vs-CW flip cancel. Achiral polygons
/// survive a chirality quotient as a single class.
pub fn is_achiral(canonical: &[i8]) -> bool {
    let mut reflected: Vec<i8> = canonical.iter().rev().copied().collect();
    let offset = lex_min_rot(&reflected);
    reflected.rotate_left(offset);
    reflected == canonical
}

/// Print per-boundary-length statistics over the enumerated canonical
/// sequences: total count, achiral count, and a histogram of cyclic
/// rotational symmetry orders.
pub fn print_stats(rats: &[Vec<i8>]) {
    use std::collections::BTreeMap;
    // length -> (total, achiral, rep_factor -> count)
    let mut per_len: BTreeMap<usize, (usize, usize, BTreeMap<usize, usize>)> = BTreeMap::new();
    for seq in rats {
        let n = seq.len();
        let entry = per_len.entry(n).or_default();
        entry.0 += 1;
        if is_achiral(seq) {
            entry.1 += 1;
        }
        let rf = repetition_factor(seq);
        *entry.2.entry(rf).or_insert(0) += 1;
    }
    println!("statistics by boundary length:");
    println!("  length |  total | achiral (%) | rotational symmetry histogram (rep_factor: count)");
    println!("  -------+--------+-------------+--------------------------------------------------");
    for (n, (total, achiral, hist)) in &per_len {
        let pct = if *total > 0 {
            100.0 * (*achiral as f64) / (*total as f64)
        } else {
            0.0
        };
        let mut hist_str = String::new();
        for (rf, cnt) in hist {
            if !hist_str.is_empty() {
                hist_str.push_str("  ");
            }
            hist_str.push_str(&format!("x{rf}={cnt}"));
        }
        println!("  {n:>6} | {total:>6} | {achiral:>5} ({pct:>5.1}%) | {hist_str}");
    }
}
