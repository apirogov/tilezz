//! Sequences of common tiles.

use super::angles::upscale_angles;
use super::snake::Snake;
use crate::cyclotomic::{HasZZ4, HasZZ6, HasZZ10, HasZZ12};

/// Return sequence of a square tile over a compatible ring (divisible by 4).
pub fn square<T: HasZZ4>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(4, &[1, 1, 1, 1]).as_slice()).unwrap()
}

/// Return sequence of a equilateral triangle tile over a compatible ring (divisible by 6).
pub fn triangle<T: HasZZ6>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(6, &[2, 2, 2]).as_slice()).unwrap()
}

/// Return sequence of a hexagon tile over a compatible ring (divisible by 6).
pub fn hexagon<T: HasZZ6>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(6, &[1, 1, 1, 1, 1, 1]).as_slice()).unwrap()
}

// ----

/// Return sequence of the spectre tile over a compatible ring (divisible by 12).
pub fn spectre<T: HasZZ12>() -> Snake<T> {
    Snake::try_from(
        upscale_angles::<T>(12, &[3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2]).as_slice(),
    )
    .unwrap()
}

// ----

// NOTE: Penrose P2 (kite and dart) cannot be represented with only unit steps,
// but could with an extension of Snakes to take real-valued positive elements of the cyclotomic ring
// (see https://github.com/apirogov/tilezz/issues/1)
// However, P3 (pair of rhombs) have edges of the same length, so it works.
//
// The matching rules (see https://commons.wikimedia.org/wiki/File:Penrose_rhombs_matching_rules.svg)
// are baked into the edge geometry as bump/dent gadgets. Crucially there are
// TWO distinct gadget shapes -- the analog of Penrose's single- vs double-arrow
// markings -- so that the two arrow types can never be confused (see
// https://github.com/apirogov/tilezz/issues/29):
//   * type 1: a shallow +-1 spike  ([-1,2,-1] bump / [1,-2,1] dent), a 4-step edge;
//   * type 2: a taller  +-2 spike  ([-2,4,-2] bump / [2,-4,2] dent), a 5-step edge.
// A bump (+) mates ONLY the dent (-) of the SAME type (their step runs are exact
// reverse-complements), so a type-1 edge can never glue to a type-2 edge. The two
// edge types share one length because 2cos36 = 1 + 2cos72, so the +-2 edge just
// needs one extra straight step; this keeps each tile an equilateral rhombus.

/// Return the sequence of the Penrose P3 "narrow" (thin) rhomb tile
/// over a compatible 10-fold cyclotomic ring.
///
/// Four unit-length sides; CCW edge markings (sign = bump(+)/dent(-),
/// magnitude = arrow type): -1, +1, -2, +2. Enumeration starts at an acute
/// (36deg) vertex, so the corner turns are 4,1,4,1. The +-2 (type-2) spikes are
/// buffered away from the sharp 36deg corners by the lone extra straight step
/// so the outline stays simple.
pub fn penrose_p3_narrow<T: HasZZ10>() -> Snake<T> {
    let seq: &[i8] = &[
        4, 1, -2, 1, //     -1  type-1 dent  (acute corner 4)
        1, -1, 2, -1, //    +1  type-1 bump  (obtuse corner 1)
        4, 0, 2, -4, 2, //  -2  type-2 dent  (acute corner 4, pad before spike)
        1, -2, 4, -2, 0, // +2  type-2 bump  (obtuse corner 1, pad after spike)
    ];
    Snake::try_from(upscale_angles::<T>(10, seq).as_slice()).unwrap()
}

/// Return the sequence of the Penrose P3 "wide" (fat) rhomb tile
/// over a compatible 10-fold cyclotomic ring.
///
/// Four unit-length sides; CCW edge markings: -2, -1, +1, +2. Enumeration
/// starts at an acute (72deg) vertex, so the corner turns are 3,2,3,2. The +-2
/// (type-2) spikes are buffered away from the acute corners.
pub fn penrose_p3_wide<T: HasZZ10>() -> Snake<T> {
    let seq: &[i8] = &[
        3, 0, 2, -4, 2, //  -2  type-2 dent  (acute corner 3, pad before spike)
        2, 1, -2, 1, //     -1  type-1 dent  (obtuse corner 2)
        3, -1, 2, -1, //    +1  type-1 bump  (acute corner 3)
        2, -2, 4, -2, 0, // +2  type-2 bump  (obtuse corner 2, pad after spike)
    ];
    Snake::try_from(upscale_angles::<T>(10, seq).as_slice()).unwrap()
}

// ---- Tetrominoes ----
//
// The seven standard 4-cell polyomino shapes over a 4-fold cyclotomic
// ring. Each is a connected union of unit squares with the canonical
// letter-shaped outline.

/// The square / O-tetromino: a 2×2 block of unit cells.
#[allow(non_snake_case)]
pub fn tetromino_O<T: HasZZ4>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(4, &[0, 1, 0, 1, 0, 1, 0, 1]).as_slice()).unwrap()
}

/// The straight / I-tetromino: a 1×4 strip.
#[allow(non_snake_case)]
pub fn tetromino_I<T: HasZZ4>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(4, &[0, 0, 0, 1, 1, 0, 0, 0, 1, 1]).as_slice()).unwrap()
}

/// The T-tetromino: three cells in a row with one centered above.
#[allow(non_snake_case)]
pub fn tetromino_T<T: HasZZ4>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(4, &[-1, 1, 1, -1, 1, 1, 0, 0, 1, 1]).as_slice()).unwrap()
}

/// The S-tetromino (right-handed skew): like an S laid on its side.
#[allow(non_snake_case)]
pub fn tetromino_S<T: HasZZ4>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(4, &[-1, 1, 1, 0, 1, -1, 1, 1, 0, 1]).as_slice()).unwrap()
}

/// The Z-tetromino: the mirror image of S.
#[allow(non_snake_case)]
pub fn tetromino_Z<T: HasZZ4>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(4, &[-1, 1, 0, 1, 1, -1, 1, 0, 1, 1]).as_slice()).unwrap()
}

/// The J-tetromino: three cells in a row with one attached to the
/// upper-left of the leftmost cell.
#[allow(non_snake_case)]
pub fn tetromino_J<T: HasZZ4>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(4, &[-1, 1, 1, 0, 1, 0, 0, 1, 1, 0]).as_slice()).unwrap()
}

/// The L-tetromino: the mirror image of J.
#[allow(non_snake_case)]
pub fn tetromino_L<T: HasZZ4>() -> Snake<T> {
    Snake::try_from(upscale_angles::<T>(4, &[-1, 0, 1, 1, 0, 0, 1, 0, 1, 1]).as_slice()).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ10, ZZ24};
    use crate::geom::rat::Rat;

    #[test]
    fn test_upscaling() {
        // test adaptive scaling dependent on chosen ring
        let seq_zz24: Vec<i8> = vec![6, 4, 0, 4, -6, 4, 6, 4, -6, 4, 6, -4, 6, -4];
        let s: Snake<ZZ24> = spectre();
        assert_eq!(s.angles(), seq_zz24);
    }

    #[test]
    fn test_penrose() {
        // The decorated P3 rhombi: 4 unit-edge sides, each a meta-edge carrying
        // ONE of two distinct, non-interchangeable gadgets (the +-1 "single
        // arrow" / +-2 "double arrow" spikes). A bump mates only the dent of the
        // same type.
        let fat: Rat<ZZ10> = Rat::try_from(&penrose_p3_wide()).unwrap();
        let thin: Rat<ZZ10> = Rat::try_from(&penrose_p3_narrow()).unwrap();

        // perimeter = two type-1 edges (4 steps) + two type-2 edges (5 steps).
        assert_eq!(fat.len(), 18);
        assert_eq!(thin.len(), 18);

        // The CCW edge markings must match the Penrose P3 matching rules:
        //   thin : -1, +1, -2, +2        fat : -2, -1, +1, +2
        // (sign = bump(+)/dent(-), magnitude = arrow type). Read off the raw
        // (pre-canonicalization) angle sequence at the known meta-edge offsets:
        // type = max|interior turn|/2; bump iff the first nonzero interior turn
        // is negative (bump = [-t,2t,-t], dent = [t,-2t,t]).
        fn marking(seq: &[i8], starts: &[usize]) -> Vec<(char, i8)> {
            let n = seq.len();
            (0..starts.len())
                .map(|k| {
                    let interior = &seq[starts[k] + 1..*starts.get(k + 1).unwrap_or(&n)];
                    let ty = interior.iter().map(|t| t.abs()).max().unwrap() / 2;
                    let first = *interior.iter().find(|&&t| t != 0).unwrap();
                    (if first < 0 { '+' } else { '-' }, ty)
                })
                .collect()
        }
        assert_eq!(
            marking(penrose_p3_narrow::<ZZ10>().angles(), &[0, 4, 8, 13]),
            [('-', 1), ('+', 1), ('-', 2), ('+', 2)]
        );
        assert_eq!(
            marking(penrose_p3_wide::<ZZ10>().angles(), &[0, 5, 9, 13]),
            [('-', 2), ('-', 1), ('+', 1), ('+', 2)]
        );

        // Two distinct gadget shapes are present (a single-gadget encoding would
        // lack the +-2 spike's magnitude-4 turn or the +-1 [-1,2,-1] run).
        for t in [&fat, &thin] {
            assert!(t.seq().iter().any(|&a| a.abs() == 4));
            assert!(t.seq().windows(3).any(|w| matches!(w, [-1, 2, -1] | [1, -2, 1])));
        }

        // The decorated rhombi still assemble into a valid (simple, hole-free) patch.
        assert!(fat.try_glue((0, 0), &thin).is_ok());
    }
}
