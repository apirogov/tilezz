//! Sequences of common tiles.

use super::angles::upscale_angles;
use super::snake::Snake;
use crate::cyclotomic::{HasZZ10, HasZZ12, HasZZ4, HasZZ6};

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
// The needed edge markings (for proper Penrose tilings) are geometrically encoded as bump/dent pairs.

/// Return the sequence of the Penrose P3 "narrow" (thin) rhomb tile
/// over a compatible 10-fold cyclotomic ring.
pub fn penrose_p3_narrow<T: HasZZ10>() -> Snake<T> {
    // NOTE: see https://github.com/apirogov/tilezz/issues/29
    let seq: &[i8] = &[
        /* */ 1, 0, -1, 2, -1, 0, 0, // 1+ = N3  [0,7]
        /* */ 4, 0, 0, -1, 2, -1, 0, // 2+ = N11 [7,14]
        /* */ 1, 0, 1, -2, 1, 0, 0, //  2- = N17 [14,21]
        /* */ 4, 0, 0, 1, -2, 1, 0, //  1- = N25 [21,28]
    ];
    Snake::try_from(upscale_angles::<T>(10, seq).as_slice()).unwrap()
}

/// Return the sequence of the Penrose P3 "wide" (fat) rhomb tile
/// over a compatible 10-fold cyclotomic ring.
pub fn penrose_p3_wide<T: HasZZ10>() -> Snake<T> {
    // NOTE: see https://github.com/apirogov/tilezz/issues/29
    let seq: &[i8] = &[
        /* */ 2, 0, -1, 2, -1, 0, 0, // 1+ = W3   [0,7]
        /* */ 3, 0, 0, 1, -2, 1, 0, //  1- = W11  [7,14]
        /* */ 2, 0, 0, -1, 2, -1, 0, // 2+ = W18  [14,21]
        /* */ 3, 0, 1, -2, 1, 0, 0, //  2- = W24  [21,28]
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
        // Test the decorated Penrose rhombi (create a circular patch with a star inside).
        let fat_rat: Rat<ZZ10> = Rat::try_from(&penrose_p3_wide()).unwrap();
        let thin_rat: Rat<ZZ10> = Rat::try_from(&penrose_p3_narrow()).unwrap();

        let mut penta_star = fat_rat.glue((18, 24), &fat_rat);
        penta_star = penta_star.glue((3, 39), &penta_star);
        penta_star = penta_star.glue((32, 24), &fat_rat);

        let mut penta_circ = penta_star.clone();
        penta_circ = penta_circ.glue((3, 25), &thin_rat);
        penta_circ = penta_circ.glue((10, 25), &thin_rat);
        penta_circ = penta_circ.glue((10, 25), &thin_rat);
        penta_circ = penta_circ.glue((10, 25), &thin_rat);
        penta_circ = penta_circ.glue((10, 25), &thin_rat);
        assert_eq!(penta_circ.len(), 70);
    }
}
