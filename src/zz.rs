use super::traits::{Ccw, InnerIntType, IntRing};
use super::zzbase::{
    zz_partial_signum_1_sym, zz_partial_signum_2_sym, zz_partial_signum_4_sym,
    zz_partial_signum_fallback, Frac, GInt, ZZBase, ZZNum, ZZParams,
};
use crate::traits::ComplexIntRing;
use crate::{zz_base_impl, zz_ops_impl};

use num_traits::{One, Zero};
use std::f64::consts::SQRT_2;
use std::fmt;
use std::fmt::Display;
use std::marker::PhantomData;
use std::ops::{Add, Mul, Neg, Sub};

// definitions needed to derive different ZZn types
//
// NOTE: The bases here were derived in an ad-hoc manner. Apparently, a systematic approach
// for deriving integral bases is described in this paper:
// https://link.springer.com/article/10.1007/s002000050065

// numeric constants for the chosen (in general non-unique) linearly independent algebraic bases
// --------
// NOTE: sqrt(2 + sqrt(2)) = (sqrt(2)+1)sqrt(2 - sqrt(2))
// and   sqrt(2 - sqrt(2)) = (sqrt(2)-1)sqrt(2 + sqrt(2))
const ZZ16_Y: f64 = 2.0 + SQRT_2;

// NOTE: sqrt(2+sqrt(2+sqrt(2))) = ( 1 + sqrt(2) + sqrt(2(2+sqrt(2)))) sqrt(2-sqrt(2+sqrt(2)))
// and   sqrt(2-sqrt(2+sqrt(2))) = (-1 - sqrt(2) + sqrt(2(2+sqrt(2)))) sqrt(2+sqrt(2+sqrt(2)))
const ZZ32_Z: f64 = 2.0 + 1.84775906502257351225; // 2+sqrt(2+sqrt(2))

// NOTE: Let x = sqrt(5), y = sqrt(2*(5-sqrt(5))), y' = sqrt(2*(5+sqrt(5)))
// We have: y = 1/2(xy' - y') ^ y' = 1/2(xy + y)
const SQRT_5: f64 = 2.23606797749978969;
const ZZ10_Y: f64 = 2.0 * (5.0 - SQRT_5);

// NOTE: could construct ZZ120 and ZZ240 as well, using only the roots we already use:
//
// e^(i*pi/60)  = 1/16 (1+i) (2 sqrt(5+sqrt(5)) - i sqrt(2) - sqrt(2)sqrt(3)
//              + i sqrt(2)sqrt(5) + sqrt(2)sqrt(3)sqrt(5) - 2i sqrt(3)sqrt(5+sqrt(5)))

// e^(i*pi/120) = 1/16 (sqrt(2-sqrt(2)) (sqrt(3) + sqrt(15) - sqrt(2(5-sqrt(5)))) + sqrt(2+sqrt(2)) (1 + sqrt(5) + sqrt(6(5-sqrt(5)))))
//              + 1/16i(sqrt(2+sqrt(2)) (sqrt(3) + sqrt(15) - sqrt(2(5-sqrt(5)))) - sqrt(2-sqrt(2)) (1 + sqrt(5) + sqrt(6(5-sqrt(5)))))

// --------

/// NOTE:
///
/// Each ring is represented by a set of linearly independent "square root units"
/// (one of these is units is always 1, i.e. the "pure integer unit")
/// and all their symbolic products, i.e. for k distinct units there are 2^k entries.
/// A value is represented as a linear combination of these units,
/// where the scalar coefficients multiplied by Gaussian integer coefficients,
/// all of that divided by a fixed common ring-specific scaling factor.
///
/// Addition is simply component-wise addition of coefficients, multiplication however
/// is more involved and can be represented as a 3D matrix (a rank 3 tensor?) that describes
/// the expanded coefficients for each component of the outer product of the two input vectors.
///
/// For each pair of input vectors u (index i) and v (index j) and the multiplication tensor M
/// that assigns to each i,j,k the contribution along the dimension k of the value u_i * v_j,
/// the result of the multiplication w (index k) is given by w_k := sum_{i,j}(u_i * v_j * M_{i,j,k}).
///
/// This conceptual 3D symbolic vector multiplication matrix is not implemented
/// explicitly, instead the resulting coefficients are manually collected and
/// simplified for each ring separately.

/// Gauss integers
pub const ZZ4_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 4,
    sym_roots_num: 1,
    sym_roots_sqs: &[1.],
    sym_roots_lbls: &["1"],
    scaling_fac: 1,
    ccw_unit_coeffs: &[[0, 1]],
};
fn zz4_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 1), *array_ref!(y, 0, 1)] {
        [[a], [b]] => vec![a * b],
    }
}
/// Eisenstein integers
pub const ZZ6_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 6,
    sym_roots_num: 2,
    sym_roots_sqs: &[1., 3.],
    sym_roots_lbls: &["1", "3"],
    scaling_fac: 2,
    ccw_unit_coeffs: &[[1, 0], [0, 1]],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
///    c d
/// a [1 s]
/// b [s 3]
/// where s = sqrt(3)
fn zz6_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 2), *array_ref!(y, 0, 2)] {
        [[a, b], [c, d]] => vec![a * c + (b * d * 3), a * d + b * c],
    }
}
/// Compass integers
pub const ZZ8_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 8,
    sym_roots_num: 2,
    sym_roots_sqs: &[1., 2.],
    sym_roots_lbls: &["1", "2"],
    scaling_fac: 2,
    ccw_unit_coeffs: &[[0, 0], [1, 1]],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
///    c d
/// a [1 s]
/// b [s 2]
/// where s = sqrt(2)
fn zz8_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 2), *array_ref!(y, 0, 2)] {
        [[a, b], [c, d]] => vec![a * c + (b * d * 2), a * d + b * c],
    }
}
/// Halfrose integers (impractical, has no quarter turn)
pub const ZZ10_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 10,
    sym_roots_num: 4,
    sym_roots_sqs: &[1., 5., ZZ10_Y, 5. * ZZ10_Y],
    sym_roots_lbls: &["1", "5", "2(5-sqrt(5))", "10(5-sqrt(5))"],
    scaling_fac: 8,
    ccw_unit_coeffs: &[[2, 0], [2, 0], [0, 2], [0, 0]],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
///    e     f      g        h
/// a [1   ,  x   ,   y    ,   xy   ]
/// b [ x  , 5    ,  xy    ,  5 y   ]
/// c [  y ,  xy  , 10-2x  , 10(x-1)]
/// d [ xy , 5 y  , 10(x-1), 10(5-x)]
/// where x = sqrt(5), y = sqrt(2(5-x))
fn zz10_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 4), *array_ref!(y, 0, 4)] {
        [[a, b, c, d], [e, f, g, h]] => {
            let c1 = a * e + b * f * 5 + (c * g + d * h * 5 - c * h - d * g) * 10;
            let c2 = a * f + b * e - c * g * 2 + (c * h + d * g - d * h) * 10;
            let c3 = a * g + (b * h + d * f) * 5 + c * e;
            let c4 = a * h + b * g + c * f + d * e;
            vec![c1, c2, c3, c4]
        }
    }
}
/// Clock integers
pub const ZZ12_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 12,
    sym_roots_num: ZZ6_PARAMS.sym_roots_num,
    sym_roots_sqs: ZZ6_PARAMS.sym_roots_sqs,
    sym_roots_lbls: ZZ6_PARAMS.sym_roots_lbls,
    scaling_fac: ZZ6_PARAMS.scaling_fac,
    ccw_unit_coeffs: &[[0, 1], [1, 0]],
};
fn zz12_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    return zz6_mul(x, y);
}
/// Hex integers
pub const ZZ16_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 16,
    sym_roots_num: 4,
    sym_roots_sqs: &[1.0, 2.0, ZZ16_Y, 2.0 * ZZ16_Y],
    sym_roots_lbls: &["1", "2", "2+sqrt(2)", "2(2+sqrt(2))"],
    scaling_fac: 2,
    ccw_unit_coeffs: &[[0, 0], [0, 0], [1, -1], [0, 1]],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
///     e    f     g      h
/// a [1   ,  x  ,   y  ,  xy  ]
/// b [ x  , 2   ,  xy  , 2 y  ]
/// c [  y ,  xy , 2+x  , 2+2x ]
/// d [ xy , 2 y , 2+2x , 4+2x ]
/// where x = sqrt(2), y = sqrt(2+sqrt(2))
fn zz16_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 4), *array_ref!(y, 0, 4)] {
        [[a, b, c, d], [e, f, g, h]] => {
            let c1 = a * e + (b * f + c * g + c * h + d * g + d * h * 2) * 2;
            let c2 = a * f + b * e + c * g + (c * h + d * g + d * h) * 2;
            let c3 = a * g + c * e + (b * h + d * f) * 2;
            let c4 = a * h + b * g + c * f + d * e;
            vec![c1, c2, c3, c4]
        }
    }
}
/// Penrose integers
pub const ZZ20_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 20,
    sym_roots_num: ZZ10_PARAMS.sym_roots_num,
    sym_roots_sqs: ZZ10_PARAMS.sym_roots_sqs,
    sym_roots_lbls: ZZ10_PARAMS.sym_roots_lbls,
    scaling_fac: ZZ10_PARAMS.scaling_fac,
    ccw_unit_coeffs: &[[0, -2], [0, 2], [1, 0], [1, 0]],
};
fn zz20_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    zz10_mul(x, y)
}
/// Digiclock integers
pub const ZZ24_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 24,
    sym_roots_num: 4,
    sym_roots_sqs: &[1.0, 2.0, 3.0, 6.0],
    sym_roots_lbls: &["1", "2", "3", "6"],
    scaling_fac: 4,
    ccw_unit_coeffs: &[[0, 0], [1, -1], [0, 0], [1, 1]],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
///     e    f     g    h
/// a [1   ,  x  ,   y ,  xy ]
/// b [ x  , 2   ,  xy , 2 y ]
/// c [  y ,  xy , 3   , 3x  ]
/// d [ xy , 2 y , 3x  , 6   ]
fn zz24_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 4), *array_ref!(y, 0, 4)] {
        [[a, b, c, d], [e, f, g, h]] => {
            let c1 = a * e + b * f * 2 + c * g * 3 + d * h * 6;
            let c2 = a * f + b * e + (c * h + d * g) * 3;
            let c3 = a * g + c * e + (b * h + d * f) * 2;
            let c4 = a * h + b * g + c * f + d * e;
            vec![c1, c2, c3, c4]
        }
    }
}
/// Month integers
pub const ZZ30_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 30,
    sym_roots_num: 8,
    sym_roots_sqs: &[
        1.,
        3.,
        5.,
        ZZ10_Y,
        15.,
        3. * ZZ10_Y,
        5. * ZZ10_Y,
        15. * ZZ10_Y,
    ],
    sym_roots_lbls: &[
        "1",
        "3",
        "5",
        "2(5-sqrt(5))",
        "15",
        "6(5-sqrt(5))",
        "10(5-sqrt(5))",
        "30(5-sqrt(5))",
    ],
    // e^(i*pi/15) = 1/16  (-2 + 2 sqrt(5) + sqrt(6 (5 - sqrt(5))) + sqrt(30 (5 - sqrt(5))))
    //             + 1/16 i(2 sqrt(3) - 2 sqrt(15) + sqrt(2 (5 - sqrt(5))) + sqrt(10 (5 - sqrt(5))))
    scaling_fac: 16,
    ccw_unit_coeffs: &[
        [-2, 0],
        [0, 2],
        [2, 0],
        [0, 1],
        [0, -2],
        [1, 0],
        [0, 1],
        [1, 0],
    ],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
/// [1    ,  x   ,   y    ,    z     ,  xy  ,  x z     ,   yz     ,  xyz      ]
/// [ x   , 3    ,
/// [  y  ,  xy  ,  5     ,
/// [   z ,  x z ,    yz  , 10-2y    ,
/// [ xy  , 3 y  ,  5x    , xyz      , 15   ,
/// [ x z , 3  z ,   xyz  , 10x-2xy  ,  3yz , 3(10-2y) ,
/// [  yz ,  xyz ,  5  z  , 10y-10   ,  5xz , 10xy-10x , 5(10-2y)
/// [ xyz , 3 yz ,  5x z  , 10xy-10x ,  15z , 30y-30   , 50x-10xy , 15(10-2y)
/// where x = sqrt(3), y = sqrt(5), z = sqrt(2(5-y)) = sqrt(10-2y)
fn zz30_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 8), *array_ref!(y, 0, 8)] {
        [[l1, l2, l3, l4, l5, l6, l7, l8], [r1, r2, r3, r4, r5, r6, r7, r8]] => {
            // 1
            let c1_d = ((l1 * r1) + (l2 * r2) * 3 + (l3 * r3) * 5 + (l4 * r4) * 10)
                + ((l5 * r5) * 15 + (l6 * r6) * 30 + (l7 * r7) * 50 + (l8 * r8) * 150);
            let c1_a = l4 * r7 * (-10) + l6 * r8 * (-30);
            let c1_b = r4 * l7 * (-10) + r6 * l8 * (-30);
            let c1 = c1_d + c1_a + c1_b;

            // x
            let c2_a = l1 * r2 + l3 * r5 * 5 + (l4 * r6 - l4 * r8 - l6 * r7) * 10 + l7 * r8 * 50;
            let c2_b = r1 * l2 + r3 * l5 * 5 + (r4 * l6 - r4 * l8 - r6 * l7) * 10 + r7 * l8 * 50;
            let c2 = c2_a + c2_b;

            // y
            let c3_d = r4 * l4 * (-2) + r6 * l6 * (-6) + r7 * l7 * (-10) + r8 * l8 * (-30);
            let c3_a = l1 * r3 + l2 * r5 * 3 + l4 * r7 * 10 + l6 * r8 * 30;
            let c3_b = r1 * l3 + r2 * l5 * 3 + r4 * l7 * 10 + r6 * l8 * 30;
            let c3 = c3_d + c3_a + c3_b;

            // z
            let c4_a = l1 * r4 + l2 * r6 * 3 + l3 * r7 * 5 + l5 * r8 * 15;
            let c4_b = r1 * l4 + r2 * l6 * 3 + r3 * l7 * 5 + r5 * l8 * 15;
            let c4 = c4_a + c4_b;

            // xy
            let c5_a = l1 * r5 + l2 * r3 - l4 * r6 * 2 + (l4 * r8 + l6 * r7 - l7 * r8) * 10;
            let c5_b = r1 * l5 + r2 * l3 - r4 * l6 * 2 + (r4 * l8 + r6 * l7 - r7 * l8) * 10;
            let c5 = c5_a + c5_b;

            // xz
            let c6_a = l1 * r6 + l2 * r4 + (l3 * r8 + l5 * r7) * 5;
            let c6_b = r1 * l6 + r2 * l4 + (r3 * l8 + r5 * l7) * 5;
            let c6 = c6_a + c6_b;

            // yz
            let c7_a = l1 * r7 + l2 * r8 * 3 + l3 * r4 + l5 * r6 * 3;
            let c7_b = r1 * l7 + r2 * l8 * 3 + r3 * l4 + r5 * l6 * 3;
            let c7 = c7_a + c7_b;

            // xyz
            let c8 = r1 * l8 + r2 * l7 + r3 * l6 + r4 * l5 + r5 * l4 + r6 * l3 + r7 * l2 + r8 * l1;

            vec![c1, c2, c3, c4, c5, c6, c7, c8]
        }
    }
}
/// Hex integers
pub const ZZ32_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 32,
    sym_roots_num: 8,
    sym_roots_sqs: &[
        1.,
        2.,
        ZZ16_Y,
        ZZ32_Z,
        2. * ZZ16_Y,
        2. * ZZ32_Z,
        ZZ16_Y * ZZ32_Z,
        2. * ZZ16_Y * ZZ32_Z,
    ],
    sym_roots_lbls: &[
        "1",
        "2",
        "2+sqrt(2)",
        "2+sqrt(2+sqrt(2))",
        "2(2+sqrt(2))",
        "2(2+sqrt(2+sqrt(2)))",
        "(2+sqrt(2))(2+sqrt(2+sqrt(2)))",
        "2(2+sqrt(2))(2+sqrt(2+sqrt(2)))",
    ],
    // e^(i*pi/16) = 1/2((1-i)*sqrt(2+sqrt(2+sqrt(2)))
    //             - i*sqrt(2)sqrt(2+sqrt(2+sqrt(2)))
    //             + i*sqrt(2)sqrt(2+sqrt(2))sqrt(2+sqrt(2+sqrt(2))))
    scaling_fac: 2,
    ccw_unit_coeffs: &[
        [0, 0],
        [0, 0],
        [0, 0],
        [1, -1],
        [0, 0],
        [0, -1],
        [0, 0],
        [0, 1],
    ],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
/// [1    ,  x   ,   y  ,    z   ,  xy  ,  x z    ,   yz        ,  xyz ]
/// [ x   , 2    ,
/// [  y  ,  xy  , 2+x  ,
/// [   z ,  x z ,   yz , 2+y    ,
/// [ xy  , 2 y  , 2+2x ,  xyz   , 4+2x ,
/// [ x z , 2 z  ,  xyz ,2x+xy   , 2 yz , 4+2y    ,
/// [  yz ,  xyz ,2z+ xz,2+ x+2y ,2z+2xz, 2+2x+2xy, 4+2x+2y+ xy
/// [ xyz , 2 yz ,2z+2xz,2+2x+2xy,4z+2xz, 4+2x+4y , 4+4x+2y+2xy , 8+4x+4y+2xy
/// where x = sqrt(2), y = sqrt(2+x), z = sqrt(2+y)
fn zz32_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 8), *array_ref!(y, 0, 8)] {
        [[l1, l2, l3, l4, l5, l6, l7, l8], [r1, r2, r3, r4, r5, r6, r7, r8]] => {
            // 1
            let c1_d = (l1 * r1)
                + (l2 * r2 + l3 * r3 + l4 * r4) * 2
                + (l5 * r5 + l6 * r6 + l7 * r7) * 4
                + l8 * r8 * 8;
            let c1_a = (l3 * r5 + l4 * r7 + l4 * r8 + l6 * r7) * 2 + (l6 * r8 + l7 * r8) * 4;
            let c1_b = (r3 * l5 + r4 * l7 + r4 * l8 + r6 * l7) * 2 + (r6 * l8 + r7 * l8) * 4;
            let c1 = c1_d + c1_a + c1_b;

            // x
            let c2_d = r3 * l3 + (r5 * l5 + r7 * l7) * 2 + r8 * l8 * 4;
            let c2_a = (l1 * r2 + l4 * r7)
                + (l3 * r5 + l4 * r6 + l4 * r8 + l6 * r7 + l6 * r8) * 2
                + l7 * r8 * 4;
            let c2_b = (r1 * l2 + r4 * l7)
                + (r3 * l5 + r4 * l6 + r4 * l8 + r6 * l7 + r6 * l8) * 2
                + r7 * l8 * 4;
            let c2 = c2_d + c2_a + c2_b;

            // y
            let c3_d = l4 * r4 + (l6 * r6 + l7 * r7) * 2 + (l8 * r8) * 4;
            let c3_a = l1 * r3 + (l2 * r5 + l4 * r7 + l7 * r8) * 2 + (l6 * r8) * 4;
            let c3_b = r1 * l3 + (r2 * l5 + r4 * l7 + r7 * l8) * 2 + (r6 * l8) * 4;
            let c3 = c3_d + c3_a + c3_b;

            // z
            let c4_a = l1 * r4 + (l2 * r6 + l3 * r7 + l3 * r8 + l5 * r7) * 2 + l5 * r8 * 4;
            let c4_b = r1 * l4 + (r2 * l6 + r3 * l7 + r3 * l8 + r5 * l7) * 2 + r5 * l8 * 4;
            let c4 = c4_a + c4_b;

            // xy
            let c5_d = l7 * r7 + l8 * r8 * 2;
            let c5_a = l1 * r5 + l2 * r3 + l4 * r6 + (l4 * r8 + l6 * r7 + l7 * r8) * 2;
            let c5_b = r1 * l5 + r2 * l3 + r4 * l6 + (r4 * l8 + r6 * l7 + r7 * l8) * 2;
            let c5 = c5_d + c5_a + c5_b;

            // xz
            let c6_a = l1 * r6 + l2 * r4 + l3 * r7 + (l3 * r8 + l5 * r7 + l5 * r8) * 2;
            let c6_b = r1 * l6 + r2 * l4 + r3 * l7 + (r3 * l8 + r5 * l7 + r5 * l8) * 2;
            let c6 = c6_a + c6_b;

            // yz
            let c7_a = l1 * r7 + l3 * r4 + (l2 * r8 + l5 * r6) * 2;
            let c7_b = r1 * l7 + r3 * l4 + (r2 * l8 + r5 * l6) * 2;
            let c7 = c7_a + c7_b;

            // xyz
            let c8 = l1 * r8 + l2 * r7 + l3 * r6 + l4 * r5 + l5 * r4 + l6 * r3 + l7 * r2 + l8 * r1;

            vec![c1, c2, c3, c4, c5, c6, c7, c8]
        }
    }
}
pub const ZZ60_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 60,
    sym_roots_num: ZZ30_PARAMS.sym_roots_num,
    sym_roots_sqs: ZZ30_PARAMS.sym_roots_sqs,
    sym_roots_lbls: ZZ30_PARAMS.sym_roots_lbls,
    scaling_fac: ZZ30_PARAMS.scaling_fac,
    // e^(i*pi/30) = 1/8 (sqrt(3) + sqrt(15) + sqrt(2 (5 - sqrt(5))) + i(-1 - sqrt(5) + sqrt(3)sqrt(2 (5 - sqrt(5)))))
    ccw_unit_coeffs: &[
        [0, -2],
        [2, 0],
        [0, -2],
        [2, 0],
        [2, 0],
        [0, 2],
        [0, 0],
        [0, 0],
    ],
};
fn zz60_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    zz30_mul(x, y)
}
// --------

// generate boilerplate implementations
zz_base_impl!(ZZ4, ZZ4_PARAMS, zz4_mul, zz_partial_signum_1_sym);
zz_base_impl!(ZZ6, ZZ6_PARAMS, zz6_mul, zz_partial_signum_2_sym);
zz_base_impl!(ZZ8, ZZ8_PARAMS, zz8_mul, zz_partial_signum_2_sym);
zz_base_impl!(ZZ10, ZZ10_PARAMS, zz10_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ12, ZZ12_PARAMS, zz12_mul, zz_partial_signum_2_sym);
zz_base_impl!(ZZ16, ZZ16_PARAMS, zz16_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ20, ZZ20_PARAMS, zz20_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ24, ZZ24_PARAMS, zz24_mul, zz_partial_signum_4_sym);
zz_base_impl!(ZZ30, ZZ30_PARAMS, zz30_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ32, ZZ32_PARAMS, zz32_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ60, ZZ60_PARAMS, zz60_mul, zz_partial_signum_fallback);
zz_ops_impl!(ZZ4 ZZ6 ZZ8 ZZ10 ZZ12 ZZ16 ZZ20 ZZ24 ZZ30 ZZ32 ZZ60);

pub mod constants {
    use super::*;

    pub fn zz_units_sum<T: ZZNum>() -> T {
        let mut p = T::zero();
        for i in 1..T::turn() {
            p = p + T::unit(i).scale(i as i64);
        }
        p
    }

    // NOTE: as we can get the real-valued square roots represented,
    // it means that we can represent any linear combination
    // in a ring that supports quarter turn rotation.

    pub fn zz8_sqrt2() -> ZZ8 {
        ZZ8::unit(1) + ZZ8::unit(-1)
    }
    pub fn zz16_sqrt2() -> ZZ16 {
        ZZ16::unit(2) + ZZ16::unit(-2)
    }
    pub fn zz24_sqrt2() -> ZZ24 {
        ZZ24::unit(3) + ZZ24::unit(-3)
    }
    pub fn zz32_sqrt2() -> ZZ32 {
        ZZ32::unit(4) + ZZ32::unit(-4)
    }

    pub fn zz6_isqrt3() -> ZZ6 {
        ZZ6::unit(1) + ZZ6::unit(2)
    }
    pub fn zz12_sqrt3() -> ZZ12 {
        ZZ12::unit(1) + ZZ12::unit(-1)
    }
    pub fn zz24_sqrt3() -> ZZ24 {
        ZZ24::unit(2) + ZZ24::unit(-2)
    }

    pub fn zz10_isqrt_penta() -> ZZ10 {
        ZZ10::unit(1) * ZZ10::from(4) - ZZ10::one() - zz10_sqrt5()
    }
    pub fn zz20_half_sqrt_penta() -> ZZ20 {
        ZZ20::unit(3) + ZZ20::unit(-3)
    }

    pub fn zz10_sqrt5() -> ZZ10 {
        (ZZ10::unit(1) + ZZ10::unit(-1)) * ZZ10::from(2) - ZZ10::one()
    }
    pub fn zz20_sqrt5() -> ZZ20 {
        (ZZ20::unit(2) + ZZ20::unit(-2)) * ZZ20::from(2) - ZZ20::one()
    }

    pub fn zz24_sqrt6() -> ZZ24 {
        (ZZ24::unit(1) + ZZ24::unit(-1)) * ZZ24::from(2) - zz24_sqrt2()
    }
}

#[cfg(test)]
mod tests {
    use constants::zz_units_sum;

    use super::*;
    use crate::zzbase::{signum_sum_sqrt_expr_2, signum_sum_sqrt_expr_4};

    #[test]
    fn test_constants() {
        use super::constants::*;
        use std::f64::consts::SQRT_2;

        let sq2 = SQRT_2;
        let sq3 = 3.0_f64.sqrt();
        let sq_penta = ZZ10_Y.sqrt();
        let hsq_penta = 0.5 * ZZ10_Y.sqrt();
        let sq5 = 5.0_f64.sqrt();
        let sq6 = 6.0_f64.sqrt();

        assert_eq!(zz8_sqrt2().complex().re, sq2);
        assert_eq!(zz16_sqrt2().complex().re, sq2);
        assert_eq!(zz24_sqrt2().complex().re, sq2);
        assert_eq!(zz32_sqrt2().complex().re, sq2);

        assert_eq!(zz6_isqrt3().complex().im, sq3);
        assert_eq!(zz12_sqrt3().complex().re, sq3);
        assert_eq!(zz24_sqrt3().complex().re, sq3);

        assert_eq!(zz10_isqrt_penta().complex().im, sq_penta);
        assert_eq!(zz20_half_sqrt_penta().complex().re, hsq_penta);

        assert_eq!(zz10_sqrt5().complex().re, sq5);
        assert_eq!(zz20_sqrt5().complex().re, sq5);

        assert_eq!(zz24_sqrt6().complex().re, sq6);
    }

    #[test]
    fn test_sum_root_expr_sign_2() {
        assert_eq!(signum_sum_sqrt_expr_2(0, 2, 0, 3), 0);
        assert_eq!(signum_sum_sqrt_expr_2(1, 2, 0, 3), 1);
        assert_eq!(signum_sum_sqrt_expr_2(0, 2, -1, 3), -1);
        assert_eq!(signum_sum_sqrt_expr_2(2, 2, -1, 3), 1);
        assert_eq!(signum_sum_sqrt_expr_2(-5, 2, 4, 3), -1);
        assert_eq!(signum_sum_sqrt_expr_2(-5, 2, 5, 3), 1);
    }

    #[test]
    fn test_sum_root_expr_sign_4() {
        let sign_zz24 = |a, b, c, d| signum_sum_sqrt_expr_4(a, 1, b, 2, c, 3, d, 6);
        // trivial sanity-checks
        assert_eq!(sign_zz24(0, 0, 0, 0), 0);
        assert_eq!(sign_zz24(1, 1, 1, 1), 1);
        assert_eq!(sign_zz24(-1, -1, -1, -1), -1);
        assert_eq!(sign_zz24(1, 0, 0, 0), 1);
        assert_eq!(sign_zz24(0, -1, 0, 0), -1);
        assert_eq!(sign_zz24(0, 0, 1, 0), 1);
        assert_eq!(sign_zz24(0, 0, 0, -1), -1);
        // non-trivial tests
        assert_eq!(sign_zz24(5, 7, 11, -13), 1);
        assert_eq!(sign_zz24(5, 7, 11, -14), -1);
        assert_eq!(sign_zz24(17, -11, 9, -7), -1);
        assert_eq!(sign_zz24(18, -11, 9, -7), 1);
        assert_eq!(sign_zz24(18, -11, 8, -7), -1);
        assert_eq!(sign_zz24(18, -11, 8, -6), 1);

        // try with parameters where terms are all really close
        {
            let (a, b, c, d) = (130, 92, 75, 53);
            assert_eq!(sign_zz24(-a, -b, c, d), -1);
            assert_eq!(sign_zz24(-a, b, -c, d), 1);
            assert_eq!(sign_zz24(-a, b, c, -d), 1);
            assert_eq!(sign_zz24(a, -b, -c, d), -1);
            assert_eq!(sign_zz24(a, -b, c, -d), -1);
            assert_eq!(sign_zz24(a, b, -c, -d), 1);
        }
        {
            let (a, b, c, d) = (485, 343, 280, 198);
            assert_eq!(sign_zz24(-a, -b, c, d), -1);
            assert_eq!(sign_zz24(-a, b, -c, d), 1);
            assert_eq!(sign_zz24(-a, b, c, -d), 1);
            assert_eq!(sign_zz24(a, -b, -c, d), -1);
            assert_eq!(sign_zz24(a, -b, c, -d), -1);
            assert_eq!(sign_zz24(a, b, -c, -d), 1);
        }
    }

    macro_rules! zz_tests {
    ($($name:ident: $type:ty,)*) => {$(
mod $name {
    use super::*;

    type ZZi = $type;

    #[test]
    fn test_units_sum_is_complex() {
        let p: ZZi = zz_units_sum();
        assert!(p.is_complex());
    }

    #[test]
    fn test_basic() {
        // # of full turn steps is an even natural (so a half turn is possible)
        assert!(ZZi::turn() > 0);
        assert!(ZZi::turn() % 2 == 0);

        // check vector sizes
        let roots_num = ZZi::zz_params().sym_roots_num;
        assert_eq!(ZZi::zz_params().sym_roots_sqs.len(), roots_num);
        assert_eq!(ZZi::zz_params().ccw_unit_coeffs.len(), roots_num);

        // check zero-vector
        let z = ZZi::zero();
        let cs = z.zz_coeffs();
        assert_eq!(cs.len(), roots_num);
        for i in 0..roots_num {
            assert_eq!(cs[i], GInt::zero());
        }

        // check one-vector
        let o = ZZi::one();
        let cs = o.zz_coeffs();
        assert_eq!(cs.len(), roots_num);
        assert_eq!(cs[0], GInt::one());
        for i in 1..roots_num {
            assert_eq!(cs[i], GInt::zero());
        }
    }

    #[test]
    fn test_add_sub() {
        // test addition / subtraction
        assert_eq!(ZZi::zero() + ZZi::zero(), ZZi::zero());
        assert_eq!(ZZi::one() + ZZi::zero(), ZZi::one());
        assert_eq!(ZZi::zero() + ZZi::one(), ZZi::one());
        assert_eq!(-ZZi::one() + ZZi::one(), ZZi::zero());
        assert_eq!(ZZi::one() - ZZi::one(), ZZi::zero());
    }

    #[test]
    fn test_mul() {
        // test scalar multiplication
        assert_eq!(ZZi::zero().scale(2), ZZi::zero());
        assert_eq!(ZZi::ccw().scale(3), ZZi::ccw() + ZZi::ccw() + ZZi::ccw());
        assert_eq!(ZZi::one().scale(-1), -ZZi::one());
        assert_eq!(ZZi::one().scale(-42), ZZi::from(-42));

        // test multiplication
        assert_eq!(ZZi::zero() * ZZi::zero(), ZZi::zero());
        assert_eq!(ZZi::one() * ZZi::zero(), ZZi::zero());
        assert_eq!(ZZi::zero() * ZZi::one(), ZZi::zero());
        assert_eq!(ZZi::one() * ZZi::one(), ZZi::one());
        assert_eq!(-ZZi::one() * ZZi::one(), -ZZi::one());
        assert_eq!(ZZi::one() * (-ZZi::one()), -ZZi::one());
        assert_eq!((-ZZi::one()) * (-ZZi::one()), ZZi::one());
    }

    #[test]
    fn test_rotations() {
        // test ccw()
        assert_eq!(ZZi::ccw() * ZZi::ccw().conj(), ZZi::one());
        assert_eq!(-(-(ZZi::one()) * ZZi::ccw()), ZZi::ccw());

        // test going around the unit circle step by step
        let mut x = ZZi::one();
        for _ in 0..ZZi::turn() {
            x = x * ZZi::ccw();
        }
        assert_eq!(x, ZZi::one());

        // test unit()
        assert_eq!(ZZi::unit(0), ZZi::one());
        assert_eq!(ZZi::unit(-1), ZZi::unit(ZZi::turn() - 1));
        assert_eq!(ZZi::unit(1), ZZi::unit(ZZi::turn() + 1));
        assert_eq!(ZZi::unit(-ZZi::hturn()), ZZi::unit(ZZi::hturn()));
        assert_eq!(ZZi::unit(ZZi::hturn()), -ZZi::one());
        if ZZi::turn() % 4 == 0 {
            assert_eq!(ZZi::one_i().zz_coeffs()[0], GInt::from((0, 1)));
        }

        // test powi()
        assert_eq!(ZZi::ccw().powi(ZZi::hturn()), -ZZi::one());
        assert_eq!(ZZi::ccw().powi(ZZi::turn()), ZZi::one());
        assert_eq!(ZZi::ccw().powi(ZZi::hturn()).powi(2), ZZi::one());
    }

    #[test]
    #[should_panic]
    fn test_neg_powi() {
        ZZi::one().powi(-1);
    }

    #[test]
    fn test_scaling_fac() {
        // test scaling fac is correct by checking denom. of coeffs of all units
        // (that the denom. always can be expressed as multple of scaling factor)
        // and that the chosen constant factor is indeed minimal
        let sc_fac = ZZi::zz_params().scaling_fac;
        let mut max_fac: i64 = 0;
        for i in 0..ZZi::turn() {
            let x = ZZi::unit(i);
            println!("{x}");
            for c in x.coeffs {
                assert_eq!(sc_fac % c.real.denom(), 0);
                assert_eq!(sc_fac % c.imag.denom(), 0);
                max_fac = max_fac.max(*c.real.denom());
                max_fac = max_fac.max(*c.imag.denom());
            }
        }
        assert_eq!(sc_fac, max_fac);
    }

    #[test]
    fn test_xy() {
        if !ZZi::has_qturn() {
            return;
        }

        // test correctness of splitting and reconstruction
        for a in 0..ZZi::hturn() {
            let p = ZZi::unit(a);
            let (x, y) = p.xy();
            assert_eq!(p, x + y * ZZi::one_i());
        }
    }

    #[test]
    fn test_dot() {
        // get a non-trivial point
        let p: ZZi = zz_units_sum();

        // all rotations of the same point around origin
        // have the same squared distance, i.e. quadrance
        // and it is real-valued.
        let q = p.norm_sq();
        assert!(q.is_real());
        for i in 1..ZZi::turn() {
            let pi = p * ZZi::unit(i);
            let qi = pi.norm_sq();
            assert_eq!(qi, q);
        }
    }

    #[test]
    fn test_colinear() {
        if !ZZi::has_qturn() {
            return;
        }

        // Test point layout:
        // -------
        // E F
        //
        // A B C D
        // -------
        let a: ZZi = ZZi::zero();
        let b: ZZi = ZZi::one();
        let c: ZZi = ZZi::from(2);
        let d: ZZi = ZZi::from(3);
        let e: ZZi = ZZi::one_i();
        let f: ZZi = b + e;

        let l_ab = a.line_through(&b);
        let l_ac = a.line_through(&c);
        let l_af = a.line_through(&f);

        // colinear, overlap
        assert!(b.is_colinear(&l_ac));
        assert!(d.is_colinear(&l_ac));
        // colinear, no overlap
        assert!(c.is_colinear(&l_ab));
        assert!(d.is_colinear(&l_ab));
        // parallel (not colinear)
        assert!(!e.is_colinear(&l_ab));
        assert!(!f.is_colinear(&l_ab));
        // perpendicular (touches in one point, not in the other)
        assert!(!(a.is_colinear(&l_ab) && e.is_colinear(&l_ab)));
        // general case
        assert!(!b.is_colinear(&l_af));
        assert!(!d.is_colinear(&l_af));
    }

    #[test]
    fn test_is_real_imag_complex() {
        assert!(ZZi::zero().is_real());
        assert!(ZZi::zero().is_imag());

        assert!(ZZi::one().is_real());
        assert!((-ZZi::one()).is_real());
        assert!(!ZZi::one().is_imag());
        assert!(!ZZi::one().is_complex());

        if ZZi::has_qturn() && ZZi::qturn() > 1 {
            let p = ZZi::ccw();
            assert!(!p.is_real());
            assert!(!p.is_imag());
            assert!(p.is_complex());
        }

        if ZZi::has_qturn() {
            assert!(!ZZi::one_i().is_real());
            assert!(ZZi::one_i().is_imag());
            assert!((-ZZi::one_i()).is_imag());
            assert!(!ZZi::one_i().is_complex());
        }
    }

    #[test]
    fn test_complex() {
        use num_complex::Complex64;

        let x = ZZi::zero();
        assert_eq!(x.complex(), Complex64::zero());
        let x = ZZi::one();
        assert_eq!(x.complex(), Complex64::one());
        let x = -ZZi::one();
        assert_eq!(x.complex(), -Complex64::one());
        let x = ZZi::one() + ZZi::one();
        assert_eq!(x.complex(), Complex64::new(2.0, 0.0));

        let x = ZZi::ccw();
        let c = x.complex();
        println!("{c} = {x}");
    }

    #[test]
    fn test_hashable() {
        use std::collections::HashSet;

        let mut s: HashSet<ZZi> = HashSet::new();
        s.insert(ZZi::zero());
        s.insert(ZZi::one());
        s.insert(ZZi::ccw());
        assert!(s.contains(&ZZi::ccw()));
        assert!(s.contains(&(ZZi::ccw() + ZZi::zero())));
        assert!(!s.contains(&(ZZi::ccw() + ZZi::one())));
    }
}
    )*}
}
    zz_tests! {
        zz4: ZZ4,
        zz6: ZZ6,
        zz8: ZZ8,
        zz10: ZZ10,
        zz12: ZZ12,
        zz16: ZZ16,
        zz20: ZZ20,
        zz24: ZZ24,
        zz30: ZZ30,
        zz32: ZZ32,
        zz60: ZZ60,
    }

    #[test]
    fn test_re_signum() {
        use super::constants::*;

        let sq2 = zz24_sqrt2();
        let sq3 = zz24_sqrt3();
        let sq6 = zz24_sqrt6();

        let z = Frac::zero();
        let p = Frac::one();
        let m = -p;

        // use same test as above
        let sign_zz24 = |a, b, c, d| {
            (ZZ24::from(a) + ZZ24::from(b) * sq2 + ZZ24::from(c) * sq3 + ZZ24::from(d) * sq6)
                .re_signum()
        };

        let (a, b, c, d) = (485, 343, 280, 198);
        assert_eq!(sign_zz24(0, 0, 0, 0), z);
        assert_eq!(sign_zz24(-a, -b, c, d), m);
        assert_eq!(sign_zz24(-a, b, -c, d), p);
        assert_eq!(sign_zz24(-a, b, c, -d), p);
        assert_eq!(sign_zz24(a, -b, -c, d), m);
        assert_eq!(sign_zz24(a, -b, c, -d), m);
        assert_eq!(sign_zz24(a, b, -c, -d), p);
    }

    #[test]
    fn test_display() {
        let x = ZZ24::zero();
        assert_eq!(format!("{x}"), "0");

        let x = ZZ24::one();
        assert_eq!(format!("{x}"), "1");

        let x = ZZ24::one() + ZZ24::one();
        assert_eq!(format!("{x}"), "2");

        let x = -ZZ24::one();
        assert_eq!(format!("{x}"), "-1");

        let x = ZZ24::one() + (ZZ24::ccw()).powi(2);
        assert_eq!(format!("{x}"), "1+1/2i + (1/2)*sqrt(3)");

        let x: ZZ10 = zz_units_sum();
        assert_eq!(
            format!("{x}"),
            "-5 + (-15/4i)*sqrt(2(5-sqrt(5))) + (-5/4i)*sqrt(10(5-sqrt(5)))"
        );
    }

    #[test]
    fn test_dot_extra() {
        let p1 = ZZ12::one();
        let p2 = ZZ12::from(2);
        let p3 = ZZ12::from(3);
        let pi = ZZ12::unit(3); // i
        let p60 = ZZ12::unit(2);
        let pm60 = ZZ12::unit(-2);

        assert_eq!(p1.dot(&pi), ZZ12::zero());
        assert_eq!(ZZ12::zero().dot(&pi), ZZ12::zero());
        assert_eq!(p2.dot(&p3), ZZ12::from(6));

        // {0, 1} dot {1/2, sqrt(3)/2} = sqrt(3) / 2
        // => dot^2 = 3/4
        let d1 = p60.dot(&pi).powi(2).complex();
        assert_eq!(d1.re, 0.75);
        assert_eq!(d1.im, 0.0);

        // same but with negative sign (check indirectly)
        let d2 = pm60.dot(&pi).complex();
        assert!(d2.re < 0.0);
        assert_eq!(d2.im, 0.0);
        assert_eq!(pm60.dot(&pi).powi(2).complex().re, 0.75);
    }

    #[test]
    fn test_is_between() {
        type ZZi = ZZ24;

        let a: ZZi = ZZi::zero();
        let b: ZZi = ZZi::one();
        let c: ZZi = ZZi::from(2);
        let e: ZZi = ZZi::unit(ZZi::hturn() / 2);
        let f: ZZi = b + e;
        let g: ZZi = ZZi::unit(1) + ZZi::unit(-1) - ZZi::one();

        // is actually in between
        assert!(b.is_between(&a, &c));
        assert!(g.is_between(&a, &b));
        // is an endpoint
        assert!(!b.is_between(&a, &b));
        assert!(!a.is_between(&a, &b));
        // colinear, but not between
        assert!(!c.is_between(&a, &b));
        // not colinear
        assert!(!f.is_between(&a, &b));
    }
}
