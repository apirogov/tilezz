//! In this module the different ring types are generated using macros.

use std::f64::consts::SQRT_2;

use num_traits::FromPrimitive;

use crate::traits::IntRing;
use crate::zzbase::ZZParams;

// ----------------

// definitions needed to derive different ZZn types
//
// numeric constants for the chosen (in general non-unique) linearly independent algebraic bases
// --------
// NOTE: sqrt(2 + sqrt(2)) = (sqrt(2)+1)sqrt(2 - sqrt(2))
// and   sqrt(2 - sqrt(2)) = (sqrt(2)-1)sqrt(2 + sqrt(2))
pub const ZZ16_Y: f64 = 2.0 + SQRT_2;

// NOTE: sqrt(2+sqrt(2+sqrt(2))) = ( 1 + sqrt(2) + sqrt(2(2+sqrt(2)))) sqrt(2-sqrt(2+sqrt(2)))
// and   sqrt(2-sqrt(2+sqrt(2))) = (-1 - sqrt(2) + sqrt(2(2+sqrt(2)))) sqrt(2+sqrt(2+sqrt(2)))
pub const ZZ32_Z: f64 = 2.0 + 1.84775906502257351225; // 2+sqrt(2+sqrt(2))

// NOTE: Let x = sqrt(5), y = sqrt(2*(5-sqrt(5))), y' = sqrt(2*(5+sqrt(5)))
// We have: y = 1/2(xy' - y') ^ y' = 1/2(xy + y)
pub const SQRT_5: f64 = 2.23606797749978969;
pub const ZZ10_Y: f64 = 2.0 * (5.0 - SQRT_5);

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
/// simplified for each ring separately, which is more efficient.

// ----------------
/// Gauss integers
pub const ZZ4_PARAMS: ZZParams = ZZParams {
    full_turn_steps: 4,
    sym_roots_num: 1,
    sym_roots_sqs: &[1.],
    sym_roots_lbls: &["1"],
    scaling_fac: 1,
    ccw_unit_coeffs: &[[0, 1]],
};
pub fn zz4_mul<T: IntRing>(x: &[T], y: &[T]) -> Vec<T> {
    match [*array_ref!(x, 0, 1), *array_ref!(y, 0, 1)] {
        [[a], [b]] => vec![a * b],
    }
}

// ----------------
/// Eisenstein integers
pub const ZZ6_PARAMS: ZZParams = ZZParams {
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
pub fn zz6_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    let int3 = T::from_i64(3).unwrap();
    match [*array_ref!(x, 0, 2), *array_ref!(y, 0, 2)] {
        [[a, b], [c, d]] => vec![a * c + (b * d * int3), a * d + b * c],
    }
}

// ----------------
/// Compass integers
pub const ZZ8_PARAMS: ZZParams = ZZParams {
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
pub fn zz8_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    let int2 = T::from_i64(2).unwrap();
    match [*array_ref!(x, 0, 2), *array_ref!(y, 0, 2)] {
        [[a, b], [c, d]] => vec![a * c + (b * d * int2), a * d + b * c],
    }
}

// ----------------
/// Halfrose integers (impractical, has no quarter turn)
pub const ZZ10_PARAMS: ZZParams = ZZParams {
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
pub fn zz10_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    let i = |n: i8| T::from_i8(n).unwrap();
    let (int2, int5, int10) = (i(2), i(5), i(10));
    match [*array_ref!(x, 0, 4), *array_ref!(y, 0, 4)] {
        [[a, b, c, d], [e, f, g, h]] => {
            let c1 = a * e + b * f * int5 + (c * g + d * h * int5 - c * h - d * g) * int10;
            let c2 = a * f + b * e - c * g * int2 + (c * h + d * g - d * h) * int10;
            let c3 = a * g + (b * h + d * f) * int5 + c * e;
            let c4 = a * h + b * g + c * f + d * e;
            vec![c1, c2, c3, c4]
        }
    }
}

// ----------------
/// Clock integers
pub const ZZ12_PARAMS: ZZParams = ZZParams {
    full_turn_steps: 12,
    sym_roots_num: ZZ6_PARAMS.sym_roots_num,
    sym_roots_sqs: ZZ6_PARAMS.sym_roots_sqs,
    sym_roots_lbls: ZZ6_PARAMS.sym_roots_lbls,
    scaling_fac: ZZ6_PARAMS.scaling_fac,
    ccw_unit_coeffs: &[[0, 1], [1, 0]],
};
pub fn zz12_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    return zz6_mul(x, y);
}

// ----------------
/// Hex integers
pub const ZZ16_PARAMS: ZZParams = ZZParams {
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
pub fn zz16_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    let int2 = T::from_i64(2).unwrap();
    match [*array_ref!(x, 0, 4), *array_ref!(y, 0, 4)] {
        [[a, b, c, d], [e, f, g, h]] => {
            let c1 = a * e + (b * f + c * g + c * h + d * g + d * h * int2) * int2;
            let c2 = a * f + b * e + c * g + (c * h + d * g + d * h) * int2;
            let c3 = a * g + c * e + (b * h + d * f) * int2;
            let c4 = a * h + b * g + c * f + d * e;
            vec![c1, c2, c3, c4]
        }
    }
}

// ----------------
/// Penrose integers
pub const ZZ20_PARAMS: ZZParams = ZZParams {
    full_turn_steps: 20,
    sym_roots_num: ZZ10_PARAMS.sym_roots_num,
    sym_roots_sqs: ZZ10_PARAMS.sym_roots_sqs,
    sym_roots_lbls: ZZ10_PARAMS.sym_roots_lbls,
    scaling_fac: ZZ10_PARAMS.scaling_fac,
    ccw_unit_coeffs: &[[0, -2], [0, 2], [1, 0], [1, 0]],
};
pub fn zz20_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    zz10_mul(x, y)
}

// ----------------
/// Digiclock integers
pub const ZZ24_PARAMS: ZZParams = ZZParams {
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
pub fn zz24_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    let i = |n: i8| T::from_i8(n).unwrap();
    let (int2, int3, int6) = (i(2), i(3), i(6));
    match [*array_ref!(x, 0, 4), *array_ref!(y, 0, 4)] {
        [[a, b, c, d], [e, f, g, h]] => {
            let c1 = a * e + b * f * int2 + c * g * int3 + d * h * int6;
            let c2 = a * f + b * e + (c * h + d * g) * int3;
            let c3 = a * g + c * e + (b * h + d * f) * int2;
            let c4 = a * h + b * g + c * f + d * e;
            vec![c1, c2, c3, c4]
        }
    }
}

// ----------------
/// Month integers
pub const ZZ30_PARAMS: ZZParams = ZZParams {
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
pub fn zz30_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    let i = |n: i64| T::from_i64(n).unwrap();
    let int2 = i(2);
    let int3 = i(3);
    let int5 = i(5);
    let int6 = i(6);
    let int10 = i(10);
    let int15 = i(15);
    let int30 = i(30);
    let int50 = i(50);
    let int150 = i(150);
    match [*array_ref!(x, 0, 8), *array_ref!(y, 0, 8)] {
        [[l1, l2, l3, l4, l5, l6, l7, l8], [r1, r2, r3, r4, r5, r6, r7, r8]] => {
            // 1
            let c1_d = ((l1 * r1) + (l2 * r2) * int3 + (l3 * r3) * int5 + (l4 * r4) * int10)
                + ((l5 * r5) * int15 + (l6 * r6) * int30 + (l7 * r7) * int50 + (l8 * r8) * int150);
            let c1_a = l4 * r7 * (-int10) + l6 * r8 * (-int30);
            let c1_b = r4 * l7 * (-int10) + r6 * l8 * (-int30);
            let c1 = c1_d + c1_a + c1_b;

            // x
            let c2_a =
                l1 * r2 + l3 * r5 * int5 + (l4 * r6 - l4 * r8 - l6 * r7) * int10 + l7 * r8 * int50;
            let c2_b =
                r1 * l2 + r3 * l5 * int5 + (r4 * l6 - r4 * l8 - r6 * l7) * int10 + r7 * l8 * int50;
            let c2 = c2_a + c2_b;

            // y
            let c3_d =
                r4 * l4 * (-int2) + r6 * l6 * (-int6) + r7 * l7 * (-int10) + r8 * l8 * (-int30);
            let c3_a = l1 * r3 + l2 * r5 * int3 + l4 * r7 * int10 + l6 * r8 * int30;
            let c3_b = r1 * l3 + r2 * l5 * int3 + r4 * l7 * int10 + r6 * l8 * int30;
            let c3 = c3_d + c3_a + c3_b;

            // z
            let c4_a = l1 * r4 + l2 * r6 * int3 + l3 * r7 * int5 + l5 * r8 * int15;
            let c4_b = r1 * l4 + r2 * l6 * int3 + r3 * l7 * int5 + r5 * l8 * int15;
            let c4 = c4_a + c4_b;

            // xy
            let c5_a = l1 * r5 + l2 * r3 - l4 * r6 * int2 + (l4 * r8 + l6 * r7 - l7 * r8) * int10;
            let c5_b = r1 * l5 + r2 * l3 - r4 * l6 * int2 + (r4 * l8 + r6 * l7 - r7 * l8) * int10;
            let c5 = c5_a + c5_b;

            // xz
            let c6_a = l1 * r6 + l2 * r4 + (l3 * r8 + l5 * r7) * int5;
            let c6_b = r1 * l6 + r2 * l4 + (r3 * l8 + r5 * l7) * int5;
            let c6 = c6_a + c6_b;

            // yz
            let c7_a = l1 * r7 + l2 * r8 * int3 + l3 * r4 + l5 * r6 * int3;
            let c7_b = r1 * l7 + r2 * l8 * int3 + r3 * l4 + r5 * l6 * int3;
            let c7 = c7_a + c7_b;

            // xyz
            let c8 = r1 * l8 + r2 * l7 + r3 * l6 + r4 * l5 + r5 * l4 + r6 * l3 + r7 * l2 + r8 * l1;

            vec![c1, c2, c3, c4, c5, c6, c7, c8]
        }
    }
}

// ----------------
/// ??? integers
pub const ZZ32_PARAMS: ZZParams = ZZParams {
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
pub fn zz32_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    let i = |n: i8| T::from_i8(n).unwrap();
    let int2 = i(2);
    let int4 = i(4);
    match [*array_ref!(x, 0, 8), *array_ref!(y, 0, 8)] {
        [[l1, l2, l3, l4, l5, l6, l7, l8], [r1, r2, r3, r4, r5, r6, r7, r8]] => {
            // 1
            let c1_d = (l1 * r1)
                + (l2 * r2 + l3 * r3 + l4 * r4) * int2
                + (l5 * r5 + l6 * r6 + l7 * r7) * int4
                + l8 * r8 * i(8);
            let c1_a = (l3 * r5 + l4 * r7 + l4 * r8 + l6 * r7) * int2 + (l6 * r8 + l7 * r8) * int4;
            let c1_b = (r3 * l5 + r4 * l7 + r4 * l8 + r6 * l7) * int2 + (r6 * l8 + r7 * l8) * int4;
            let c1 = c1_d + c1_a + c1_b;

            // x
            let c2_d = r3 * l3 + (r5 * l5 + r7 * l7) * int2 + r8 * l8 * int4;
            let c2_a = (l1 * r2 + l4 * r7)
                + (l3 * r5 + l4 * r6 + l4 * r8 + l6 * r7 + l6 * r8) * int2
                + l7 * r8 * int4;
            let c2_b = (r1 * l2 + r4 * l7)
                + (r3 * l5 + r4 * l6 + r4 * l8 + r6 * l7 + r6 * l8) * int2
                + r7 * l8 * int4;
            let c2 = c2_d + c2_a + c2_b;

            // y
            let c3_d = l4 * r4 + (l6 * r6 + l7 * r7) * int2 + (l8 * r8) * int4;
            let c3_a = l1 * r3 + (l2 * r5 + l4 * r7 + l7 * r8) * int2 + (l6 * r8) * int4;
            let c3_b = r1 * l3 + (r2 * l5 + r4 * l7 + r7 * l8) * int2 + (r6 * l8) * int4;
            let c3 = c3_d + c3_a + c3_b;

            // z
            let c4_a = l1 * r4 + (l2 * r6 + l3 * r7 + l3 * r8 + l5 * r7) * int2 + l5 * r8 * int4;
            let c4_b = r1 * l4 + (r2 * l6 + r3 * l7 + r3 * l8 + r5 * l7) * int2 + r5 * l8 * int4;
            let c4 = c4_a + c4_b;

            // xy
            let c5_d = l7 * r7 + l8 * r8 * int2;
            let c5_a = l1 * r5 + l2 * r3 + l4 * r6 + (l4 * r8 + l6 * r7 + l7 * r8) * int2;
            let c5_b = r1 * l5 + r2 * l3 + r4 * l6 + (r4 * l8 + r6 * l7 + r7 * l8) * int2;
            let c5 = c5_d + c5_a + c5_b;

            // xz
            let c6_a = l1 * r6 + l2 * r4 + l3 * r7 + (l3 * r8 + l5 * r7 + l5 * r8) * int2;
            let c6_b = r1 * l6 + r2 * l4 + r3 * l7 + (r3 * l8 + r5 * l7 + r5 * l8) * int2;
            let c6 = c6_a + c6_b;

            // yz
            let c7_a = l1 * r7 + l3 * r4 + (l2 * r8 + l5 * r6) * int2;
            let c7_b = r1 * l7 + r3 * l4 + (r2 * l8 + r5 * l6) * int2;
            let c7 = c7_a + c7_b;

            // xyz
            let c8 = l1 * r8 + l2 * r7 + l3 * r6 + l4 * r5 + l5 * r4 + l6 * r3 + l7 * r2 + l8 * r1;

            vec![c1, c2, c3, c4, c5, c6, c7, c8]
        }
    }
}

// ----------------
/// Babylonian integers
pub const ZZ60_PARAMS: ZZParams = ZZParams {
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
pub fn zz60_mul<T: IntRing + FromPrimitive>(x: &[T], y: &[T]) -> Vec<T> {
    zz30_mul(x, y)
}
