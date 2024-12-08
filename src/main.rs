use num_traits::{One, Zero};
use rational_tiles::gaussint::GaussInt32;
use rational_tiles::traits::Ccw;
use rational_tiles::zz::{ZZ12, ZZ4};

fn main() {
    let v0 = GaussInt32::zero();
    let v1 = GaussInt32::one();
    let v2 = GaussInt32::ccw();
    let v3 = v2 * v2;
    let v4 = v3 * v2;
    let v5 = v4 * v2;
    let v6 = v1 + v2;
    let v7 = v1 - v2;
    println!("{v0}\n{v1}\n{v2}\n{v3}\n{v4}\n{v5}\n{v6}\n{v7}\n");

    let z0 = ZZ4::zero();
    println!("bla");
    let z1 = ZZ4::one();
    let z2 = ZZ4::ccw();
    let z3 = z2 * z2;
    let z4 = z3 * z2;
    let z5 = z4 * z2;
    let z6 = z1 * z2;
    let z7 = z1 - z2;
    println!("{z0:?}\n{z1:?}\n{z2:?}\n{z3:?}\n{z4:?}\n{z5:?}\n{z6:?}\n{z7:?}\n");

    let tmp = ZZ12::ccw() + ZZ12::one();
    println!("{tmp}");
}
