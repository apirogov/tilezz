use rational_tiles::zz::ZZ12;
use rational_tiles::zzbase::ZZBase;
use std::collections::HashSet;

fn main() {
    let mut vs: HashSet<ZZ12> = HashSet::new();
    for i in 0..ZZ12::turn() {
        vs.insert(ZZ12::unit(i));
    }
    for _ in 0..3 {
        let mut vs_new: HashSet<ZZ12> = HashSet::new();
        for v1 in vs.iter() {
            for v2 in vs.iter() {
                vs_new.insert(*v1 + *v2);
            }
        }
        vs = vs_new;
    }
}
