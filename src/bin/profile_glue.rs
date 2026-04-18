use std::collections::{BTreeSet, HashSet};
use std::time::Instant;

use tilezz::cyclotomic::ZZ12;
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;

fn main() {
    let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
    let ts1 = TileSet::new(vec![hex.clone()]);
    let size2: Vec<Rat<ZZ12>> = ts1
        .all_valid_glues()
        .iter()
        .map(|g| g.result.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    let ts2 = TileSet::new(size2.clone());
    let size4: Vec<Rat<ZZ12>> = ts2
        .all_valid_glues()
        .iter()
        .map(|g| g.result.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    eprintln!("size4 patches: {}", size4.len());
    profile_self_match(&size4);
}

fn is_single_edge_candidate(a: &[i8], ia: usize, b: &[i8], ib: usize) -> bool {
    let na = a.len();
    let nb = b.len();
    let left = a[ia] as i32 + b[ib] as i32;
    let right = a[(ia + 1) % na] as i32 + b[(ib + nb - 1) % nb] as i32;
    left > 0 && right > 0
}

fn junction_gap_positive(a: &[i8], ns: usize, mlen: usize, b: &[i8], ne: usize) -> bool {
    let na = a.len();
    let nb = b.len();
    let left = a[(ns + mlen) % na] as i32 + b[(ne + nb - mlen) % nb] as i32;
    let right = a[ns] as i32 + b[ne] as i32;
    left > 0 && right > 0
}

fn profile_self_match(patches: &[Rat<ZZ12>]) {
    let n = patches.len();
    let ts = TileSet::new(patches.to_vec());

    let mut total_cmi_query_ns: u64 = 0;
    let mut total_get_match_ns: u64 = 0;
    let mut total_snake_ns: u64 = 0;
    let mut total_cmi_seeds = 0usize;
    let mut total_single_seeds = 0usize;
    let mut total_dedup_saved = 0usize;
    let mut total_unique_tried = 0usize;
    let mut total_glue_ok = 0usize;
    let mut total_glue_fail = 0usize;
    let mut total_try_glue_ns: u64 = 0;
    let mut total_junction_reject = 0usize;

    let t_all = Instant::now();

    for i in 0..n {
        for j in 0..n {
            let a = ts.rat(i);
            let b = ts.rat(j);
            let n_a = a.len();
            let n_b = b.len();
            if n_a == 0 || n_b == 0 {
                continue;
            }

            let t = Instant::now();
            let cmi_matches = ts.shared_boundaries(i, j);
            total_cmi_query_ns += t.elapsed().as_nanos() as u64;

            let mut seen: HashSet<(i64, usize, i64)> = HashSet::new();

            let pair_cmi_seeds: Vec<(i64, i64)> = {
                let mut seeds = Vec::new();
                for m in &cmi_matches {
                    let pa = m.pos_a as i64;
                    let pb = m.pos_b as i64;
                    for k in 0..=m.len {
                        let ka = (pa + k as i64).rem_euclid(n_a as i64);
                        let kb = (pb - k as i64).rem_euclid(n_b as i64);
                        seeds.push((ka, kb));
                    }
                    seeds.push((
                        (pa - 1).rem_euclid(n_a as i64),
                        (pb + 1).rem_euclid(n_b as i64),
                    ));
                }
                seeds
            };

            let seq_a = a.seq();
            let seq_b = b.seq();
            let pair_single_seeds: Vec<(i64, i64)> = (0..n_a)
                .flat_map(|ia| {
                    (0..n_b).filter_map(move |ib| {
                        is_single_edge_candidate(seq_a, ia, seq_b, ib)
                            .then_some((ia as i64, ib as i64))
                    })
                })
                .collect();

            total_cmi_seeds += pair_cmi_seeds.len();
            total_single_seeds += pair_single_seeds.len();

            for &(ia, ib) in pair_cmi_seeds.iter().chain(pair_single_seeds.iter()) {
                let t1 = Instant::now();
                let (ns, len, ne) = a.get_match((ia, ib), b);
                total_get_match_ns += t1.elapsed().as_nanos() as u64;

                if len == 0 || !seen.insert((ns, len, ne)) {
                    total_dedup_saved += 1;
                    continue;
                }
                total_unique_tried += 1;

                if !junction_gap_positive(a.seq(), ns as usize, len, b.seq(), ne as usize) {
                    total_junction_reject += 1;
                    continue;
                }

                let t2 = Instant::now();
                let glue_ok = a.try_glue((ia, ib), b).is_ok();
                let glue_elapsed = t2.elapsed().as_nanos() as u64;
                total_try_glue_ns += glue_elapsed;

                if glue_ok {
                    total_glue_ok += 1;
                } else {
                    total_glue_fail += 1;
                }
            }
        }
    }

    let dt_all = t_all.elapsed();
    let total_seeds = total_cmi_seeds + total_single_seeds;

    eprintln!("\n=== Profile: {} x {} tiles ===", n, n);
    eprintln!("Total wall time:       {:.2?}", dt_all);
    eprintln!("CMI seeds:             {}", total_cmi_seeds);
    eprintln!("Single-edge seeds:     {}", total_single_seeds);
    eprintln!("Total seeds:           {}", total_seeds);
    eprintln!(
        "Dedup skipped:         {} ({:.1}% of seeds)",
        total_dedup_saved,
        100.0 * total_dedup_saved as f64 / total_seeds as f64
    );
    eprintln!("Unique matches tried:  {}", total_unique_tried);
    eprintln!(
        "Junction gap reject:   {} ({:.1}% of unique)",
        total_junction_reject,
        100.0 * total_junction_reject as f64 / total_unique_tried as f64
    );
    eprintln!("Glue OK:               {}", total_glue_ok);
    eprintln!("Glue FAIL (snake):     {}", total_glue_fail);
    eprintln!();
    eprintln!(
        "CMI query:             {:.2?} ({:.1}%)",
        std::time::Duration::from_nanos(total_cmi_query_ns),
        100.0 * total_cmi_query_ns as f64 / dt_all.as_nanos() as f64
    );
    eprintln!(
        "get_match (explicit):  {:.2?} ({:.1}%)",
        std::time::Duration::from_nanos(total_get_match_ns),
        100.0 * total_get_match_ns as f64 / dt_all.as_nanos() as f64
    );
    eprintln!(
        "try_glue (full):       {:.2?} ({:.1}%)",
        std::time::Duration::from_nanos(total_try_glue_ns),
        100.0 * total_try_glue_ns as f64 / dt_all.as_nanos() as f64
    );
    let accounted = total_cmi_query_ns + total_get_match_ns + total_try_glue_ns;
    eprintln!(
        "Other (seeds/iter):    {:.2?} ({:.1}%)",
        dt_all - std::time::Duration::from_nanos(accounted),
        100.0 * (dt_all.as_nanos() as u64 - accounted) as f64 / dt_all.as_nanos() as f64
    );
    eprintln!();
    eprintln!("NOTE: try_glue includes a redundant get_match (pub API).",);
    eprintln!(
        "  Pure Snake cost est: {:.2?}",
        std::time::Duration::from_nanos(total_try_glue_ns.saturating_sub(total_get_match_ns)),
    );
}
