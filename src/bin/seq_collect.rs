use std::collections::BTreeSet;
use std::sync::Arc;
use std::time::Instant;

use clap::{Parser, Subcommand};
use tilezz::cyclotomic::{ZZ10, ZZ12};
use tilezz::intgeom::matchtypes::MatchFinder;
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::seq_explorer::{replay_glue, Provenance, SeqExplorer};
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;

#[derive(Parser)]
#[command(name = "seq_collect", about = "SeqExplorer: collect + validate")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Collect {
        #[arg(long, value_enum, default_value = "hex")]
        tile: TileSetKind,
        #[arg(long)]
        output: String,
    },
    Validate {
        #[arg(long)]
        input: String,
    },
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum TileSetKind {
    Hex,
    Square,
    Mixed,
    Spectre,
    Penrose,
}

fn make_ts_12(kind: &TileSetKind) -> Arc<TileSet<ZZ12>> {
    match kind {
        TileSetKind::Hex => {
            let rat = Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![rat]))
        }
        TileSetKind::Square => {
            let rat = Rat::try_from(&tiles::square::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![rat]))
        }
        TileSetKind::Mixed => {
            let sq = Rat::try_from(&tiles::square::<ZZ12>()).unwrap();
            let hex = Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![sq, hex]))
        }
        TileSetKind::Spectre => {
            let rat = Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![rat]))
        }
        TileSetKind::Penrose => panic!("penrose requires ZZ10"),
    }
}

fn make_ts_10() -> Arc<TileSet<ZZ10>> {
    let n = Rat::try_from(&tiles::penrose_p3_narrow()).unwrap();
    let w = Rat::try_from(&tiles::penrose_p3_wide()).unwrap();
    Arc::new(TileSet::new(vec![n, w]))
}

fn write_collection<
    T: tilezz::cyclotomic::IsComplex + tilezz::cyclotomic::IsRingOrField + tilezz::cyclotomic::Units,
>(
    explorer: &SeqExplorer<T>,
    ring_name: &str,
    path: &str,
) {
    let t0 = Instant::now();
    let mut out = String::new();

    out.push_str(&format!("TILESET {}\n", ring_name));
    for i in 0..explorer.tileset().num_tiles() {
        let seq = explorer.tileset().rat(i).seq();
        let angles: Vec<String> = seq.iter().map(|a| a.to_string()).collect();
        out.push_str(&format!("TILE {} {}\n", i, angles.join(" ")));
    }

    for rat_id in 0..explorer.num_known_rats() {
        let prov = explorer.provenance(rat_id);
        match prov {
            Provenance::Seed { tile_idx } => {
                out.push_str(&format!("RAT {} seed {}\n", rat_id, tile_idx));
            }
            Provenance::Glue {
                source_rat_id,
                tile_idx,
                start_a,
                start_b,
                len,
            } => {
                out.push_str(&format!(
                    "RAT {} glue {} {} {} {} {}\n",
                    rat_id, source_rat_id, tile_idx, start_a, start_b, len
                ));
            }
        }
    }

    let collected = explorer.export_collected();
    let mut subseq_count = 0usize;
    for (rat_id, seqs) in &collected {
        for seq in seqs {
            let angles: Vec<String> = seq.iter().map(|a| a.to_string()).collect();
            out.push_str(&format!("SUBSEQ {} {}\n", rat_id, angles.join(" ")));
            subseq_count += 1;
        }
    }

    std::fs::write(path, &out).unwrap();
    eprintln!(
        "  Written {} bytes, {} subseqs in {:.2?}",
        out.len(),
        subseq_count,
        t0.elapsed()
    );
}

struct ParsedFile {
    ring_name: String,
    tile_seqs: Vec<Vec<i8>>,
    rat_entries: Vec<RatEntry>,
    subseqs: Vec<(usize, Vec<i8>)>,
}

enum RatEntry {
    Seed {
        tile_idx: usize,
    },
    Glue {
        source_rat_id: usize,
        tile_idx: usize,
        start_a: usize,
        start_b: usize,
        len: usize,
    },
}

fn parse_file(path: &str) -> Result<ParsedFile, String> {
    let content = std::fs::read_to_string(path).map_err(|e| format!("read error: {}", e))?;
    let lines: Vec<&str> = content.lines().collect();

    let mut ring_name = String::new();
    let mut tile_seqs: Vec<Vec<i8>> = Vec::new();
    let mut rat_entries: Vec<RatEntry> = Vec::new();
    let mut subseqs: Vec<(usize, Vec<i8>)> = Vec::new();

    for line in &lines {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }
        match parts[0] {
            "TILESET" => {
                ring_name = parts[1].to_string();
            }
            "TILE" => {
                let _id: usize = parts[1].parse().unwrap();
                let angles: Vec<i8> = parts[2..].iter().map(|s| s.parse().unwrap()).collect();
                tile_seqs.push(angles);
            }
            "RAT" => {
                let id: usize = parts[1].parse().unwrap();
                while rat_entries.len() <= id {
                    rat_entries.push(RatEntry::Seed { tile_idx: 0 });
                }
                match parts[2] {
                    "seed" => {
                        let tile_idx: usize = parts[3].parse().unwrap();
                        rat_entries[id] = RatEntry::Seed { tile_idx };
                    }
                    "glue" => {
                        let source_rat_id: usize = parts[3].parse().unwrap();
                        let tile_idx: usize = parts[4].parse().unwrap();
                        let start_a: usize = parts[5].parse().unwrap();
                        let start_b: usize = parts[6].parse().unwrap();
                        let len: usize = parts[7].parse().unwrap();
                        rat_entries[id] = RatEntry::Glue {
                            source_rat_id,
                            tile_idx,
                            start_a,
                            start_b,
                            len,
                        };
                    }
                    _ => return Err(format!("unknown RAT provenance: {}", parts[2])),
                }
            }
            "SUBSEQ" => {
                let rat_id: usize = parts[1].parse().unwrap();
                let angles: Vec<i8> = parts[2..].iter().map(|s| s.parse().unwrap()).collect();
                subseqs.push((rat_id, angles));
            }
            _ => {}
        }
    }

    Ok(ParsedFile {
        ring_name,
        tile_seqs,
        rat_entries,
        subseqs,
    })
}

fn is_cyclic_substring(haystack: &[i8], needle: &[i8]) -> bool {
    let n = haystack.len();
    let m = needle.len();
    if m == 0 || m > n {
        return false;
    }
    for start in 0..n {
        let mut ok = true;
        for l in 0..m {
            if haystack[(start + l) % n] != needle[l] {
                ok = false;
                break;
            }
        }
        if ok {
            return true;
        }
    }
    false
}

fn validate_common<
    T: tilezz::cyclotomic::IsComplex + tilezz::cyclotomic::IsRingOrField + tilezz::cyclotomic::Units,
>(
    pf: &ParsedFile,
    tile_ts: &Arc<TileSet<T>>,
) -> Result<(), String> {
    let t0 = Instant::now();
    let k = tile_ts.rats().iter().map(|r| r.len()).max().unwrap_or(0);

    eprintln!(
        "  Parsed: ring={}, tiles={}, rats={}, subseqs={}, k={}",
        pf.ring_name,
        pf.tile_seqs.len(),
        pf.rat_entries.len(),
        pf.subseqs.len(),
        k,
    );

    eprintln!("  Replaying glues...");
    let mut rats: Vec<Rat<T>> = Vec::with_capacity(pf.rat_entries.len());
    for (id, entry) in pf.rat_entries.iter().enumerate() {
        let rat = match entry {
            RatEntry::Seed { tile_idx } => tile_ts.rat(*tile_idx).clone(),
            RatEntry::Glue {
                source_rat_id,
                tile_idx,
                start_a,
                start_b,
                len,
            } => {
                let source = &rats[*source_rat_id];
                let tile = tile_ts.rat(*tile_idx);
                replay_glue(source, tile, *start_a, *start_b, *len)
                    .map_err(|e| format!("RAT {}: glue failed: {}", id, e))?
            }
        };
        rats.push(rat);
    }
    eprintln!("  Replayed {} rats in {:.2?}", rats.len(), t0.elapsed());

    eprintln!("  Checking subseq presence...");
    let t1 = Instant::now();
    let mut subseq_errors = 0usize;
    for (rat_id, seq) in &pf.subseqs {
        if !is_cyclic_substring(rats[*rat_id].seq(), seq) {
            subseq_errors += 1;
            if subseq_errors <= 5 {
                eprintln!(
                    "  ERROR: subseq {:?} not found in rat {} (len={})",
                    seq,
                    rat_id,
                    rats[*rat_id].len()
                );
            }
        }
    }
    if subseq_errors > 0 {
        return Err(format!("{} subseq presence errors", subseq_errors));
    }
    eprintln!(
        "  All {} subseqs verified present in {:.2?}",
        pf.subseqs.len(),
        t1.elapsed()
    );

    eprintln!("  Checking completeness (witnesses × tiles, batched)...");
    let t2 = Instant::now();

    let witness_ids: Vec<usize> = {
        let mut ids: BTreeSet<usize> = pf.subseqs.iter().map(|(id, _)| *id).collect();
        ids.into_iter().collect()
    };
    let witness_rats: Vec<Rat<T>> = witness_ids.iter().map(|&id| rats[id].clone()).collect();

    let mut known_subseqs: BTreeSet<Vec<i8>> = BTreeSet::new();
    for (_, seq) in &pf.subseqs {
        known_subseqs.insert(seq.clone());
    }

    let batch_size = 500;
    let mut new_found = 0usize;
    let mut total_matches = 0usize;
    let num_batches = witness_rats.len().div_ceil(batch_size);

    for (batch_num, chunk) in witness_rats.chunks(batch_size).enumerate() {
        let chunk_ts = Arc::new(TileSet::new(chunk.to_vec()));
        let mf = MatchFinder::crossing(Arc::clone(&chunk_ts), Arc::clone(tile_ts));

        for rat_idx in 0..chunk_ts.num_tiles() {
            for tile_idx in 0..tile_ts.num_tiles() {
                for mt in mf.valid_matches(rat_idx, tile_idx) {
                    total_matches += 1;
                    let glued = mt.apply(chunk_ts.rats(), tile_ts.rats());
                    let seq = glued.seq();
                    let n = seq.len();
                    let max_len = k.min(n);
                    for start in 0..n {
                        for l in 1..=max_len {
                            let subseq: Vec<i8> = (0..l).map(|j| seq[(start + j) % n]).collect();
                            if !known_subseqs.contains(&subseq) {
                                new_found += 1;
                                if new_found <= 10 {
                                    eprintln!(
                                        "  MISSING: subseq {:?} from glued rat (len={}) not in collection",
                                        subseq, n
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }

        if (batch_num + 1) % 10 == 0 || batch_num + 1 == num_batches {
            eprintln!(
                "    batch {}/{}: matches={} new={} elapsed={:.1}s",
                batch_num + 1,
                num_batches,
                total_matches,
                new_found,
                t2.elapsed().as_secs_f64()
            );
        }
    }

    eprintln!(
        "  Completeness: {} matches checked, {} new subseqs found in {:.2?}",
        total_matches,
        new_found,
        t2.elapsed()
    );

    if new_found > 0 {
        return Err(format!(
            "Completeness check FAILED: {} new subseqs not in collection",
            new_found
        ));
    }

    eprintln!("  Validation PASSED in {:.2?}", t0.elapsed());
    Ok(())
}

fn main() {
    #[cfg(feature = "pprof")]
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
        .unwrap();

    let cli = Cli::parse();

    match &cli.command {
        Commands::Collect { tile, output } => match tile {
            TileSetKind::Penrose => {
                let ts = make_ts_10();
                let t0 = Instant::now();
                let explorer = SeqExplorer::new(ts);
                let elapsed = t0.elapsed();
                eprintln!(
                    "[penrose] subseqs={} witnesses={} rats={} k={} time={:.2?}",
                    explorer.num_subseqs(),
                    explorer.num_witnesses(),
                    explorer.num_known_rats(),
                    explorer.max_subseq_len(),
                    elapsed,
                );
                write_collection(&explorer, "ZZ10", output);
            }
            _ => {
                let ts = make_ts_12(tile);
                let t0 = Instant::now();
                let explorer = SeqExplorer::new(ts);
                let elapsed = t0.elapsed();
                eprintln!(
                    "[{:?}] subseqs={} witnesses={} rats={} k={} time={:.2?}",
                    tile,
                    explorer.num_subseqs(),
                    explorer.num_witnesses(),
                    explorer.num_known_rats(),
                    explorer.max_subseq_len(),
                    elapsed,
                );
                write_collection(&explorer, "ZZ12", output);
            }
        },
        Commands::Validate { input } => {
            eprintln!("=== Validating: {} ===", input);
            let pf = parse_file(input).unwrap_or_else(|e| {
                eprintln!("Parse error: {}", e);
                std::process::exit(1);
            });

            match pf.ring_name.as_str() {
                "ZZ12" => {
                    let tile_ts = Arc::new(TileSet::new(
                        pf.tile_seqs
                            .iter()
                            .map(|s| Rat::<ZZ12>::from_slice_unchecked(s))
                            .collect(),
                    ));
                    if let Err(e) = validate_common(&pf, &tile_ts) {
                        eprintln!("FAIL: {}", e);
                        std::process::exit(1);
                    }
                }
                "ZZ10" => {
                    let tile_ts = Arc::new(TileSet::new(
                        pf.tile_seqs
                            .iter()
                            .map(|s| Rat::<ZZ10>::from_slice_unchecked(s))
                            .collect(),
                    ));
                    if let Err(e) = validate_common(&pf, &tile_ts) {
                        eprintln!("FAIL: {}", e);
                        std::process::exit(1);
                    }
                }
                _ => {
                    eprintln!("Unsupported ring: {}", pf.ring_name);
                    std::process::exit(1);
                }
            }
            println!("OK");
        }
    }

    #[cfg(feature = "pprof")]
    {
        if let Ok(report) = guard.report().build() {
            let path = "flamegraph_seq.svg";
            let mut file = std::fs::File::create(path).unwrap();
            report.flamegraph(&mut file).unwrap();
            eprintln!("Flame graph written to {}", path);

            let mut text_report = Vec::new();
            use std::io::Write;
            report.text_detailed(&mut text_report).unwrap();
            let text_path = "pprof_seq.txt";
            std::fs::write(text_path, &text_report).unwrap();
            eprintln!("Text report written to {}", text_path);
        }
    }
}
