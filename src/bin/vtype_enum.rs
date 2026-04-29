use std::collections::{BTreeSet, HashMap};
use std::sync::Arc;
use std::time::Instant;

use clap::{Parser, Subcommand};
use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ10, ZZ12};
use tilezz::intgeom::matchtypes::MatchTypeIndex;
use tilezz::intgeom::patch::{
    candidates_from_flat, EdgeInfo, GrowingPatch, PatchMatch, VertexType,
};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;
use tilezz::intgeom::vertextypes::VertexTypeIndex;

#[cfg(feature = "pprof")]
use pprof::ProfilerGuardBuilder;

#[derive(Parser)]
#[command(
    name = "vtype_enum",
    about = "Vertex type enumeration: collect + validate"
)]
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
        TileSetKind::Spectre => {
            let rat = Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![rat]))
        }
        TileSetKind::Mixed => {
            let sq = Rat::try_from(&tiles::square::<ZZ12>()).unwrap();
            let hex = Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![sq, hex]))
        }
        TileSetKind::Penrose => panic!("penrose requires ZZ10"),
    }
}

fn make_ts_10() -> Arc<TileSet<ZZ10>> {
    let n = Rat::try_from(&tiles::penrose_p3_narrow()).unwrap();
    let w = Rat::try_from(&tiles::penrose_p3_wide()).unwrap();
    Arc::new(TileSet::new(vec![n, w]))
}

struct ParsedVtype {
    id: usize,
    angle: i8,
    cw_tile_id: usize,
    cw_offset: usize,
    inner: Vec<EdgeInfo>,
    ccw_tile_id: usize,
    ccw_offset: usize,
    kind: String,
    cursed: bool,
}

struct ParsedWitness {
    vtype_id: usize,
    pos: usize,
    angles: Vec<i8>,
    edges: Vec<EdgeInfo>,
    candidates: Vec<PatchMatch>,
    inner_chains: Vec<Vec<EdgeInfo>>,
}

struct ParsedTransition {
    src_id: usize,
    dst_id: usize,
    start_a: usize,
    len: usize,
    start_b: usize,
    tile_id: usize,
}

struct ParsedFile {
    ring_name: String,
    tile_seqs: Vec<Vec<i8>>,
    vtypes: Vec<ParsedVtype>,
    witnesses: Vec<ParsedWitness>,
    transitions: Vec<ParsedTransition>,
}

fn write_collection<T: IsComplex + IsRingOrField + Units>(
    idx: &VertexTypeIndex<T>,
    ring_name: &str,
    path: &str,
) {
    let t0 = Instant::now();
    let mut out = String::new();
    let mut witness_keys: BTreeSet<String> = BTreeSet::new();

    out.push_str(&format!("TILESET {}\n", ring_name));
    for i in 0..idx.tileset().num_tiles() {
        let seq = idx.tileset().rat(i).seq();
        let angles: Vec<String> = seq.iter().map(|a| a.to_string()).collect();
        out.push_str(&format!("TILE {} {}\n", i, angles.join(" ")));
    }

    for id in 1..=idx.num_types() {
        let info = idx.get_info(id);
        let vt = info.vtype();
        let kind_str = if info.is_dead() { "dead" } else { "open" };
        out.push_str(&format!(
            "VTYPE {} {} {} {} {} {} {} {} {} {}\n",
            id,
            info.gap_angle(),
            vt.cw.tile_id,
            vt.cw.tile_offset,
            vt.inner.len(),
            vt.inner
                .iter()
                .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
                .collect::<Vec<_>>()
                .join(" "),
            vt.ccw.tile_id,
            vt.ccw.tile_offset,
            kind_str,
            if info.is_cursed() { 1 } else { 0 },
        ));

        let w = info.witness();
        let angles = w.angles();
        let edges = w.edges();
        let n = angles.len();
        let angle_strs: Vec<String> = angles.iter().map(|a| a.to_string()).collect();
        let edge_strs: Vec<String> = edges
            .iter()
            .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
            .collect();
        let all_matches = w.get_all_matches();
        let cand_strs: Vec<String> = all_matches
            .iter()
            .map(|m| format!("{}.{}.{}.{}", m.start_a, m.len, m.start_b, m.tile_id))
            .collect();
        let inner_strs: Vec<String> = w
            .inner_chains()
            .iter()
            .map(|chain| {
                if chain.is_empty() {
                    "-".to_string()
                } else {
                    chain
                        .iter()
                        .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
                        .collect::<Vec<_>>()
                        .join(",")
                }
            })
            .collect();

        let witness_key = format!(
            "{}|{}|{}|{}",
            angle_strs.join(" "),
            edge_strs.join(" "),
            cand_strs.join(" "),
            inner_strs.join(" "),
        );
        witness_keys.insert(witness_key);

        out.push_str(&format!(
            "WITNESS {} {} {} {} {} {} {} {}\n",
            id,
            info.witness_pos(),
            n,
            angle_strs.join(" "),
            edge_strs.join(" "),
            all_matches.len(),
            cand_strs.join(" "),
            inner_strs.join(" "),
        ));
    }

    for tr in idx.transitions() {
        let side_str = match tr.side {
            tilezz::intgeom::patch::TransitionSide::Cw => "cw",
            tilezz::intgeom::patch::TransitionSide::Ccw => "ccw",
        };
        out.push_str(&format!(
            "TRANS {} {} {} {} {} {} {}\n",
            tr.src_id,
            tr.dst_id,
            side_str,
            tr.patch_match.start_a,
            tr.patch_match.len,
            tr.patch_match.start_b,
            tr.patch_match.tile_id,
        ));
    }

    std::fs::write(path, &out).unwrap();
    let n_alive = idx.entries().iter().filter(|e| e.is_alive()).count();
    let n_dead = idx.entries().iter().filter(|e| e.is_dead()).count();
    let n_cursed = idx.entries().iter().filter(|e| e.is_cursed()).count();
    eprintln!(
        "  Written {} bytes, {} types ({} alive, {} dead, {} cursed), {} unique witnesses, {} transitions in {:.2?}",
        out.len(),
        idx.num_types(),
        n_alive,
        n_dead,
        n_cursed,
        witness_keys.len(),
        idx.transitions().len(),
        t0.elapsed(),
    );
}

fn parse_file(path: &str) -> Result<ParsedFile, String> {
    let content = std::fs::read_to_string(path).map_err(|e| format!("read error: {}", e))?;
    let mut ring_name = String::new();
    let mut tile_seqs: Vec<Vec<i8>> = Vec::new();
    let mut vtypes: Vec<ParsedVtype> = Vec::new();
    let mut witnesses: Vec<ParsedWitness> = Vec::new();
    let mut transitions: Vec<ParsedTransition> = Vec::new();

    for line in content.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }
        match parts[0] {
            "TILESET" => {
                ring_name = parts[1].to_string();
            }
            "TILE" => {
                let angles: Vec<i8> = parts[2..].iter().map(|s| s.parse().unwrap()).collect();
                tile_seqs.push(angles);
            }
            "VTYPE" => {
                let id: usize = parts[1].parse().map_err(|e| format!("VTYPE id: {}", e))?;
                let angle: i8 = parts[2]
                    .parse()
                    .map_err(|e| format!("VTYPE angle: {}", e))?;
                let cw_tile_id: usize = parts[3]
                    .parse()
                    .map_err(|e| format!("VTYPE cw_tid: {}", e))?;
                let cw_offset: usize = parts[4]
                    .parse()
                    .map_err(|e| format!("VTYPE cw_off: {}", e))?;
                let inner_len: usize = parts[5]
                    .parse()
                    .map_err(|e| format!("VTYPE inner_len: {}", e))?;
                let mut inner = Vec::with_capacity(inner_len);
                for k in 0..inner_len {
                    let entry = &parts[6 + k];
                    let mut parts_inner = entry.split('.');
                    let tid: usize = parts_inner
                        .next()
                        .unwrap()
                        .parse()
                        .map_err(|e| format!("VTYPE inner_tid: {}", e))?;
                    let off: usize = parts_inner
                        .next()
                        .unwrap()
                        .parse()
                        .map_err(|e| format!("VTYPE inner_off: {}", e))?;
                    inner.push(EdgeInfo {
                        tile_id: tid,
                        tile_offset: off,
                    });
                }
                let base = 6 + inner_len;
                let ccw_tile_id: usize = parts[base]
                    .parse()
                    .map_err(|e| format!("VTYPE ccw_tid: {}", e))?;
                let ccw_offset: usize = parts[base + 1]
                    .parse()
                    .map_err(|e| format!("VTYPE ccw_off: {}", e))?;
                let kind = parts[base + 2].to_string();
                let cursed: bool = parts[base + 3] != "0";
                vtypes.push(ParsedVtype {
                    id,
                    angle,
                    cw_tile_id,
                    cw_offset,
                    inner,
                    ccw_tile_id,
                    ccw_offset,
                    kind,
                    cursed,
                });
            }
            "WITNESS" => {
                let vtype_id: usize = parts[1].parse().map_err(|e| format!("WITNESS id: {}", e))?;
                let pos: usize = parts[2]
                    .parse()
                    .map_err(|e| format!("WITNESS pos: {}", e))?;
                let n: usize = parts[3].parse().map_err(|e| format!("WITNESS n: {}", e))?;
                let angles: Vec<i8> = parts[4..4 + n]
                    .iter()
                    .map(|s| s.parse())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| format!("WITNESS angles: {}", e))?;
                let edges: Vec<EdgeInfo> = parts[4 + n..4 + 2 * n]
                    .iter()
                    .map(|s| {
                        let mut sp = s.split('.');
                        let tile_id: usize = sp.next().unwrap().parse().unwrap();
                        let tile_offset: usize = sp.next().unwrap().parse().unwrap();
                        EdgeInfo {
                            tile_id,
                            tile_offset,
                        }
                    })
                    .collect();
                let num_cands: usize = parts[4 + 2 * n]
                    .parse()
                    .map_err(|e| format!("WITNESS num_cands: {}", e))?;
                let candidates: Vec<PatchMatch> = parts[4 + 2 * n + 1..4 + 2 * n + 1 + num_cands]
                    .iter()
                    .map(|s| {
                        let mut sp = s.split('.');
                        let start_a: usize = sp.next().unwrap().parse().unwrap();
                        let len: usize = sp.next().unwrap().parse().unwrap();
                        let start_b: usize = sp.next().unwrap().parse().unwrap();
                        let tile_id: usize = sp.next().unwrap().parse().unwrap();
                        PatchMatch {
                            start_a,
                            len,
                            start_b,
                            tile_id,
                        }
                    })
                    .collect();
                let inner_base = 4 + 2 * n + 1 + num_cands;
                let inner_chains: Vec<Vec<EdgeInfo>> = if inner_base < parts.len() {
                    parts[inner_base..]
                        .iter()
                        .map(|s| {
                            if *s == "-" {
                                vec![]
                            } else {
                                s.split(',')
                                    .map(|e| {
                                        let mut sp = e.split('.');
                                        let tid: usize = sp.next().unwrap().parse().unwrap();
                                        let off: usize = sp.next().unwrap().parse().unwrap();
                                        EdgeInfo {
                                            tile_id: tid,
                                            tile_offset: off,
                                        }
                                    })
                                    .collect()
                            }
                        })
                        .collect()
                } else {
                    vec![vec![]; n]
                };
                witnesses.push(ParsedWitness {
                    vtype_id,
                    pos,
                    angles,
                    edges,
                    candidates,
                    inner_chains,
                });
            }
            "TRANS" => {
                let src_id: usize = parts[1].parse().unwrap();
                let dst_id: usize = parts[2].parse().unwrap();
                let _side: &str = parts[3];
                let start_a: usize = parts[4].parse().unwrap();
                let len: usize = parts[5].parse().unwrap();
                let start_b: usize = parts[6].parse().unwrap();
                let tile_id: usize = parts[7].parse().unwrap();
                transitions.push(ParsedTransition {
                    src_id,
                    dst_id,
                    start_a,
                    len,
                    start_b,
                    tile_id,
                });
            }
            _ => {}
        }
    }

    Ok(ParsedFile {
        ring_name,
        tile_seqs,
        vtypes,
        witnesses,
        transitions,
    })
}

fn validate_common<T: IsComplex + IsRingOrField + Units>(
    pf: &ParsedFile,
    tile_ts: &Arc<TileSet<T>>,
) -> Result<(), String> {
    let t0 = Instant::now();

    eprintln!(
        "  Parsed: ring={}, tiles={}, vtypes={}, witnesses={}, transitions={}",
        pf.ring_name,
        pf.tile_seqs.len(),
        pf.vtypes.len(),
        pf.witnesses.len(),
        pf.transitions.len(),
    );

    let num_types = pf.vtypes.len();
    let mut witnesses_by_id: HashMap<usize, &ParsedWitness> = HashMap::new();
    for w in &pf.witnesses {
        witnesses_by_id.insert(w.vtype_id, w);
    }

    eprintln!("  Phase 1: Reconstructing witnesses...");
    let t1 = Instant::now();
    let mi = Arc::new(MatchTypeIndex::new(Arc::clone(tile_ts)));
    let mut reconstructed: HashMap<usize, GrowingPatch<T>> = HashMap::new();
    for pw in &pf.witnesses {
        let n = pw.angles.len();
        let cands = candidates_from_flat(n, pw.candidates.clone());
        let inner_chains = if pw.inner_chains.len() == n {
            pw.inner_chains.clone()
        } else {
            vec![vec![]; n]
        };
        let gp = GrowingPatch::from_parts(
            Arc::clone(&mi),
            pw.angles.clone(),
            pw.edges.clone(),
            cands,
            inner_chains,
        )
        .ok_or_else(|| format!("WITNESS {}: from_parts failed", pw.vtype_id))?;
        reconstructed.insert(pw.vtype_id, gp);
    }
    eprintln!(
        "  Reconstructed {} witnesses in {:.2?}",
        reconstructed.len(),
        t1.elapsed(),
    );

    eprintln!("  Phase 2: Verifying vertex types...");
    let t2 = Instant::now();
    let mut vtype_errors = 0usize;
    for pv in &pf.vtypes {
        let gp = reconstructed
            .get(&pv.id)
            .ok_or_else(|| format!("VTYPE {}: no witness found", pv.id))?;
        let pvt = gp.full_vertex_type_at(witnesses_by_id[&pv.id].pos);
        let expected = VertexType {
            cw: EdgeInfo {
                tile_id: pv.cw_tile_id,
                tile_offset: pv.cw_offset,
            },
            inner: pv.inner.clone(),
            ccw: EdgeInfo {
                tile_id: pv.ccw_tile_id,
                tile_offset: pv.ccw_offset,
            },
        };
        match pvt {
            Some(actual) if actual.cw == expected.cw && actual.ccw == expected.ccw => {}
            Some(actual) => {
                vtype_errors += 1;
                if vtype_errors <= 5 {
                    eprintln!(
                        "  ERROR: VTYPE {} mismatch: expected {:?}, got {:?}",
                        pv.id, expected, actual
                    );
                }
            }
            None => {
                vtype_errors += 1;
                if vtype_errors <= 5 {
                    eprintln!("  ERROR: VTYPE {}: vertex_type_at returned None", pv.id);
                }
            }
        }
    }
    if vtype_errors > 0 {
        return Err(format!("{} vertex type mismatches", vtype_errors));
    }
    eprintln!(
        "  All {} vertex types verified in {:.2?}",
        pf.vtypes.len(),
        t2.elapsed(),
    );

    eprintln!("  Phase 3: Verifying transitions...");
    let t3 = Instant::now();
    let mut trans_errors = 0usize;
    let mut known_ids: BTreeSet<usize> = pf.vtypes.iter().map(|v| v.id).collect();
    for pt in &pf.transitions {
        if !known_ids.contains(&pt.src_id) || !known_ids.contains(&pt.dst_id) {
            trans_errors += 1;
            if trans_errors <= 5 {
                eprintln!("  ERROR: TRANS {} -> {}: unknown id", pt.src_id, pt.dst_id);
            }
            continue;
        }
        let gp = reconstructed
            .get(&pt.src_id)
            .ok_or_else(|| format!("TRANS {} -> {}: no source witness", pt.src_id, pt.dst_id))?;
        let pos = witnesses_by_id[&pt.src_id].pos;
        let old_n = gp.boundary_len();
        let mut gp2 = gp.clone();
        let pm = PatchMatch {
            start_a: pt.start_a,
            len: pt.len,
            start_b: pt.start_b,
            tile_id: pt.tile_id,
        };
        if gp2.add_tile(&pm).is_none() {
            trans_errors += 1;
            if trans_errors <= 5 {
                eprintln!(
                    "  ERROR: TRANS {} -> {}: add_tile failed",
                    pt.src_id, pt.dst_id
                );
            }
            continue;
        }
        let junction_pos = if pt.start_a == pos { old_n - pt.len } else { 0 };
        let junction_vt = gp2.full_vertex_type_at(junction_pos);
        let dst_info = reconstructed.get(&pt.dst_id);
        let dst_pos = witnesses_by_id.get(&pt.dst_id).map(|w| w.pos);
        match (junction_vt, dst_info, dst_pos) {
            (Some(actual), Some(dst_gp), Some(dst_p)) => {
                let expected = dst_gp.full_vertex_type_at(dst_p);
                let matches = match (actual, expected) {
                    (a, Some(e)) => a.cw == e.cw && a.ccw == e.ccw,
                    _ => false,
                };
                if !matches {
                    trans_errors += 1;
                    if trans_errors <= 5 {
                        eprintln!(
                            "  ERROR: TRANS {} -> {}: junction type mismatch",
                            pt.src_id, pt.dst_id
                        );
                    }
                }
            }
            (None, _, _) => {
                trans_errors += 1;
                if trans_errors <= 5 {
                    eprintln!(
                        "  ERROR: TRANS {} -> {}: junction has no vertex type",
                        pt.src_id, pt.dst_id
                    );
                }
            }
            _ => {
                trans_errors += 1;
                if trans_errors <= 5 {
                    eprintln!(
                        "  ERROR: TRANS {} -> {}: missing dst witness",
                        pt.src_id, pt.dst_id
                    );
                }
            }
        }
    }
    if trans_errors > 0 {
        return Err(format!("{} transition errors", trans_errors));
    }
    eprintln!(
        "  All {} transitions verified in {:.2?}",
        pf.transitions.len(),
        t3.elapsed(),
    );

    eprintln!("  Phase 4: Completeness check...");
    let t4 = Instant::now();
    let mut missing_types: BTreeSet<usize> = BTreeSet::new();
    let mut total_match_checks = 0usize;
    for pv in &pf.vtypes {
        if pv.kind == "dead" {
            continue;
        }
        let gp = match reconstructed.get(&pv.id) {
            Some(g) => g,
            None => continue,
        };
        let pos = witnesses_by_id[&pv.id].pos;
        for pm in gp.get_matches_touching_vertex(pos) {
            total_match_checks += 1;
            let old_n = gp.boundary_len();
            let mut gp2 = gp.clone();
            if gp2.add_tile(pm).is_none() || !gp2.is_growing() {
                continue;
            }
            let junction_pos = if pm.start_a == pos { old_n - pm.len } else { 0 };
            if let Some(jvt) = gp2.full_vertex_type_at(junction_pos) {
                let found = pf.vtypes.iter().any(|pv| {
                    let gp = match reconstructed.get(&pv.id) {
                        Some(g) => g,
                        None => return false,
                    };
                    let wp = match witnesses_by_id.get(&pv.id) {
                        Some(w) => w.pos,
                        None => return false,
                    };
                    match gp.full_vertex_type_at(wp) {
                        Some(stored) => stored == jvt,
                        None => false,
                    }
                });
                if !found {
                    missing_types.insert(pv.id);
                }
            }
        }
    }
    if !missing_types.is_empty() {
        return Err(format!(
            "Completeness FAILED: {} vertex types produce unknown junction types: {:?}",
            missing_types.len(),
            missing_types.iter().take(10).collect::<Vec<_>>(),
        ));
    }
    eprintln!(
        "  Completeness: {} match checks passed in {:.2?}",
        total_match_checks,
        t4.elapsed(),
    );

    eprintln!("  Validation PASSED in {:.2?}", t0.elapsed());
    Ok(())
}

fn collect_generic<T: IsComplex + IsRingOrField + Units>(
    ts: Arc<TileSet<T>>,
    ring_name: &str,
    output: &str,
    label: &str,
) {
    eprintln!("[{}] Running BFS...", label);
    let t0 = Instant::now();
    let idx = VertexTypeIndex::new(ts);
    let elapsed = t0.elapsed();

    let n_alive = idx.entries().iter().filter(|e| e.is_alive()).count();
    let n_dead = idx.entries().iter().filter(|e| e.is_dead()).count();
    let n_cursed = idx.entries().iter().filter(|e| e.is_cursed()).count();
    eprintln!(
        "[{}] types={} (alive={}, dead={}, cursed={}) transitions={} time={:.2?}",
        label,
        idx.num_types(),
        n_alive,
        n_dead,
        n_cursed,
        idx.transitions().len(),
        elapsed,
    );

    write_collection(&idx, ring_name, output);
}

fn main() {
    #[cfg(feature = "pprof")]
    let guard = ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
        .unwrap();

    let cli = Cli::parse();

    match &cli.command {
        Commands::Collect { tile, output } => match tile {
            TileSetKind::Penrose => {
                let ts = make_ts_10();
                collect_generic(ts, "ZZ10", output, "penrose");
            }
            _ => {
                let ts = make_ts_12(tile);
                let label = format!("{:?}", tile).to_lowercase();
                collect_generic(ts, "ZZ12", output, &label);
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
            let path = "flamegraph_vtype.svg";
            let mut file = std::fs::File::create(path).unwrap();
            report.flamegraph(&mut file).unwrap();
            eprintln!("Flame graph written to {}", path);
        }
    }
}
