#!/usr/bin/env rust-script
//! Independent fast verifier for a `tilezz-rat-dafsa-blocks` asset.
//!
//! This is the Rust counterpart to the shipped Python tools and the
//! one the project's pipeline runs by default. It decodes every stored
//! rat straight from the gzipped binary blocks -- reimplementing the
//! block format from `schemas/blocks_schema.txt`, with no dependency on
//! the `tilezz` crate that produced the asset -- and then runs the same
//! decode-heavy checks the Python tools do, plus per-block integrity:
//!
//!   * per-block SHA-256 against the manifest (cf. verify_sha256.py),
//!   * every stored sequence is its own dihedral-canonical CCW form
//!     (cf. verify_canonical.py): turn-sum > 0 and equal to the lex-min
//!     over all rotations of the sequence AND of its reverse,
//!   * the seven per-perimeter family series (free, oneSided, achiral,
//!     rotationSymmetric, symmetric, subring, coset) re-derived from the
//!     decoded rats and checked against the RO-Crate `variableMeasured`
//!     and `n_sequences` (cf. count.py --verify).
//!
//! It is orders of magnitude faster than the Python sweep (compiled +
//! multi-threaded), which is why it is the default pipeline check. The
//! Python tools stay in `tools/` as the zero-dependency, no-toolchain
//! fallback; CI asserts the two agree.
//!
//! Run (needs a Rust toolchain + rust-script):
//!
//!     cargo install rust-script
//!     rust-script tools/verify.rs [path/to/asset]      # default: .
//!
//! Exits 0 if every check passes, 1 on any mismatch / violation, 2 on
//! usage or I/O error.
//!
//! ```cargo
//! [dependencies]
//! flate2 = "1"
//! sha2 = "0.10"
//! serde_json = "1"
//! ```

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use flate2::read::GzDecoder;
use serde_json::Value;
use sha2::{Digest, Sha256};

// ---- binary layout (schemas/blocks_schema.txt); all LE ----
const BLOCK_MAGIC: &[u8; 4] = b"TRB1";
const HEADER_BYTES: usize = 16;
const STATE_REC_BYTES: usize = 16;
const EDGE_REC_BYTES: usize = 8;

const FAMILIES: [&str; 7] = [
    "free",
    "oneSided",
    "achiral",
    "rotationSymmetric",
    "symmetric",
    "subring",
    "coset",
];

fn rd_u32(b: &[u8], o: usize) -> u32 {
    u32::from_le_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]])
}

/// Global, flat CSR view of the whole DAFSA (root + every block).
/// State 0 is the manifest root; states 1..n_states come from blocks.
struct Dafsa {
    /// edges[edge_off[s]..edge_off[s+1]] are state `s`'s out-edges.
    edge_off: Vec<u64>,
    labels: Vec<i8>,
    targets: Vec<u32>,
    is_accept: Vec<bool>,
}

/// Index of the lexicographically least rotation of the length-`n`
/// sequence whose doubled form (`seq ++ seq`) is `dd`. Comparing
/// rotations as `dd[k..k+n]` slices avoids the per-element modulo that
/// dominates the cost at billions of rats. O(n^2), allocation-free.
fn least_rot_doubled(dd: &[i8], n: usize) -> usize {
    let mut best = 0usize;
    for k in 1..n {
        if dd[k..k + n] < dd[best..best + n] {
            best = k;
        }
    }
    best
}

/// Per-perimeter family tallies plus canonical-check bookkeeping.
struct Acc {
    free: Vec<u64>,
    achiral: Vec<u64>,
    rotsym: Vec<u64>,
    symm: Vec<u64>,
    subring: Vec<u64>,
    coset: Vec<u64>,
    checked: u64,
    violations: u64,
    examples: Vec<String>,
}
impl Acc {
    fn new(maxp: usize) -> Self {
        let z = || vec![0u64; maxp + 1];
        Acc {
            free: z(),
            achiral: z(),
            rotsym: z(),
            symm: z(),
            subring: z(),
            coset: z(),
            checked: 0,
            violations: 0,
            examples: Vec::new(),
        }
    }
    fn merge(&mut self, o: &Acc) {
        for i in 0..self.free.len() {
            self.free[i] += o.free[i];
            self.achiral[i] += o.achiral[i];
            self.rotsym[i] += o.rotsym[i];
            self.symm[i] += o.symm[i];
            self.subring[i] += o.subring[i];
            self.coset[i] += o.coset[i];
        }
        self.checked += o.checked;
        self.violations += o.violations;
        for e in &o.examples {
            if self.examples.len() < 10 {
                self.examples.push(e.clone());
            }
        }
    }

    /// Classify one decoded rat (no length prefix) into the families and
    /// run the canonical-CCW check. Mirrors count.py / verify_canonical.py
    /// exactly but allocation-free.
    fn record(&mut self, rat: &[i8]) {
        let n = rat.len();
        if n == 0 || n >= self.free.len() {
            self.violations += 1;
            if self.examples.len() < 10 {
                self.examples.push(format!("{rat:?} -- length {n} out of range"));
            }
            return;
        }
        self.free[n] += 1;
        if rat.iter().all(|&a| a % 2 == 0) {
            self.subring[n] += 1;
        }
        if rat.iter().all(|&a| a % 2 != 0) {
            self.coset[n] += 1;
        }

        // Doubled stack buffers (no heap, no modulo): dd = rat ++ rat,
        // rr = reverse(rat) ++ reverse(rat). Every rotation is then a
        // contiguous length-n slice. n <= 31 here (2n <= 64).
        let mut dd = [0i8; 64];
        let mut rr = [0i8; 64];
        for i in 0..n {
            dd[i] = rat[i];
            dd[n + i] = rat[i];
            rr[i] = rat[n - 1 - i];
            rr[n + i] = rat[n - 1 - i];
        }
        let ia = least_rot_doubled(&dd, n);
        let ib = least_rot_doubled(&rr, n);
        // achiral: least rotation of rat == least rotation of its reverse
        let achiral = dd[ia..ia + n] == rr[ib..ib + n];
        if achiral {
            self.achiral[n] += 1;
        }
        // rotational symmetry: some non-zero rotation reproduces rat
        let rot = (1..n).any(|d| dd[d..d + n] == dd[0..n]);
        if rot {
            self.rotsym[n] += 1;
        }
        if achiral || rot {
            self.symm[n] += 1;
        }

        // canonical CCW: turn-sum > 0, and rat equals the min over
        // rotations of itself and its reverse. The stored rat is
        // rotation-0 of itself, so it is canonical iff (a) it is its own
        // least rotation and (b) that least rotation is <= the reverse's
        // least rotation.
        let sum: i64 = rat.iter().map(|&a| a as i64).sum();
        let is_own_min = dd[0..n] == dd[ia..ia + n];
        let le_reverse = dd[ia..ia + n] <= rr[ib..ib + n];
        if !(sum > 0 && is_own_min && le_reverse) {
            self.violations += 1;
            if self.examples.len() < 10 {
                let reason = if sum <= 0 {
                    format!("not CCW (turn-sum {sum} <= 0)")
                } else {
                    "not the dihedral-canonical (lex-min rotation/reverse) form".to_string()
                };
                self.examples.push(format!("{rat:?} -- {reason}"));
            }
        }
        self.checked += 1;
    }
}

/// DFS from `start` (already reached via `prefix`, which includes the
/// leading length byte), recording every accepted rat into `acc`. The
/// rat is the path labels with the length byte stripped.
fn dfs(d: &Dafsa, start: u32, mut prefix: Vec<i8>, acc: &mut Acc) {
    if d.is_accept[start as usize] {
        acc.record(&prefix[1..]);
    }
    let mut stack: Vec<(u32, u64)> = vec![(start, d.edge_off[start as usize])];
    while let Some(&mut (state, ref mut cur)) = stack.last_mut() {
        let end = d.edge_off[state as usize + 1];
        if *cur < end {
            let e = *cur as usize;
            *cur += 1;
            let label = d.labels[e];
            let tgt = d.targets[e];
            prefix.push(label);
            if d.is_accept[tgt as usize] {
                acc.record(&prefix[1..]);
            }
            stack.push((tgt, d.edge_off[tgt as usize]));
        } else {
            stack.pop();
            // The bottom (start) frame's label lives in the seed prefix,
            // not pushed by this loop -- don't pop it.
            if !stack.is_empty() {
                prefix.pop();
            }
        }
    }
}

/// Cannot perform the verification (missing / unreadable / unparseable
/// inputs, bad usage): exit 2, matching the Python tools' convention.
fn die(msg: impl AsRef<str>) -> ! {
    eprintln!("verify.rs: {}", msg.as_ref());
    std::process::exit(2);
}

/// The verification ran and found the asset invalid (block hash mismatch
/// or a structural inconsistency): exit 1, distinct from the exit-2
/// "couldn't run" so callers can tell a failed dataset from a failed
/// invocation.
fn fail(msg: impl AsRef<str>) -> ! {
    eprintln!("verify.rs: VERIFICATION FAILED: {}", msg.as_ref());
    std::process::exit(1);
}

fn main() {
    let arg = std::env::args().nth(1).unwrap_or_else(|| ".".to_string());
    let dir = PathBuf::from(arg);
    let manifest_path = dir.join("block_index.json");
    let crate_path = dir.join("ro-crate-metadata.json");
    if !manifest_path.exists() || !crate_path.exists() {
        die(format!(
            "need block_index.json + ro-crate-metadata.json in {}",
            dir.display()
        ));
    }

    let manifest = read_json(&manifest_path);
    let n_sequences = manifest["n_sequences"]
        .as_u64()
        .unwrap_or_else(|| die("n_sequences"));
    let maxp = manifest["max_indexed_length"]
        .as_u64()
        .unwrap_or_else(|| die("max_indexed_length")) as usize;

    // Decode + integrity-check the whole DAFSA (root from the manifest,
    // states from the gzipped blocks).
    let dafsa = load_dafsa(&dir, &manifest);

    // Split the tree into enough independent subtrees to balance all
    // cores, re-derive the families + canonical-check every rat in
    // parallel, then fold in the accepts found at interior frontier nodes.
    let hw = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let (tasks, interior) = build_frontier(&dafsa, maxp, hw);
    let mut acc = run_workers(&dafsa, &tasks, hw, maxp);
    acc.merge(&interior);

    let emitted = load_emitted(&crate_path);
    if !report(&acc, &emitted, n_sequences, maxp) {
        std::process::exit(1);
    }
}

/// Read + parse a JSON file, or exit 2 (cannot run).
fn read_json(path: &Path) -> Value {
    let bytes =
        std::fs::read(path).unwrap_or_else(|e| die(format!("read {}: {e}", path.display())));
    serde_json::from_slice(&bytes).unwrap_or_else(|e| die(format!("parse {}: {e}", path.display())))
}

/// Build the global CSR view of the whole DAFSA: state 0 (root) from the
/// manifest, states 1..n_states from the gzipped blocks. Integrity-checks
/// each block's SHA-256 against the manifest and the structural counts;
/// exits via `fail` (exit 1) on any invalid asset.
fn load_dafsa(dir: &Path, manifest: &Value) -> Dafsa {
    let n_states = manifest["n_states"].as_u64().unwrap_or_else(|| die("n_states")) as usize;
    let n_edges = manifest["n_edges"].as_u64().unwrap_or_else(|| die("n_edges")) as usize;

    let mut edge_off = vec![0u64; n_states + 1];
    let mut labels: Vec<i8> = Vec::with_capacity(n_edges);
    let mut targets: Vec<u32> = Vec::with_capacity(n_edges);
    let mut is_accept = vec![false; n_states];
    let mut running: u64 = 0;

    // state 0 = root, served from the manifest
    let root = &manifest["root"];
    is_accept[0] = root["is_accept"].as_bool().unwrap_or(false);
    for e in root["edges"].as_array().unwrap_or(&vec![]).iter() {
        labels.push(e["label"].as_i64().unwrap_or_else(|| die("root label")) as i8);
        targets.push(e["target"].as_u64().unwrap_or_else(|| die("root target")) as u32);
        running += 1;
    }

    // states 1..n_states from the content-addressed blocks, in id order
    let mut blocks: Vec<(u32, String)> = manifest["blocks"]
        .as_array()
        .unwrap_or_else(|| die("blocks"))
        .iter()
        .map(|b| {
            (
                b["first_state"]
                    .as_u64()
                    .unwrap_or_else(|| die("first_state")) as u32,
                b["sha256"]
                    .as_str()
                    .unwrap_or_else(|| die("sha256"))
                    .to_string(),
            )
        })
        .collect();
    blocks.sort_by_key(|(fs, _)| *fs);

    let mut next_gid: u32 = 1;
    for (first_state, sha_want) in &blocks {
        if *first_state != next_gid {
            fail(format!(
                "block first_state {first_state} != expected {next_gid} (gap/overlap)"
            ));
        }
        let path = dir.join("blocks").join(format!("{sha_want}.bin"));
        let gz =
            std::fs::read(&path).unwrap_or_else(|e| die(format!("read {}: {e}", path.display())));
        // integrity: sha256 of the gzipped bytes == manifest entry
        let got = {
            let mut h = Sha256::new();
            h.update(&gz);
            hex(&h.finalize())
        };
        if &got != sha_want {
            fail(format!(
                "block sha256 mismatch (file {got}, manifest {sha_want})"
            ));
        }
        let mut buf = Vec::new();
        std::io::Read::read_to_end(&mut GzDecoder::new(&gz[..]), &mut buf)
            .unwrap_or_else(|e| die(format!("gunzip {}: {e}", path.display())));
        decode_block(
            &buf,
            *first_state,
            &mut edge_off,
            &mut labels,
            &mut targets,
            &mut is_accept,
            &mut running,
        );
        next_gid += block_state_count(&buf);
    }
    edge_off[n_states] = running;
    if next_gid as usize != n_states {
        fail(format!(
            "blocks cover {next_gid} states but manifest says n_states={n_states}"
        ));
    }
    if labels.len() != n_edges {
        fail(format!(
            "decoded {} edges but manifest n_edges={n_edges}",
            labels.len()
        ));
    }
    Dafsa {
        edge_off,
        labels,
        targets,
        is_accept,
    }
}

/// Split the DAFSA into independent subtrees for parallel processing.
/// BFS-expand from the root until there are at least `8 * hw` task roots
/// (or a depth cap is hit), so the work balances even when one length
/// class holds most of the mass under a few shallow edges (ZZ8 n=20: most
/// of a(10) sits under one or two first turns -- a fixed shallow split
/// would starve cores). Returns the surviving frontier as task roots plus
/// an `Acc` holding accepts found at the interior nodes we expanded past.
/// Every rat is counted exactly once: an accept at an expanded node is
/// recorded here; all others are recorded by the worker that owns the
/// task root at/above them.
fn build_frontier(d: &Dafsa, maxp: usize, hw: usize) -> (Vec<(u32, Vec<i8>)>, Acc) {
    let target = (hw * 8).max(64);
    let edges_of = |s: u32| -> std::ops::Range<usize> {
        d.edge_off[s as usize] as usize..d.edge_off[s as usize + 1] as usize
    };
    let mut interior = Acc::new(maxp);
    let mut frontier: Vec<(u32, Vec<i8>)> = edges_of(0)
        .map(|i| (d.targets[i], vec![d.labels[i]]))
        .collect();
    let mut depth = 1usize;
    while frontier.len() < target && depth < 12 {
        let mut next_front: Vec<(u32, Vec<i8>)> = Vec::with_capacity(frontier.len() * 2);
        let mut grew = false;
        for (state, prefix) in frontier.drain(..) {
            if d.is_accept[state as usize] {
                interior.record(&prefix[1..]); // interior accept: count now
            }
            let r = edges_of(state);
            if r.is_empty() {
                continue; // terminal node, already accounted above
            }
            for i in r {
                let mut p = prefix.clone();
                p.push(d.labels[i]);
                next_front.push((d.targets[i], p));
            }
            grew = true;
        }
        frontier = next_front;
        depth += 1;
        if !grew {
            break;
        }
    }
    (frontier, interior)
}

/// Run every task root through `dfs` on a pool of up to `hw` threads,
/// merging the per-thread tallies. Tasks are pulled from a shared atomic
/// cursor, so uneven subtree sizes still balance across the cores.
fn run_workers(d: &Dafsa, tasks: &[(u32, Vec<i8>)], hw: usize, maxp: usize) -> Acc {
    let n_threads = hw.min(tasks.len().max(1));
    let next = std::sync::atomic::AtomicUsize::new(0);
    let merged = Mutex::new(Acc::new(maxp));
    std::thread::scope(|scope| {
        for _ in 0..n_threads {
            scope.spawn(|| {
                let mut local = Acc::new(maxp);
                loop {
                    let idx = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if idx >= tasks.len() {
                        break;
                    }
                    let (s, ref prefix) = tasks[idx];
                    dfs(d, s, prefix.clone(), &mut local);
                }
                merged.lock().unwrap().merge(&local);
            });
        }
    });
    merged.into_inner().unwrap()
}

/// Print the seven re-derived family series, compare each against the
/// emitted RO-Crate metadata, and print the canonical-check + totals
/// summary. Returns true iff everything matched: all seven series equal
/// metadata, `free` total == `n_sequences`, and no canonical violations.
fn report(acc: &Acc, emitted: &HashMap<String, Vec<i64>>, n_sequences: u64, maxp: usize) -> bool {
    // oneSided = 2*free - achiral; every series computed here.
    let series = |f: &str| -> Vec<i64> {
        let pick = |v: &Vec<u64>| (1..=maxp).map(|n| v[n] as i64).collect::<Vec<i64>>();
        match f {
            "free" => pick(&acc.free),
            "achiral" => pick(&acc.achiral),
            "rotationSymmetric" => pick(&acc.rotsym),
            "symmetric" => pick(&acc.symm),
            "subring" => pick(&acc.subring),
            "coset" => pick(&acc.coset),
            "oneSided" => (1..=maxp)
                .map(|n| 2 * acc.free[n] as i64 - acc.achiral[n] as i64)
                .collect(),
            _ => unreachable!(),
        }
    };

    let mut bad = 0;
    for f in FAMILIES {
        let got = series(f);
        let ok = emitted.get(f).map(|e| e == &got).unwrap_or(false);
        if !ok {
            bad += 1;
        }
        let body = got
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "{f}: {body}{}",
            if ok { "" } else { "   *** MISMATCH vs metadata ***" }
        );
    }

    let free_total: u64 = acc.free.iter().sum();
    println!(
        "# n_sequences {n_sequences}; free total {free_total}; checked {}; mode rust-verify (all re-derived from rats)",
        acc.checked
    );

    let total_ok = free_total == n_sequences;
    if acc.violations > 0 {
        println!("# canonical: {} violation(s):", acc.violations);
        for e in acc.examples.iter().take(10) {
            println!("#   NON-CANONICAL: {e}");
        }
    } else {
        println!(
            "# canonical: OK ({} rat(s) verified canonical CCW)",
            acc.checked
        );
    }

    if bad > 0 || !total_ok || acc.violations > 0 {
        println!(
            "# FAIL: {bad} family series disagree with metadata{}{}",
            if total_ok { "" } else { "; free total != n_sequences" },
            if acc.violations == 0 {
                String::new()
            } else {
                format!("; {} canonical violation(s)", acc.violations)
            }
        );
        return false;
    }
    println!(
        "# OK ({} families re-derived + verified, all canonical CCW, to perimeter {maxp})",
        FAMILIES.len()
    );
    true
}

fn hex(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

/// Number of states a decoded block body holds (its header `n_states`).
fn block_state_count(buf: &[u8]) -> u32 {
    rd_u32(buf, 8)
}

/// Decode a block body into the growing global CSR. `running` is the
/// global edge cursor; states are appended in id order.
#[allow(clippy::too_many_arguments)]
fn decode_block(
    buf: &[u8],
    first_state: u32,
    edge_off: &mut [u64],
    labels: &mut Vec<i8>,
    targets: &mut Vec<u32>,
    is_accept: &mut [bool],
    running: &mut u64,
) {
    if buf.len() < HEADER_BYTES || &buf[0..4] != BLOCK_MAGIC {
        fail("block: bad magic / short header");
    }
    let fs = rd_u32(buf, 4);
    if fs != first_state {
        fail(format!("block header first_state {fs} != manifest {first_state}"));
    }
    let nb_states = rd_u32(buf, 8) as usize;
    let nb_edges = rd_u32(buf, 12) as usize;
    let states_off = HEADER_BYTES;
    let edges_off = states_off + nb_states * STATE_REC_BYTES;
    if buf.len() < edges_off + nb_edges * EDGE_REC_BYTES {
        fail("block: truncated body");
    }
    // local edge offsets per state
    let local_eo = |j: usize| -> usize {
        if j < nb_states {
            rd_u32(buf, states_off + j * STATE_REC_BYTES) as usize
        } else {
            nb_edges
        }
    };
    for j in 0..nb_states {
        let gid = first_state as usize + j;
        let o = states_off + j * STATE_REC_BYTES;
        // edges_offset @ o (u32), count @ o+4 (u64), is_accept @ o+12 (u8)
        is_accept[gid] = buf[o + 12] != 0;
        edge_off[gid] = *running;
        let lo = local_eo(j);
        let hi = local_eo(j + 1);
        for e in lo..hi {
            let eo = edges_off + e * EDGE_REC_BYTES;
            labels.push(buf[eo] as i8);
            targets.push(rd_u32(buf, eo + 4));
            *running += 1;
        }
    }
}

/// Parse the seven `variableMeasured` series out of the RO-Crate.
fn load_emitted(crate_path: &Path) -> HashMap<String, Vec<i64>> {
    let v: Value = serde_json::from_slice(
        &std::fs::read(crate_path).unwrap_or_else(|e| die(format!("read crate: {e}"))),
    )
    .unwrap_or_else(|e| die(format!("parse crate: {e}")));
    let mut out = HashMap::new();
    if let Some(graph) = v["@graph"].as_array() {
        for e in graph {
            let is_dataset = e["@id"] == Value::from("./")
                || e["@type"] == Value::from("Dataset")
                || e["@type"]
                    .as_array()
                    .map(|a| a.iter().any(|t| t == &Value::from("Dataset")))
                    .unwrap_or(false);
            if !is_dataset {
                continue;
            }
            if let Some(vm) = e["variableMeasured"].as_array() {
                for pv in vm {
                    if let (Some(name), Some(val)) = (pv["name"].as_str(), pv["value"].as_str()) {
                        let series: Vec<i64> =
                            val.split(',').filter_map(|x| x.trim().parse().ok()).collect();
                        out.insert(name.to_string(), series);
                    }
                }
            }
        }
    }
    out
}
