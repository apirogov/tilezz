//! [RO-Crate 1.2](https://www.researchobject.org/ro-crate/specification/1.2/)
//! metadata emitter for a `tilezz-rat-dafsa-blocks` asset directory.
//!
//! Pipeline: after [`RatDafsa::write_blocks`] populates `dir/` with
//! `block_index.json` and `blocks/<sha256>.bin`, call
//! [`write_archival_extras`] to add the self-description files, then
//! [`write_ro_crate`] to write the RO-Crate 1.2 entry-point JSON.
//!
//! The full archival-extras set:
//!
//! - `schemas/` (`block_index.schema.json`, `blocks_schema.txt`,
//!   `rat_schema.txt`) -- the prose and machine-checkable schemas
//!   linked from the manifest's `conformsTo` chain.
//! - `tools/decode.py` -- standalone Python 3 decoder (no external
//!   deps) that walks the blocked DAFSA and prints every rat as a
//!   line of plain text. Lets an archivist extract sequences
//!   without having to compile Rust.
//! - `tools/verify_sha256.py` -- standalone Python 3 script that
//!   checks every sha256 recorded in `ro-crate-metadata.json`
//!   against the on-disk file bytes. Same dependency-free shape
//!   as decode.py.
//! - `README.md` -- human-readable entry point: dataset summary,
//!   author and copyright notice (required by CC-BY-SA's
//!   attribution clause), the CC-BY-SA 4.0 license terms, an exact
//!   rebuild recipe (repo URL, commit, CLI invocation), and a
//!   pointer at the verify script.
//! - `ro-crate-metadata.json` -- the RO-Crate 1.2 entry point. The
//!   single starting point an archivist needs: it links every file
//!   above with sha256 / conformsTo / encodingFormat properties, and
//!   the `CreateAction` records the verbatim CLI invocation so the
//!   asset can be regenerated bit-for-bit.
//!
//! The crate stays a leaf consumer of `serde_json::Value` -- no
//! `ro-crate-rs` dependency. The reasoning: that crate pulls
//! `reqwest`/`tokio`/`chrono`/`uuid`/`url`/`zip`/`walkdir` for read
//! / fetch / validate paths we don't need. A 200-line writer over
//! `serde_json::Value` keeps the dep surface flat.
//!
//! Reproducibility: every output (manifest, schemas, blocks) is
//! deterministic given the input rats. The one timestamp we have
//! to choose -- `CreateAction.endTime` -- honours
//! [`SOURCE_DATE_EPOCH`](https://reproducible-builds.org/specs/source-date-epoch/)
//! if set, so two builds from the same commit + same enumeration
//! parameters produce bit-identical `ro-crate-metadata.json`.

use std::collections::BTreeMap;
use std::io;
use std::path::Path;

use serde_json::{Value, json};
use sha2::{Digest, Sha256};

/// The three prose / JSON schema files we ship inside every
/// self-describing asset. `include_str!`-baked at compile time so
/// the build doesn't depend on the source tree.
const BLOCKS_SCHEMA_TXT: &str = include_str!("schemas/blocks_schema.txt");
const RAT_SCHEMA_TXT: &str = include_str!("schemas/rat_schema.txt");
const BLOCK_INDEX_SCHEMA_JSON: &str = include_str!("schemas/block_index.schema.json");

/// Standalone Python decoder (no external deps) for the blocked
/// DAFSA format. Shipped at `tools/decode.py` inside every asset
/// so a future archivist can extract sequences to plain text
/// without compiling Rust.
const DECODER_PY: &str = include_str!("extras/decode.py");

/// Standalone Python SHA-256 verifier. Shipped at
/// `tools/verify_sha256.py` inside every asset so an archivist can
/// check the recorded hashes match on-disk bytes without copy-
/// pasting a heredoc out of the README. Exits 0 on full match, 1
/// on mismatch.
const VERIFY_SHA256_PY: &str = include_str!("extras/verify_sha256.py");

/// Standalone Python per-perimeter family counter + verifier (one
/// tool, two modes). Shipped at `tools/count.py` inside every asset.
/// `--print` (default, fast) reads the `free` series off the DAFSA
/// rank index and cross-checks it against the RO-Crate, echoing the
/// other six families from `variableMeasured` for display. `--verify`
/// (slow) decodes every rat, independently re-derives all seven
/// families, and checks each against the RO-Crate -- the
/// metadata-integrity gate. Merges the former count_by_length.py +
/// verify_counts.py.
const COUNT_PY: &str = include_str!("extras/count.py");

/// Standalone Python canonical-form verifier. Shipped at
/// `tools/verify_canonical.py` inside every asset. Independently
/// recomputes each stored rat's dihedral-canonical CCW form (lex-min
/// over rotations of the sequence and its reverse, turn-sum > 0) and
/// checks the stored sequence already equals it -- catching any
/// non-canonical entry and, since the DAFSA stores each string once,
/// excluding the "non-canonical duplicate" failure mode that a hash /
/// count check cannot see.
const VERIFY_CANONICAL_PY: &str = include_str!("extras/verify_canonical.py");

/// How the dataset was produced. Determines what `reproduce.sh`
/// emits (one-step in-memory build vs. three-stage streaming
/// pipeline) and what the `CreateAction.description` records.
#[derive(Debug, Clone, Copy)]
pub enum ProducedVia {
    /// Single-process, in-memory build via `rat_enum --mode
    /// dafsa-blocks`. Suitable for moderately-sized assets that
    /// fit the HashSet-of-Vec<i8> peak working set.
    InMemory,
    /// Three-stage streaming pipeline (`--mode stream` ->
    /// `--mode merge` -> `--mode build`). Bounded memory; needed
    /// for assets too large to materialise in one HashSet (e.g.
    /// ZZ12 perim >= 14, ZZ10 perim >= 15). Records all three
    /// commands in the reproduce script so an archivist can
    /// rebuild on a memory-constrained host.
    StreamingPipeline,
}

/// Parameters describing one enumeration run, used to fill in the
/// human-readable Dataset properties (`name`, `description`,
/// `keywords`) and any external cross-references (OEIS sequence
/// shorthand).
#[derive(Debug, Clone)]
pub struct AssetParams<'a> {
    /// Cyclotomic ring index `n` in ZZn.
    pub ring: u8,
    /// DFS depth bound (maximum polygon perimeter).
    pub max_steps: usize,
    /// DFS direction step (1 = full ring; >1 = subring slice).
    pub step: i8,
    /// `true` if the canonicalization was free (= full dihedral
    /// symmetry reduction), `false` for plain rotation-canonical.
    pub free: bool,
    /// `--target-block-bytes` used when writing the blocks.
    /// Structurally observable (changes how many block files appear
    /// and which states each contains), so a faithful reproduction
    /// must use the same value -- echoed back in the reproduce
    /// script.
    pub target_block_bytes: u32,
    /// Final rat count, surfaced into the Dataset `description`.
    pub n_sequences: u64,
    /// Optional OEIS reference (e.g. `"A316192"`). Becomes a
    /// `subjectOf` contextual entity in the graph so an archivist
    /// can pivot from the asset to the underlying sequence.
    pub oeis_a_number: Option<&'a str>,
    /// Which CLI mode (or pipeline of modes) produced the asset.
    /// Selects the reproduction recipe in `reproduce.sh` and the
    /// CreateAction.description recorded in `ro-crate-metadata.json`.
    pub produced_via: ProducedVia,
}

/// Compile-time package metadata, all sourced from `Cargo.toml`
/// (no hard-coded duplication).
struct PkgMeta {
    name: &'static str,
    version: &'static str,
    repository: &'static str,
    description: &'static str,
    /// Colon-separated `Cargo.toml` authors field, e.g.
    /// `"Anton Pirogov <a@b>:Someone Else <c@d>"`. Use
    /// [`format_authors`] to render for display.
    authors: &'static str,
    commit: &'static str,
}

/// ORCID iD per author, keyed by the email in CARGO_PKG_AUTHORS (which
/// carries no ORCID). When present, the author's RO-Crate Person entity
/// uses the ORCID URI as its @id -- the resolvable, globally-unique form
/// -- and `creator` references it by that @id.
const AUTHOR_ORCIDS: &[(&str, &str)] = &[(
    "apirogov@users.noreply.github.com",
    "https://orcid.org/0000-0002-5077-7497",
)];

/// Look up the ORCID URI for an author email, if known. Returns
/// `None` when the email is absent or not in [`AUTHOR_ORCIDS`].
fn orcid_for(email: Option<&str>) -> Option<&'static str> {
    let email = email?;
    AUTHOR_ORCIDS
        .iter()
        .find(|(e, _)| *e == email)
        .map(|(_, orcid)| *orcid)
}

const PKG: PkgMeta = PkgMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    repository: env!("CARGO_PKG_REPOSITORY"),
    description: env!("CARGO_PKG_DESCRIPTION"),
    authors: env!("CARGO_PKG_AUTHORS"),
    // Provided by `build.rs`; "unknown" outside a git checkout.
    commit: env!("TILEZZ_GIT_COMMIT"),
};

/// RFC-3339 timestamp string for `CreateAction.endTime`. Honours
/// [`SOURCE_DATE_EPOCH`] when set; otherwise records the current
/// system time. Both branches use UTC and second precision (RO-Crate
/// 1.2's CreateAction example uses sub-day precision, so seconds is
/// well within spec).
fn build_endtime() -> String {
    let secs: i64 = match std::env::var("SOURCE_DATE_EPOCH")
        .ok()
        .and_then(|s| s.parse().ok())
    {
        Some(s) => s,
        None => std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs() as i64)
            .unwrap_or(0),
    };
    format_iso_utc(secs)
}

/// Calendar date `YYYY-MM-DD` of the build, used as the dataset
/// `version`. The leading 10 chars of [`build_endtime`]'s RFC-3339
/// string, so it honours `SOURCE_DATE_EPOCH` the same way.
fn build_date() -> String {
    build_endtime()[..10].to_string()
}

/// Format a POSIX-epoch second count as `YYYY-MM-DDTHH:MM:SSZ`. No
/// external chrono dep: the formula is a stable
/// civil-from-days algorithm (Hinnant 2013). Tested against fixed
/// epoch values in this module's unit tests.
fn format_iso_utc(secs: i64) -> String {
    let days = secs.div_euclid(86_400);
    let tod = secs.rem_euclid(86_400);
    let hour = (tod / 3600) as u32;
    let minute = ((tod % 3600) / 60) as u32;
    let second = (tod % 60) as u32;

    // Civil-from-days (Howard Hinnant). Days since 1970-01-01.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = (yoe as i64) + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = if month <= 2 { y + 1 } else { y };

    format!("{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}Z")
}

/// Streaming-SHA-256 of a file, returned as lowercase hex. The
/// `sha2` crate processes the buffer in fixed chunks, so even
/// multi-GB block files stay within a small constant memory budget.
fn sha256_hex(path: &Path) -> io::Result<String> {
    let mut file = std::fs::File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = io::Read::read(&mut file, &mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(super::hex_lower(&hasher.finalize()))
}

/// One file's RO-Crate File-entity properties.
struct FileMeta {
    /// Relative path inside the crate root. Becomes the JSON-LD `@id`.
    rel_path: String,
    /// Bytes on disk; becomes `contentSize`.
    size: u64,
    /// Lowercase hex SHA-256; becomes the schema.org `sha256` property.
    sha256: String,
}

/// Snapshot the asset directory: every file under `dir` (excluding
/// the in-progress ro-crate-metadata.json itself, which doesn't
/// exist yet) is hashed and its relative path recorded. Result is
/// sorted lex by path for reproducibility.
fn snapshot_dir(dir: &Path) -> io::Result<Vec<FileMeta>> {
    let mut out: Vec<FileMeta> = Vec::new();
    let mut stack: Vec<std::path::PathBuf> = vec![dir.to_path_buf()];
    while let Some(d) = stack.pop() {
        for entry in std::fs::read_dir(&d)? {
            let entry = entry?;
            let path = entry.path();
            let file_type = entry.file_type()?;
            if file_type.is_dir() {
                stack.push(path);
                continue;
            }
            let rel = path
                .strip_prefix(dir)
                .map_err(|e| io::Error::other(e.to_string()))?;
            // Skip the metadata file itself (this avoids a
            // chicken-and-egg on re-runs that pass over an existing
            // bundle) and any platform .DS_Store / .gitkeep noise.
            let name = rel.file_name().and_then(|s| s.to_str()).unwrap_or("");
            if name == "ro-crate-metadata.json" || name == ".DS_Store" || name == ".gitkeep" {
                continue;
            }
            let meta = entry.metadata()?;
            let sha256 = sha256_hex(&path)?;
            out.push(FileMeta {
                rel_path: rel
                    .to_str()
                    .ok_or_else(|| io::Error::other("non-UTF-8 path"))?
                    // Cross-platform: store forward slashes in the JSON-LD
                    // even on Windows, since the same metadata may be served
                    // over HTTP from a POSIX host.
                    .replace('\\', "/"),
                size: meta.len(),
                sha256,
            });
        }
    }
    out.sort_by(|a, b| a.rel_path.cmp(&b.rel_path));
    Ok(out)
}

/// Guess a MIME type from a file extension. Limited to the
/// extensions our writer actually produces; everything else falls
/// back to `application/octet-stream`.
fn encoding_format_of(rel_path: &str) -> &'static str {
    if rel_path.ends_with(".json") {
        "application/json"
    } else if rel_path.ends_with(".bin") {
        // Block files are gzip-encoded; declare the wire encoding,
        // not the payload type. RO-Crate consumers can chain on
        // `conformsTo` to find the inner format.
        "application/gzip"
    } else if rel_path.ends_with(".txt") {
        "text/plain"
    } else if rel_path.ends_with(".md") {
        "text/markdown"
    } else if rel_path.ends_with(".py") {
        "text/x-python"
    } else if rel_path.ends_with(".sh") {
        "application/x-sh"
    } else {
        "application/octet-stream"
    }
}

/// Pick the `conformsTo` chain for a file in the asset directory.
/// The block_index.json conforms to BOTH the formal JSON Schema
/// AND the prose schema; block files conform to the prose schema
/// (the JSON Schema can't cover binary). Schema files themselves
/// don't claim conformance to anything inside the crate.
fn conforms_to_of(rel_path: &str) -> Vec<Value> {
    if rel_path == "block_index.json" {
        vec![
            json!({"@id": "schemas/block_index.schema.json"}),
            json!({"@id": "schemas/blocks_schema.txt"}),
        ]
    } else if rel_path.starts_with("blocks/") && rel_path.ends_with(".bin") {
        vec![json!({"@id": "schemas/blocks_schema.txt"})]
    } else {
        vec![]
    }
}

/// Build a short, archivist-friendly description of one file. Used
/// for both `File` entities' `name`/`description` and the schema
/// entities' `name`/`description` (since RO-Crate 1.2's File entity
/// SHOULD have a `name`).
fn human_label_for(rel_path: &str) -> (&'static str, &'static str) {
    match rel_path {
        "block_index.json" => (
            "block_index.json",
            "Manifest for the tilezz-rat-dafsa-blocks asset (state / edge counts, root state record, content-addressed block index).",
        ),
        "schemas/block_index.schema.json" => (
            "block_index.schema.json (JSON Schema)",
            "Formal JSON Schema (draft 2020-12) for block_index.json. Machine-validatable; covers the manifest JSON only -- see blocks_schema.txt for the binary block file format.",
        ),
        "schemas/blocks_schema.txt" => (
            "blocks_schema.txt (prose)",
            "Prose specification of the tilezz-rat-dafsa-blocks wire format: manifest fields, block file binary layout, the read algorithm.",
        ),
        "schemas/rat_schema.txt" => (
            "rat_schema.txt (prose)",
            "Length-prefix convention used inside the DAFSA's accepted sequences: stored_sequence = [len, rat...].",
        ),
        "tools/decode.py" => (
            "decode.py (Python 3 decoder)",
            "Standalone, dependency-free Python 3 script that walks the blocked DAFSA in this directory and prints every rat as a line of space-separated signed integers. Run as `python3 tools/decode.py > rats.txt`.",
        ),
        "tools/verify_sha256.py" => (
            "verify_sha256.py (Python 3 hash verifier)",
            "Standalone, dependency-free Python 3 script that checks every sha256 recorded in ro-crate-metadata.json against the on-disk file bytes, flags any file on disk not in the manifest, and reports the producing-tool commit. Run as `python3 tools/verify_sha256.py`; add `--strict` to also fail on non-pristine provenance (`-dirty`/`unknown`). Exits 0 on full match, 1 on mismatch.",
        ),
        "tools/count.py" => (
            "count.py (Python 3 per-perimeter family counter + verifier)",
            "Standalone Python 3 script (stdlib + sibling decode.py) that prints the seven per-perimeter family series (free, oneSided, achiral, rotationSymmetric, symmetric, subring, coset) as OEIS-style terms. `--print` (default, fast) reads `free` off the DAFSA rank index without decoding and cross-checks it against ro-crate-metadata.json, echoing the other six from variableMeasured; `--verify` (slow) decodes every rat, independently re-derives all seven, and checks each against the metadata. Exits 0 on success, 1 on mismatch.",
        ),
        "tools/verify_canonical.py" => (
            "verify_canonical.py (Python 3 canonical-form verifier)",
            "Standalone Python 3 script (stdlib + sibling decode.py) that independently recomputes each stored rat's dihedral-canonical CCW form (lex-min over rotations of the sequence and its reverse, turn-sum > 0) and checks the stored sequence equals it. Catches any non-canonical entry and excludes non-canonical duplicates. Run as `python3 tools/verify_canonical.py`; `--stride N` / `--max N` sample large assets. Exits 0 if all canonical, 1 on violation.",
        ),
        "README.md" => (
            "README.md",
            "Human-readable entry point: dataset overview, author and copyright notice, CC-BY-SA 4.0 license summary, contents map, reproduction recipe (repo URL, exact commit, CLI invocation), and a SHA-256 verification snippet.",
        ),
        "tools/reproduce.sh" => (
            "reproduce.sh (executable rebuild script)",
            "Shell script that re-runs the exact pipeline that produced this dataset: clones the source repo at the recorded commit, builds rat_enum, and runs the original CLI invocation(s). Run as `bash tools/reproduce.sh`. Honours REPO, COMMIT, SRC_DIR env vars for offline / pre-cloned scenarios; honours SOURCE_DATE_EPOCH for bit-identical metadata.",
        ),
        _ if rel_path.starts_with("blocks/") => (
            "DAFSA block file",
            "One block of the rat-DAFSA in tilezz-rat-block binary format, gzip-compressed.",
        ),
        _ => ("", ""),
    }
}

/// Dataset name surfaced as `Dataset.name`. Combines ring, perim
/// bound, step, and canonicalization mode into a stable string.
fn dataset_name(p: &AssetParams) -> String {
    let canon = if p.free { "free" } else { "rotation-canonical" };
    // Presented as the effective ring (ZZ7 for a ZZ14 step-2 run); the
    // ring/step it was computed with live in additionalProperty + the
    // reproduce recipe, not the human-facing name.
    format!(
        "tilezz simple matchstick polygons on Z[zeta_{ring}], perimeter <= {n}, {canon}",
        ring = effective_ring(p),
        n = p.max_steps,
    )
}

/// Stable, machine-friendly identifier for the Dataset.
fn dataset_identifier(p: &AssetParams) -> String {
    let canon = if p.free { "free" } else { "onesided" };
    // Effective-ring slug (zz7, not zz14-step2): the step is an
    // implementation detail. Deployed stepped datasets are the odd rings
    // (effectiveRing odd), so this never collides with a native even ring.
    format!(
        "tilezz-rat-zz{ring}-n{n}-{canon}",
        ring = effective_ring(p),
        n = p.max_steps,
    )
}

/// Euclidean GCD on `u8`, for the `effectiveRing = ring/gcd(step, ring)`
/// derivation. `gcd(x, 0) = x`.
fn gcd_u8(mut a: u8, mut b: u8) -> u8 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// The order of the subring the rats actually live in:
/// `ring / gcd(step, ring)`. `--step 1` gives `ring`; `--step 2` on an
/// even ring ZZ_{2m} gives the odd ring ZZ_m. This is the ring the
/// dataset is PRESENTED as (name, identifier, keywords, description,
/// README); the `ring` / `step` it was computed with are an
/// implementation detail kept in `additionalProperty` and the reproduce
/// recipe.
fn effective_ring(p: &AssetParams) -> u8 {
    p.ring / gcd_u8(p.step.unsigned_abs().max(1), p.ring).max(1)
}

/// One-sentence note on the stored turn-angle units for a stepped
/// (sub-ring) dataset, appended to the RO-Crate description. Empty for
/// `step == 1` datasets, where the stored units already match the
/// presented ring. The dataset presents as Z[zeta_effective] but the
/// turns are physically stored in the parent Z[zeta_ring] units, so a
/// consumer reading the raw sequences needs the conversion factor.
fn turn_units_note_short(p: &AssetParams) -> String {
    let eff = effective_ring(p);
    if eff == p.ring {
        return String::new();
    }
    format!(
        " It is the order-{eff} sub-ring of Z[zeta_{parent}] (Z[zeta_{eff}] has no native lattice here, \
         so it is enumerated as the directions that are multiples of {step}, via --step {step}); each \
         stored turn is a Z[zeta_{parent}] unit (always even) equal to {step} times the Z[zeta_{eff}] turn \
         -- divide stored values by {step} to read Z[zeta_{eff}] turns. See README.md (Turn-angle units).",
        eff = eff,
        parent = p.ring,
        step = p.step,
    )
}

/// Full "Turn-angle units" README section for a stepped (sub-ring)
/// dataset, with a worked example. Empty for `step == 1`.
fn turn_units_note_readme(p: &AssetParams) -> String {
    let eff = effective_ring(p);
    if eff == p.ring {
        return String::new();
    }
    let n = eff as usize;
    let stored = vec![p.step.to_string(); n].join(" ");
    let halved = vec!["1"; n].join(" ");
    format!(
        "\n\n## Turn-angle units\n\n\
         Z[zeta_{eff}] has no native lattice in tilezz, so this dataset is the order-{eff} \
         sub-ring of Z[zeta_{parent}] -- the turn directions that are multiples of {step} -- \
         enumerated with `--step {step}`. Each stored turn is therefore a Z[zeta_{parent}] turn \
         (an integer multiple of `2*pi/{parent}`, always even) and is `{step}` times the \
         corresponding Z[zeta_{eff}] turn (a multiple of `2*pi/{eff}`).\n\n\
         To read the sequences as Z[zeta_{eff}] turns, divide every stored value by `{step}`. \
         For example the regular {eff}-gon is stored as `{stored}` and is `{halved}` in \
         Z[zeta_{eff}]. The web explorer shows the halved Z[zeta_{eff}] form; `tools/decode.py` \
         prints the raw stored Z[zeta_{parent}] values.",
        eff = eff,
        parent = p.ring,
        step = p.step,
        stored = stored,
        halved = halved,
    )
}

/// Structured, machine-readable parameters of an asset, rendered
/// as a schema.org `PropertyValue` array. Used both inside the
/// per-dataset RO-Crate's root Dataset and -- copied verbatim --
/// inside the collection-level RO-Crate's per-dataset stubs.
///
/// Keys: `ring`, `step`, `effectiveRing`, `maxPerimeter`,
/// `canonicalization`, `nSequences`, `maxIndexedLength`.
/// `canonicalization` is a string (`"free"` or `"onesided"`); the
/// rest are integers.
fn dataset_additional_properties(p: &AssetParams, max_indexed_length: usize) -> Value {
    // The order of the actual subring the rats live in. A `--step k` run on
    // ZZn enumerates the order-(n/gcd(k,n)) subring; with step 1 this is just
    // `ring`. The web explorer keys datasets on this so a ZZ14-step2 asset
    // (the odd ring ZZ7) presents as ZZ7 rather than ZZ14.
    let effective_ring = effective_ring(p);
    json!([
        {"@type": "PropertyValue", "name": "ring", "value": p.ring},
        {"@type": "PropertyValue", "name": "step", "value": p.step},
        {"@type": "PropertyValue", "name": "effectiveRing", "value": effective_ring},
        {"@type": "PropertyValue", "name": "maxPerimeter", "value": p.max_steps},
        {
            "@type": "PropertyValue",
            "name": "canonicalization",
            "value": if p.free { "free" } else { "onesided" },
        },
        {"@type": "PropertyValue", "name": "nSequences", "value": p.n_sequences},
        // The deepest perimeter that actually closes a polygon. Usually
        // equals maxPerimeter, but can be smaller when the top perimeter
        // admits no closing polygon (e.g. odd perimeters on an
        // even-perimeter-only ring), so derive it from the data rather
        // than from max_steps. Matches the block manifest.
        {"@type": "PropertyValue", "name": "maxIndexedLength", "value": max_indexed_length},
    ])
}

/// Write `ro-crate-metadata.json` into `dir`, describing every
/// file already present in the directory tree (manifest + schemas
/// + blocks).
///
/// Must be called AFTER `RatDafsa::write_blocks` and AFTER the
/// schemas have been copied into `dir/schemas/` (see
/// [`copy_schemas`]).
/// Per-exact-perimeter counts of the sub-families derivable from a
/// free DAFSA by filtering/bucketing one rat at a time. Emitted into
/// the RO-Crate as `variableMeasured` so each dataset self-documents
/// the integer sequences it realises; re-derived and checked by
/// `tools/count.py --verify`.
#[derive(Default)]
pub struct SequenceCounts {
    free: BTreeMap<usize, u64>,
    subring: BTreeMap<usize, u64>,
    coset: BTreeMap<usize, u64>,
    achiral: BTreeMap<usize, u64>,
    rotation_symmetric: BTreeMap<usize, u64>,
    symmetric: BTreeMap<usize, u64>,
}

fn rot_min(t: &[i8]) -> Vec<i8> {
    let n = t.len();
    (0..n)
        .map(|i| {
            let mut r = Vec::with_capacity(n);
            r.extend_from_slice(&t[i..]);
            r.extend_from_slice(&t[..i]);
            r
        })
        .min()
        .unwrap_or_default()
}

impl SequenceCounts {
    /// Bucket a stream of rats (angle sequences). Streaming, O(1)
    /// memory beyond the count maps.
    pub fn from_rats<I: IntoIterator<Item = Vec<i8>>>(rats: I) -> Self {
        let mut c = Self::default();
        for t in rats {
            let l = t.len();
            *c.free.entry(l).or_default() += 1;
            if t.iter().all(|a| a % 2 == 0) {
                *c.subring.entry(l).or_default() += 1;
            }
            if t.iter().all(|a| a % 2 != 0) {
                *c.coset.entry(l).or_default() += 1;
            }
            // Achiral: the polygon equals its mirror image. The
            // stored rep is free-canonical; it is achiral iff its
            // reverse (the orientation-preserving mirror) is in the
            // same rotation class -- i.e. their rotation-minima match.
            let mut rev = t.clone();
            rev.reverse();
            let achiral = rot_min(&t) == rot_min(&rev);
            // Rotationally symmetric: some nontrivial rotation fixes it.
            let rot_sym = (1..l).any(|d| t[d..].iter().chain(&t[..d]).eq(t.iter()));
            if achiral {
                *c.achiral.entry(l).or_default() += 1;
            }
            if rot_sym {
                *c.rotation_symmetric.entry(l).or_default() += 1;
            }
            if achiral || rot_sym {
                *c.symmetric.entry(l).or_default() += 1;
            }
        }
        c
    }

    /// Render one family as an OEIS-`%S`-style comma-separated string
    /// indexed by perimeter `1..=max_perim` (0 where none close).
    fn series(map: &BTreeMap<usize, u64>, max_perim: usize) -> String {
        (1..=max_perim)
            .map(|n| map.get(&n).copied().unwrap_or(0).to_string())
            .collect::<Vec<_>>()
            .join(",")
    }

    /// The `variableMeasured` JSON array: one schema.org `PropertyValue`
    /// per sub-family, value an ordered perimeter-indexed string.
    ///
    /// Series are sized to the deepest perimeter that actually closes a
    /// polygon (the max `free` key), which equals the block manifest's
    /// `max_indexed_length`. This keeps the emitted strings byte-identical
    /// to what `tools/count.py --verify` re-derives from the same blocks;
    /// sizing to the asset's `max_steps` instead would append spurious
    /// trailing zeros whenever the top perimeter has no closing polygon
    /// (e.g. odd `n` on an even-perimeter-only ring), making the verifier
    /// report a false mismatch.
    fn variable_measured(&self) -> Value {
        let max_perim = self.free.keys().next_back().copied().unwrap_or(0);
        let one_sided: BTreeMap<usize, u64> = (1..=max_perim)
            .filter_map(|n| {
                let f = self.free.get(&n).copied().unwrap_or(0);
                if f == 0 {
                    return None;
                }
                let a = self.achiral.get(&n).copied().unwrap_or(0);
                Some((n, 2 * f - a))
            })
            .collect();
        let pv = |name: &str, map: &BTreeMap<usize, u64>, desc: &str| {
            json!({
                "@type": "PropertyValue",
                "name": name,
                "unitText": "polygons",
                "value": Self::series(map, max_perim),
                "description": desc,
            })
        };
        json!([
            pv(
                "free",
                &self.free,
                "Self-avoiding polygons by exact perimeter (index = perimeter, from 1); the base count, up to rotation and reflection."
            ),
            pv(
                "oneSided",
                &one_sided,
                "Up to rotation only; mirror images counted as distinct (= 2*free - achiral)."
            ),
            pv("achiral", &self.achiral, "Equal to its own mirror image."),
            pv(
                "rotationSymmetric",
                &self.rotation_symmetric,
                "Repetition factor > 1 (nontrivial rotational symmetry)."
            ),
            pv(
                "symmetric",
                &self.symmetric,
                "Has a nontrivial symmetry: achiral OR rotationSymmetric."
            ),
            pv(
                "subring",
                &self.subring,
                "All turns even -- the polygons of the order-n/2 sub-ring."
            ),
            pv(
                "coset",
                &self.coset,
                "All turns odd -- no straight segments; even perimeter only."
            ),
        ])
    }
}

/// Compute the host-coordinate fields (identifier, optional url,
/// distribution, optional sameAs) for a deployed asset.
///
/// Picture: the blocks and provenance of a dataset are fixed history;
/// only the dataset's ADDRESS in the world changes as it is re-hosted.
/// These four fields ARE that address. Factoring them into one helper
/// makes the emit path ([`write_ro_crate`]) and the re-host path
/// ([`rehost_ro_crate`]) write byte-identical host coordinates, so a
/// re-hosted crate is indistinguishable from one emitted fresh at the
/// same host.
///
/// `base_url` is where the dataset DIRECTORY is served (any trailing
/// slash is trimmed); `doi` is a minted DOI like
/// `"10.5281/zenodo.123"`; `slug` is the stable dataset identifier used
/// for the relative fallback; `manifest_size` is `block_index.json`'s
/// byte length for the DataDownload `contentSize`.
///
/// Field rules:
/// - identifier: a DOI PropertyValue when `doi` is set (the most
///   citable id); else the bare `base_url` string when deployed; else a
///   `tilezz-dataset-id` PropertyValue carrying the slug (never a bare
///   relative string).
/// - url: the resolvable `base_url`, only when deployed.
/// - sameAs: `https://doi.org/<doi>`, only when a DOI is set.
/// - distribution: a DataDownload whose contentUrl is the absolute
///   `<base>/block_index.json` when deployed, else the relative
///   `block_index.json`.
fn host_coordinate_fields(
    base_url: Option<&str>,
    doi: Option<&str>,
    slug: &str,
    manifest_size: u64,
) -> (Value, Option<Value>, Value, Option<Value>) {
    let base = base_url.map(|b| b.trim_end_matches('/'));

    let identifier = if let Some(doi) = doi {
        json!({
            "@type": "PropertyValue",
            "propertyID": "DOI",
            "value": doi,
            "url": format!("https://doi.org/{doi}"),
        })
    } else if let Some(base) = base {
        // When deployed at a known base URL, the resolvable URL doubles
        // as the dataset's identifier (best practice).
        json!(base)
    } else {
        // Local / test build with no public home: a structured
        // PropertyValue carrying the stable slug -- never a bare string.
        json!({
            "@type": "PropertyValue",
            "propertyID": "tilezz-dataset-id",
            "value": slug,
        })
    };

    let url = base.map(|b| json!(b));
    let same_as = doi.map(|doi| json!(format!("https://doi.org/{doi}")));

    let content_url = match base {
        Some(base) => format!("{base}/block_index.json"),
        None => "block_index.json".to_string(),
    };
    let distribution = json!({
        "@type": "DataDownload",
        "contentUrl": content_url,
        "encodingFormat": "application/json",
        "contentSize": manifest_size.to_string(),
    });

    (identifier, url, distribution, same_as)
}

pub fn write_ro_crate(
    dir: &Path,
    p: &AssetParams,
    counts: &SequenceCounts,
    base_url: Option<&str>,
    doi: Option<&str>,
) -> io::Result<()> {
    let files = snapshot_dir(dir)?;

    // Dataset hasPart lists every File entity in lex order.
    let has_part: Vec<Value> = files.iter().map(|f| json!({"@id": f.rel_path})).collect();

    let result_ids: Vec<Value> = files
        .iter()
        .filter(|f| !f.rel_path.starts_with("schemas/"))
        .map(|f| json!({"@id": f.rel_path}))
        .collect();

    let mut graph: Vec<Value> = Vec::new();

    // Metadata file descriptor (RO-Crate 1.2 MUST entity).
    graph.push(json!({
        "@type": "CreativeWork",
        "@id": "ro-crate-metadata.json",
        "conformsTo": {"@id": "https://w3id.org/ro/crate/1.2"},
        "about": {"@id": "./"},
    }));

    // Root Data Entity.
    let mut root: BTreeMap<&str, Value> = BTreeMap::new();
    root.insert("@id", json!("./"));
    root.insert("@type", json!("Dataset"));
    root.insert("name", json!(dataset_name(p)));
    // Host coordinates (identifier / url / distribution / sameAs) are
    // shaped by the single source of truth shared with rehost, so an
    // emit and a re-host to the same host produce identical bytes.
    let index_size = std::fs::metadata(dir.join("block_index.json"))
        .map(|m| m.len())
        .unwrap_or(0);
    let (identifier, url, distribution, same_as) =
        host_coordinate_fields(base_url, doi, &dataset_identifier(p), index_size);
    root.insert("identifier", identifier);
    if let Some(url) = url {
        root.insert("url", url);
    }
    root.insert(
        "description",
        json!(format!(
            "Simple matchstick polygons (closed self-avoiding unit-edge polygons) on the cyclotomic ring Z[zeta_{ring}], with perimeter <= {n}, canonicalized by {canon} symmetry. Contains {count} sequences.{units} Self-describing tilezz-rat-dafsa-blocks asset (schemas alongside).",
            ring = effective_ring(p),
            n = p.max_steps,
            canon = if p.free { "free (full dihedral)" } else { "rotation only (one-sided)" },
            count = p.n_sequences,
            units = turn_units_note_short(p),
        )),
    );
    root.insert("datePublished", json!(build_endtime()));
    // Dataset version is the build DATE (YYYY-MM-DD). The exact source
    // commit lives on the #tilezz SoftwareApplication and the #build
    // CreateAction, where it belongs as build provenance; a calendar
    // version is friendlier as the dataset's public `version`.
    root.insert("version", json!(build_date()));
    root.insert(
        "license",
        json!({"@id": "https://creativecommons.org/licenses/by-sa/4.0/"}),
    );
    // Discovery keywords: standard ring notation (Z[zeta_n]) instead
    // of the internal "ZZn" codename, plus the bare OEIS A-number when
    // known so a search for e.g. "A316192" surfaces this dataset.
    let mut keywords = vec![
        json!("combinatorial enumeration"),
        json!("self-avoiding polygon"),
        json!("cyclotomic lattice"),
        json!("matchstick polygon"),
        json!(format!("Z[zeta_{}]", effective_ring(p))),
    ];
    if let Some(oeis) = p.oeis_a_number {
        keywords.push(json!(oeis));
    }
    root.insert("keywords", Value::Array(keywords));
    root.insert("mainEntity", json!({"@id": "block_index.json"}));
    // DataDownload advertising the manifest as the dataset's download
    // endpoint. With a base URL, contentUrl is the absolute deployed
    // location; otherwise it's the relative path inside the crate.
    root.insert("distribution", distribution);
    root.insert("hasPart", json!(has_part));
    // Structured, machine-readable parameters. Crawlers (Google
    // Dataset Search etc.) ignore the inner shape but display the
    // names/values. The collection-level RO-Crate emitted for the
    // hosted page copies this array verbatim into its hasPart stubs
    // so the JS explorer can drive off a single fetch without
    // string-parsing identifiers.
    let max_indexed_length = counts.free.keys().next_back().copied().unwrap_or(0);
    root.insert(
        "additionalProperty",
        dataset_additional_properties(p, max_indexed_length),
    );
    // The integer sequences this dataset realises, one schema.org
    // PropertyValue per sub-family (free / oneSided / achiral /
    // rotationSymmetric / symmetric / subring / coset), values as
    // ordered perimeter-indexed strings. Kept out of
    // `additionalProperty` (which the collection crate copies to the
    // web stub) so it doesn't bloat that; re-derived and verified by
    // `tools/count.py --verify`.
    root.insert("variableMeasured", counts.variable_measured());
    // Author(s) -- the Person(s) credited under the CC-BY-SA
    // attribution clause. Distinct from `instrument` on the
    // CreateAction (= the rat_enum tool).
    root.insert(
        "creator",
        Value::Array(
            authors_as_ids()
                .into_iter()
                .map(|id| json!({"@id": id}))
                .collect(),
        ),
    );
    if let Some(oeis) = p.oeis_a_number {
        let oeis_url = format!("https://oeis.org/{oeis}");
        root.insert("subjectOf", json!({"@id": oeis_url}));
    }
    // sameAs: the DOI's resolvable form, present only on a DOI-minted
    // (e.g. Zenodo) deployment. Kept alongside the identifier so a
    // citation tool can pivot from either.
    if let Some(same_as) = same_as {
        root.insert("sameAs", same_as);
    }
    graph.push(serde_json::Value::Object(
        root.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
    ));

    // License as a CreativeWork (CC-BY-SA 4.0).
    graph.push(json!({
        "@id": "https://creativecommons.org/licenses/by-sa/4.0/",
        "@type": "CreativeWork",
        "name": "Creative Commons Attribution-ShareAlike 4.0 International",
        "identifier": "CC-BY-SA-4.0",
    }));

    // Person contextual entities, one per author parsed from
    // CARGO_PKG_AUTHORS. The Dataset's `creator` array points at
    // these by @id, so the CC-BY-SA attribution requirement has a
    // concrete referent inside the graph.
    for (id, name, email) in parse_authors(PKG.authors) {
        let mut p_obj: BTreeMap<&str, Value> = BTreeMap::new();
        // When the @id is an ORCID URI, restate it as `identifier`
        // so the globally-unique researcher iD is explicit, not just
        // implied by the @id.
        if id.starts_with("https://orcid.org/") {
            p_obj.insert("identifier", json!(id));
        }
        p_obj.insert("@id", json!(id));
        p_obj.insert("@type", json!("Person"));
        p_obj.insert("name", json!(name));
        if let Some(e) = email {
            p_obj.insert("email", json!(e));
        }
        graph.push(serde_json::Value::Object(
            p_obj.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
        ));
    }

    // SoftwareApplication for the rat_enum binary. Repository URL +
    // version come straight from Cargo.toml via env!(); the commit
    // hash from build.rs. softwareRequirements names what an
    // archivist needs to install before running rat_enum.
    graph.push(json!({
        "@id": format!("#{}", PKG.name),
        "@type": "SoftwareApplication",
        "name": PKG.name,
        "version": PKG.version,
        "description": PKG.description,
        "url": PKG.repository,
        "codeRepository": PKG.repository,
        "softwareVersion": PKG.commit,
        "softwareRequirements": "Rust toolchain (stable, 2024-12+); standard Cargo build prerequisites (git, C linker, pkg-config).",
        "downloadUrl": format!("{}/archive/{}.tar.gz", PKG.repository, PKG.commit),
    }));

    // CreateAction recording the exact run that produced this asset.
    // `description` carries the verbatim CLI invocation, so an
    // archivist can rerun without consulting the source. README.md
    // in the directory has the human-readable version with prereqs.
    graph.push(json!({
        "@id": "#build",
        "@type": "CreateAction",
        "name": format!("Generate {}", dataset_identifier(p)),
        "description": format!(
            "Produced by: {}\nReproduce: bash tools/reproduce.sh (in this directory). See README.md for context.",
            reproduce_one_liner(p),
        ),
        "endTime": build_endtime(),
        "instrument": {"@id": format!("#{}", PKG.name)},
        "result": result_ids,
    }));

    // Optional OEIS cross-reference.
    if let Some(oeis) = p.oeis_a_number {
        let oeis_url = format!("https://oeis.org/{oeis}");
        graph.push(json!({
            "@id": oeis_url,
            "@type": "CreativeWork",
            "name": format!("OEIS {oeis}"),
            "identifier": oeis,
            "url": oeis_url.clone(),
        }));
    }

    // File entities.
    for f in &files {
        let (name, description) = human_label_for(&f.rel_path);
        let mut obj: BTreeMap<&str, Value> = BTreeMap::new();
        let is_schema = f.rel_path.starts_with("schemas/");
        // Schemas double as CreativeWork so other files can
        // `conformsTo` them by `@id`.
        obj.insert(
            "@type",
            if is_schema {
                json!(["File", "CreativeWork"])
            } else {
                json!("File")
            },
        );
        obj.insert("@id", json!(f.rel_path));
        if !name.is_empty() {
            obj.insert("name", json!(name));
            obj.insert("description", json!(description));
        }
        obj.insert("encodingFormat", json!(encoding_format_of(&f.rel_path)));
        obj.insert("contentSize", json!(f.size.to_string()));
        obj.insert("sha256", json!(f.sha256));
        let conforms = conforms_to_of(&f.rel_path);
        if !conforms.is_empty() {
            obj.insert(
                "conformsTo",
                if conforms.len() == 1 {
                    conforms.into_iter().next().unwrap()
                } else {
                    Value::Array(conforms)
                },
            );
        }
        graph.push(serde_json::Value::Object(
            obj.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
        ));
    }

    let root_obj = json!({
        "@context": "https://w3id.org/ro/crate/1.2/context",
        "@graph": graph,
    });

    let path = dir.join("ro-crate-metadata.json");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&path)?);
    serde_json::to_writer_pretty(&mut writer, &root_obj)
        .map_err(|e| io::Error::other(format!("write ro-crate-metadata.json: {e}")))?;
    io::Write::write_all(&mut writer, b"\n")?;
    Ok(())
}

/// Recover the dataset slug (the relative-fallback identifier) from an
/// already-written root Dataset entity. Prefers the value carried by an
/// existing `tilezz-dataset-id` PropertyValue identifier (the form the
/// emitter writes for a local build); falls back to the dataset `name`.
/// Used only to repopulate the relative identifier when a re-host clears
/// a DOI / base_url back to the location-independent form.
/// Recover the dataset slug (e.g. "tilezz-rat-zz7-n10-free") from the
/// `#build` CreateAction name ("Generate <slug>"). That field is part of
/// the immutable build provenance and is never touched by re-hosting --
/// unlike `identifier`, which a prior DOI re-host overwrites with the
/// DOI -- so it is the only reliable slug source when re-hosting an
/// already-re-hosted crate.
fn recover_slug(graph: &[Value]) -> String {
    for e in graph {
        if e.get("@id").and_then(|v| v.as_str()) == Some("#build")
            && let Some(slug) = e
                .get("name")
                .and_then(|v| v.as_str())
                .and_then(|n| n.strip_prefix("Generate "))
        {
            return slug.to_string();
        }
    }
    "tilezz-rat-dataset".to_string()
}

/// Re-host a dataset to a new host (or back to location-independent
/// form) by SURGICALLY patching only its host-coordinate fields.
///
/// Picture: the blocks are content-addressed and the provenance (the
/// producing commit, build time, reproduce recipe, file sha256s) is
/// fixed history. Re-hosting -- e.g. mirroring a raw.githubusercontent
/// dataset to Zenodo with a DOI -- changes only WHERE the dataset
/// lives, never WHAT it is. So this rewrites exactly `identifier`,
/// `url`, `distribution`, and `sameAs` on the root Dataset and leaves
/// every other field and entity byte-for-byte unchanged. It does NOT
/// re-run the emitter (whose build time / commit would differ and
/// corrupt the recorded provenance).
///
/// `base_url` is the new directory URL (trailing slash trimmed);
/// `doi` is a minted DOI like `"10.5281/zenodo.123"`. Passing both
/// `None` reverts the crate to relative / PropertyValue forms and
/// drops `url` / `sameAs`, so a stale absolute deployment is cleared
/// cleanly. The JSON is written back with the same pretty formatting +
/// trailing newline the emitter uses, so a diff shows only the host
/// fields.
pub fn rehost_ro_crate(dir: &Path, base_url: Option<&str>, doi: Option<&str>) -> io::Result<()> {
    let path = dir.join("ro-crate-metadata.json");
    if !path.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "{} is not an RO-Crate: no ro-crate-metadata.json found",
                dir.display()
            ),
        ));
    }
    let bytes = std::fs::read(&path)?;
    let mut crate_json: Value = serde_json::from_slice(&bytes)
        .map_err(|e| io::Error::other(format!("{} is not valid JSON: {e}", path.display())))?;

    // block_index.json byte length feeds the DataDownload contentSize;
    // it is part of the dataset payload, so it must already exist next
    // to the crate.
    let manifest_size = std::fs::metadata(dir.join("block_index.json"))
        .map_err(|e| {
            io::Error::other(format!(
                "cannot stat {}/block_index.json (re-host needs the dataset payload present): {e}",
                dir.display()
            ))
        })?
        .len();

    // Locate the root Dataset entity (@id == "./") and patch only its
    // host-coordinate keys; the @graph order and all other entities are
    // preserved untouched.
    // Recover the slug from the immutable #build provenance BEFORE taking
    // the mutable root borrow (a prior DOI re-host may have overwritten
    // the root identifier, so it is not a reliable slug source).
    let slug = crate_json
        .get("@graph")
        .and_then(|g| g.as_array())
        .map(|g| recover_slug(g))
        .unwrap_or_else(|| "tilezz-rat-dataset".to_string());

    let graph = crate_json
        .get_mut("@graph")
        .and_then(|g| g.as_array_mut())
        .ok_or_else(|| io::Error::other(format!("{} has no @graph array", path.display())))?;
    let root = graph
        .iter_mut()
        .find_map(|e| {
            let obj = e.as_object_mut()?;
            (obj.get("@id").and_then(|v| v.as_str()) == Some("./")).then_some(obj)
        })
        .ok_or_else(|| {
            io::Error::other(format!(
                "{} has no Dataset entity at @id \"./\"",
                path.display()
            ))
        })?;

    let (identifier, url, distribution, same_as) =
        host_coordinate_fields(base_url, doi, &slug, manifest_size);

    // Overwrite identifier + distribution; set-or-remove url + sameAs so
    // a re-host to the relative form clears any stale absolute values.
    root.insert("identifier".to_string(), identifier);
    root.insert("distribution".to_string(), distribution);
    match url {
        Some(v) => {
            root.insert("url".to_string(), v);
        }
        None => {
            root.remove("url");
        }
    }
    match same_as {
        Some(v) => {
            root.insert("sameAs".to_string(), v);
        }
        None => {
            root.remove("sameAs");
        }
    }

    // Write back with the exact formatting the emitter uses (pretty +
    // trailing newline) so a diff shows only the host fields.
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&path)?);
    serde_json::to_writer_pretty(&mut writer, &crate_json)
        .map_err(|e| io::Error::other(format!("write {}: {e}", path.display())))?;
    io::Write::write_all(&mut writer, b"\n")?;
    Ok(())
}

/// Emit all the archival-extras files that make the asset
/// self-describing, beyond the bare `block_index.json` +
/// `blocks/<sha256>.bin` wire format:
///
/// - `schemas/blocks_schema.txt`, `schemas/rat_schema.txt`,
///   `schemas/block_index.schema.json` -- referenced by the
///   manifest's `conformsTo` chain.
/// - `tools/decode.py` -- standalone Python 3 decoder (no deps)
///   that converts a blocked-DAFSA directory into one rat per
///   line of plain text. Lets an archivist extract sequences
///   without compiling Rust.
/// - `README.md` -- human-readable entry: dataset summary, author
///   and copyright notice, CC-BY-SA 4.0 license terms, an exact
///   rebuild recipe (prerequisites, repo + commit, CLI invocation
///   templated from `params`), and a SHA-256 verification snippet.
///
/// Idempotent and deterministic. Run before [`write_ro_crate`] so
/// it can include these files in `hasPart` / `conformsTo`.
pub fn write_archival_extras(dir: &Path, params: &AssetParams) -> io::Result<()> {
    let schemas = dir.join("schemas");
    std::fs::create_dir_all(&schemas)?;
    std::fs::write(schemas.join("blocks_schema.txt"), BLOCKS_SCHEMA_TXT)?;
    std::fs::write(schemas.join("rat_schema.txt"), RAT_SCHEMA_TXT)?;
    std::fs::write(
        schemas.join("block_index.schema.json"),
        BLOCK_INDEX_SCHEMA_JSON,
    )?;
    let tools = dir.join("tools");
    std::fs::create_dir_all(&tools)?;
    std::fs::write(tools.join("decode.py"), DECODER_PY)?;
    std::fs::write(tools.join("verify_sha256.py"), VERIFY_SHA256_PY)?;
    std::fs::write(tools.join("count.py"), COUNT_PY)?;
    std::fs::write(tools.join("verify_canonical.py"), VERIFY_CANONICAL_PY)?;
    std::fs::write(dir.join("README.md"), readme_md(params))?;
    let sh_path = tools.join("reproduce.sh");
    std::fs::write(&sh_path, reproduce_sh(params))?;
    // Best-effort chmod +x. On non-Unix targets std::fs has no
    // direct equivalent; skip rather than fail (the script is
    // still runnable via `sh tools/reproduce.sh`).
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perm = std::fs::metadata(&sh_path)?.permissions();
        perm.set_mode(0o755);
        std::fs::set_permissions(&sh_path, perm)?;
    }
    Ok(())
}

/// Build the top-level RO-Crate that describes a hosted page (e.g.
/// `web/`) and links the per-dataset crates it serves. The result
/// is a single `web_dir/ro-crate-metadata.json` whose `@graph`
/// contains:
///
/// * The page's root Dataset entity (`@id: "./"`) with `mainEntity`
///   = the WebApplication and `hasPart` = one stub per dataset
///   discovered under `web_dir/data/*/`.
/// * The WebApplication entity (the explorer itself).
/// * Person + SoftwareSourceCode contextual entities, sourced from
///   Cargo.toml (same identities as the per-dataset crates use).
/// * One Dataset stub per discovered dataset, carrying the
///   `additionalProperty` array verbatim from the child crate so
///   the JS UI can drive off a single fetch.
///
/// `web_dir` must already exist and contain `data/<name>/ro-crate-
/// metadata.json` files (produced by [`write_ro_crate`]). Missing
/// or malformed child crates are skipped with a warning to stderr;
/// the function does not fail on them so a partial deploy still
/// emits a usable top-level crate.
///
/// `page_url` is the canonical URL where the page is hosted. Used
/// only to populate the WebApplication entity's `url` property; it
/// does not affect file layout.
pub fn write_collection_ro_crate(web_dir: &Path, page_url: &str) -> io::Result<()> {
    let mut has_part_ids: Vec<Value> = Vec::new();
    let mut child_stubs: Vec<Value> = Vec::new();

    // Sorted directory scan so the resulting crate is deterministic
    // regardless of filesystem enumeration order.
    let data_dir = web_dir.join("data");
    let mut child_dirs: Vec<std::path::PathBuf> = Vec::new();
    if data_dir.is_dir() {
        for entry in std::fs::read_dir(&data_dir)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() && path.join("ro-crate-metadata.json").is_file() {
                child_dirs.push(path);
            }
        }
    }
    child_dirs.sort();

    for path in &child_dirs {
        let dir_name = match path.file_name().and_then(|s| s.to_str()) {
            Some(n) => n.to_string(),
            None => continue,
        };
        let stub_id = format!("./data/{dir_name}/");
        let manifest_id = format!("./data/{dir_name}/ro-crate-metadata.json");
        let crate_path = path.join("ro-crate-metadata.json");
        let bytes = match std::fs::read(&crate_path) {
            Ok(b) => b,
            Err(e) => {
                eprintln!("warn: cannot read {}: {e}", crate_path.display());
                continue;
            }
        };
        let crate_json: Value = match serde_json::from_slice(&bytes) {
            Ok(v) => v,
            Err(e) => {
                eprintln!("warn: {} is not valid JSON: {e}", crate_path.display());
                continue;
            }
        };
        let root = match find_root_dataset(&crate_json) {
            Some(r) => r,
            None => {
                eprintln!(
                    "warn: {} has no Dataset entity at @id \"./\"; skipping",
                    crate_path.display()
                );
                continue;
            }
        };

        // Build a stub. Keep schema.org-typed properties only;
        // ignore RO-Crate-internal fields like hasPart (the stub
        // points at the child crate's manifest, not at its files).
        let mut stub = serde_json::Map::new();
        stub.insert("@id".into(), json!(stub_id));
        stub.insert("@type".into(), json!("Dataset"));
        for key in [
            "identifier",
            "name",
            "description",
            "license",
            "datePublished",
            "version",
            "keywords",
            "creator",
            "subjectOf",
            "additionalProperty",
        ] {
            if let Some(v) = root.get(key) {
                stub.insert(key.to_string(), v.clone());
            }
        }
        stub.insert("encodingFormat".into(), json!("tilezz-rat-dafsa-blocks"));
        // Distribution: a DataDownload pointing at the child crate's
        // canonical, self-describing manifest. Crawlers follow it
        // for the full file inventory + SHA-256s.
        stub.insert(
            "distribution".into(),
            json!({
                "@type": "DataDownload",
                "encodingFormat": "application/ld+json",
                "contentUrl": manifest_id.clone(),
                "description": "RO-Crate 1.2 manifest for this dataset (canonical, self-describing).",
            }),
        );
        has_part_ids.push(json!({"@id": stub_id}));
        child_stubs.push(Value::Object(stub));
    }

    let mut graph: Vec<Value> = Vec::new();

    // RO-Crate metadata file descriptor (required entity).
    graph.push(json!({
        "@type": "CreativeWork",
        "@id": "ro-crate-metadata.json",
        "conformsTo": {"@id": "https://w3id.org/ro/crate/1.2"},
        "about": {"@id": "./"},
    }));

    // Page root Dataset.
    let mut root: BTreeMap<&str, Value> = BTreeMap::new();
    root.insert("@id", json!("./"));
    root.insert("@type", json!("Dataset"));
    root.insert("name", json!("tilezz Rat Explorer"));
    root.insert(
        "description",
        json!(
            "Interactive WebAssembly explorer for simple matchstick polygons on cyclotomic \
             lattices Z[zeta_n]. Hosts a queryable database of canonical polygons; each \
             dataset under data/ is a self-describing RO-Crate sub-dataset (follow its \
             ro-crate-metadata.json for the full file inventory)."
        ),
    );
    root.insert("identifier", json!(page_url));
    root.insert("url", json!(page_url));
    root.insert("datePublished", json!(build_endtime()));
    // Date-based version (see write_ro_crate); the source commit
    // stays on the #tilezz SoftwareSourceCode entity below.
    root.insert("version", json!(build_date()));
    // An RO-Crate root Dataset MUST carry a license. The hosted page
    // and its database share the per-dataset CC-BY-SA 4.0 terms.
    root.insert(
        "license",
        json!({"@id": "https://creativecommons.org/licenses/by-sa/4.0/"}),
    );
    root.insert("isBasedOn", json!({"@id": format!("#{}", PKG.name)}));
    root.insert(
        "creator",
        Value::Array(
            authors_as_ids()
                .into_iter()
                .map(|id| json!({"@id": id}))
                .collect(),
        ),
    );
    root.insert("mainEntity", json!({"@id": "#rat-explorer"}));
    root.insert("hasPart", json!(has_part_ids));
    graph.push(serde_json::Value::Object(
        root.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
    ));

    // WebApplication entity for the page itself.
    graph.push(json!({
        "@id": "#rat-explorer",
        "@type": "WebApplication",
        "name": "tilezz Rat Explorer",
        "description": "Interactive explorer for simple matchstick polygons on cyclotomic lattices. Builds polygons from angle sequences and looks them up in a packaged RO-Crate database of canonical forms.",
        "url": page_url,
        "applicationCategory": "BrowserApplication",
        "operatingSystem": "Any (browser with WebAssembly)",
        "isBasedOn": {"@id": format!("#{}", PKG.name)},
        "author": authors_as_ids().into_iter().map(|id| json!({"@id": id})).collect::<Vec<_>>(),
        "license": {"@id": "https://opensource.org/license/mit"},
    }));

    // Person contextual entities, one per author. Identical to the
    // per-dataset crate's #author-N identities, so a crawler that
    // resolves both crates de-duplicates them via @id.
    for (id, name, email) in parse_authors(PKG.authors) {
        let mut p_obj: BTreeMap<&str, Value> = BTreeMap::new();
        if id.starts_with("https://orcid.org/") {
            p_obj.insert("identifier", json!(id));
        }
        p_obj.insert("@id", json!(id));
        p_obj.insert("@type", json!("Person"));
        p_obj.insert("name", json!(name));
        if let Some(e) = email {
            p_obj.insert("email", json!(e));
        }
        graph.push(serde_json::Value::Object(
            p_obj.into_iter().map(|(k, v)| (k.to_string(), v)).collect(),
        ));
    }

    // SoftwareSourceCode entity. Distinct from the per-dataset
    // crate's SoftwareApplication for rat_enum: here we describe
    // the repository itself (a code artifact a developer can clone
    // and build), not a particular CLI invocation.
    graph.push(json!({
        "@id": format!("#{}", PKG.name),
        "@type": "SoftwareSourceCode",
        "name": PKG.name,
        "description": PKG.description,
        "codeRepository": PKG.repository,
        "programmingLanguage": "Rust",
        "version": PKG.version,
        "softwareVersion": PKG.commit,
        "license": {"@id": "https://opensource.org/license/mit"},
    }));

    // License contextual entity (CC-BY-SA 4.0), mirroring the
    // per-dataset crate. Guard against a duplicate @id in case a
    // child stub or other entity already introduced it.
    let cc_id = "https://creativecommons.org/licenses/by-sa/4.0/";
    let already_present = graph
        .iter()
        .any(|e| e.get("@id").and_then(|v| v.as_str()) == Some(cc_id));
    if !already_present {
        graph.push(json!({
            "@id": cc_id,
            "@type": "CreativeWork",
            "name": "Creative Commons Attribution-ShareAlike 4.0 International",
            "identifier": "CC-BY-SA-4.0",
        }));
    }

    // Per-dataset stubs.
    graph.extend(child_stubs);

    let root_obj = json!({
        "@context": "https://w3id.org/ro/crate/1.2/context",
        "@graph": graph,
    });
    let out_path = web_dir.join("ro-crate-metadata.json");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&out_path)?);
    serde_json::to_writer_pretty(&mut writer, &root_obj)
        .map_err(|e| io::Error::other(format!("write {}: {e}", out_path.display())))?;
    io::Write::write_all(&mut writer, b"\n")?;
    Ok(())
}

/// Find the `@graph` entity that describes the root directory
/// (`@id: "./"`) of a per-dataset RO-Crate. Returns `None` if the
/// document is not shaped like an RO-Crate.
fn find_root_dataset(crate_json: &Value) -> Option<&serde_json::Map<String, Value>> {
    let graph = crate_json.get("@graph")?.as_array()?;
    for entity in graph {
        let obj = entity.as_object()?;
        if obj.get("@id").and_then(|v| v.as_str()) == Some("./") {
            return Some(obj);
        }
    }
    None
}

/// Shared flag tail across all `rat_enum` invocations for a given
/// asset: ring, perim, canonicalization, optional step / OEIS ref.
/// (Block-size and threads are mode-specific and added by the
/// caller.)
fn shared_flags(p: &AssetParams) -> String {
    let canon = if p.free { " --free" } else { "" };
    let step = if p.step == 1 {
        String::new()
    } else {
        format!(" --step {}", p.step)
    };
    let oeis = match p.oeis_a_number {
        Some(a) => format!(" --oeis-a-number {a}"),
        None => String::new(),
    };
    format!(
        "--ring {ring} -n {n}{canon}{step}{oeis}",
        ring = p.ring,
        n = p.max_steps
    )
}

/// Verbatim CLI invocation(s) that produced an asset with the
/// given parameters. One string per shell command. The streaming
/// pipeline returns three commands (stream / merge / build);
/// the in-memory path returns one (`--mode dafsa-blocks`).
fn reproduce_commands(p: &AssetParams) -> Vec<String> {
    let shared = shared_flags(p);
    let ident = dataset_identifier(p);
    let tbb = p.target_block_bytes;
    match p.produced_via {
        ProducedVia::InMemory => vec![format!(
            "./target/release/rat_enum {shared} --mode dafsa-blocks \
                 --target-block-bytes {tbb} --threads 0 -o {ident}"
        )],
        ProducedVia::StreamingPipeline => vec![
            format!(
                "./target/release/rat_enum {shared} --mode stream \
                 --threads 16 -o {ident}-pipeline"
            ),
            format!(
                "./target/release/rat_enum {shared} --mode merge \
                 -o {ident}-pipeline"
            ),
            format!(
                "./target/release/rat_enum {shared} --mode build \
                 --target-block-bytes {tbb} -o {ident}-pipeline"
            ),
            // After --mode build, the asset is at {ident}-pipeline/dafsa/.
            // The final mv lines mirror that into the canonical path so
            // the reproduced output has the same top-level directory
            // name as the original.
            format!("mv {ident}-pipeline/dafsa {ident}"),
        ],
    }
}

/// Single-line summary for the CreateAction.description -- joins
/// multi-step pipelines with ` && ` so consumers that don't parse
/// the array still get a runnable shell oneliner.
fn reproduce_one_liner(p: &AssetParams) -> String {
    reproduce_commands(p).join(" && ")
}

/// Executable `reproduce.sh` body. Clones the repo at the recorded
/// commit, builds rat_enum, and runs the reproduction step(s).
/// Idempotent: safe to rerun.
fn reproduce_sh(p: &AssetParams) -> String {
    let commit = PKG.commit;
    let repo = PKG.repository;
    let name = PKG.name;
    let ident = dataset_identifier(p);
    let pipeline_label = match p.produced_via {
        ProducedVia::InMemory => "in-memory (single-step --mode dafsa-blocks)",
        ProducedVia::StreamingPipeline => "streaming pipeline (--mode stream + merge + build)",
    };
    // Each command on its own line with a `\` line-continuation
    // would be ugly here -- one line per command is simpler and
    // copy/paste-friendly out of `bash -x reproduce.sh` traces.
    let mut steps = String::new();
    for cmd in reproduce_commands(p) {
        steps.push_str(&format!("{cmd}\n"));
    }

    format!(
        "#!/bin/sh
#
# Rebuild this dataset ({ident}) from source.
# Reproduction strategy: {pipeline_label}.
#
# Tested with the same tilezz commit that produced the original.
# Re-running this script in a clean directory should yield a
# directory whose contents match the recorded sha256s in
# ro-crate-metadata.json (see the verification snippet in
# README.md).
#
# For bit-identical ro-crate-metadata.json, set
# SOURCE_DATE_EPOCH to a fixed value before running -- otherwise
# CreateAction.endTime will drift, but the block files and
# block_index.json are unaffected.

set -eu

REPO=\"${{REPO:-{repo}}}\"
COMMIT=\"${{COMMIT:-{commit}}}\"
SRC_DIR=\"${{SRC_DIR:-{name}-{commit}}}\"

# The recorded commit may carry a `-dirty` suffix (original built from
# an unclean tree); strip it for the checkout, but keep it for the
# provenance check so a dirty original is reported loudly.
BASE_COMMIT=\"${{COMMIT%-dirty}}\"
if [ \"$BASE_COMMIT\" != \"$COMMIT\" ]; then
    echo \"WARNING: recorded commit $COMMIT is marked -dirty -- the original was\" >&2
    echo \"         built from an unclean tree and may not reproduce exactly.\" >&2
fi

if [ ! -d \"$SRC_DIR/.git\" ]; then
    git clone \"$REPO\" \"$SRC_DIR\"
fi
( cd \"$SRC_DIR\" && git fetch && git checkout \"$BASE_COMMIT\" && \\
  cargo build --release --bin rat_enum --features cli )

cd \"$SRC_DIR\"

# Provenance guard (paranoid): refuse to reproduce unless the binary we
# just built self-reports the expected commit. Catches a stale/wrong
# rat_enum picked up from PATH, a build-cache mixup, or a dirty tree --
# the failure modes that silently corrupt provenance.
SELF_VERSION=\"$(./target/release/rat_enum --version)\"
case \"$SELF_VERSION\" in
    *\"$BASE_COMMIT\"*) : ;;
    *) echo \"ERROR: built rat_enum reports [$SELF_VERSION] but this dataset was\" >&2
       echo \"       produced at commit $BASE_COMMIT. Refusing to reproduce with a\" >&2
       echo \"       mismatched binary.\" >&2
       exit 1 ;;
esac
case \"$SELF_VERSION\" in
    *-dirty*) echo \"WARNING: rebuilt rat_enum is -dirty; output may not match exactly.\" >&2 ;;
esac

# Run the reproduction step(s) relative to the source tree so the
# `./target/release/rat_enum` path resolves.
{steps}\
echo
echo \"Reproduced asset: $(pwd)/{ident}\"
echo \"Compare against the recorded sha256s via the snippet in README.md.\"
",
    )
}

/// Render the human-readable `README.md`. Title, author/copyright
/// (so the CC-BY-SA attribution clause has someone to attribute),
/// license summary, contents map, reproduction recipe, and a
/// verification snippet. Templated by the asset's parameters so a
/// directory's README always describes exactly that directory.
fn readme_md(p: &AssetParams) -> String {
    let commit = PKG.commit;
    let repo = PKG.repository;
    let name = PKG.name;
    let version = PKG.version;
    let pipeline_blurb = match p.produced_via {
        ProducedVia::InMemory => "single-step in-memory build (`--mode dafsa-blocks`)",
        ProducedVia::StreamingPipeline => {
            "three-stage streaming pipeline (`--mode stream` -> `--mode merge` -> `--mode build`)"
        }
    };
    let authors = format_authors(PKG.authors);
    let year = current_year();

    let name_h = dataset_name(p);
    let ident = dataset_identifier(p);
    let canon = if p.free {
        "free (full dihedral symmetry reduction)"
    } else {
        "rotation-canonical (one-sided)"
    };
    let oeis_note = match p.oeis_a_number {
        Some(a) => format!(
            "\n\nUpstream cross-reference: this dataset realises OEIS sequence \
             [`{a}`](https://oeis.org/{a}).",
        ),
        None => String::new(),
    };

    // Turn-angle units section for stepped (sub-ring) datasets; empty otherwise.
    let units = turn_units_note_readme(p);

    format!(
        "# {name_h}

Simple matchstick polygons -- closed self-avoiding polygonal walks \
with unit-length edges and turn angles in integer multiples of \
`2*pi/{ring}` -- on the cyclotomic ring `Z[zeta_{ring}]`, with perimeter \
up to {n}, canonicalised under {canon} symmetry. {count} sequences.\
{oeis_note}{units}

## Copyright

Copyright (c) {year} {authors}. All rights reserved subject to the \
license below.

## License

This dataset is distributed under the [Creative Commons \
Attribution-ShareAlike 4.0 International \
License](https://creativecommons.org/licenses/by-sa/4.0/) \
(CC-BY-SA-4.0). In short:

- **Attribution** -- credit the original author(s) above and \
  link back to the source repository ({repo}).
- **ShareAlike** -- if you remix, transform, or build upon this \
  dataset, distribute your contributions under the same license.

The full license text is at <https://creativecommons.org/licenses/by-sa/4.0/legalcode>.

## Contents

This directory is a [RO-Crate 1.2](https://www.researchobject.org/ro-crate/specification/1.2/) \
asset. The entry point is `ro-crate-metadata.json`; every file \
listed below is also recorded there with a `sha256`, an \
`encodingFormat`, and (where applicable) a `conformsTo` pointer \
to its schema.

```
{ident}/
  README.md               this file
  ro-crate-metadata.json  RO-Crate 1.2 manifest (start here for tooling)
  block_index.json        DAFSA wire manifest (counts, root state, sha256 block index)
  schemas/
    block_index.schema.json  formal JSON Schema (draft 2020-12)
    blocks_schema.txt        prose spec covering JSON + .bin formats
    rat_schema.txt           length-prefix convention inside the DAFSA
  blocks/
    <sha256>.bin            one gzipped DAFSA block each; filename = SHA-256 of file
  tools/
    decode.py               standalone Python 3 decoder (no deps)
    verify_sha256.py        SHA-256 verifier (no deps; exits 0 on full match)
    count.py                per-perimeter family terms (--print) + re-derive/verify (--verify)
    verify_canonical.py     independently check every rat is dihedral-canonical CCW
    reproduce.sh            executable rebuild script (clones + builds + runs)
```

## Extracting sequences (no Rust toolchain needed)

`tools/decode.py` walks the blocked DAFSA and prints every \
sequence as a line of space-separated signed integers:

```sh
python3 tools/decode.py > rats.txt
```

The line count of `rats.txt` must equal `n_sequences` from \
`block_index.json`.

## Verifying SHA-256s

Every File entity in `ro-crate-metadata.json` carries a `sha256` \
that matches the on-disk bytes. `tools/verify_sha256.py` checks \
the whole set:

```sh
python3 tools/verify_sha256.py
```

Exits 0 on full match, 1 on any mismatch.

## Reproducing from source

Produced by [`{name}`]({repo}) v{version}, commit `{commit}`, via \
{pipeline_blurb}.

The `reproduce.sh` script in this directory is a self-contained \
recipe: it clones the repo at the recorded commit, builds \
`rat_enum`, and runs the exact sequence of commands that \
produced the dataset.

Prerequisites: a recent Rust toolchain (stable, 2024-12 or \
later); `git`, a C linker, and `pkg-config` (the standard Cargo \
build deps).

```sh
bash tools/reproduce.sh
```

(The script honours `REPO`, `COMMIT`, and `SRC_DIR` environment \
variables for mirroring / vendoring / pre-cloned checkouts; see \
the script header.)

The reproduced directory must match the existing one \
file-for-file. If not, the most likely cause is a different Rust \
version producing a different state ordering inside the DAFSA -- \
pin to the toolchain version this commit's `Cargo.lock` \
specifies. For bit-identical `ro-crate-metadata.json` across \
reruns, set `SOURCE_DATE_EPOCH` to a fixed POSIX-seconds value \
before running:

```sh
SOURCE_DATE_EPOCH=1780000000 bash tools/reproduce.sh
```

Otherwise the `CreateAction.endTime` will reflect the current \
wall-clock time and the metadata hash will drift; the block \
files and `block_index.json` are unaffected.
",
        ring = effective_ring(p),
        n = p.max_steps,
        count = p.n_sequences,
    )
}

/// Parse Cargo's colon-separated `authors` string into a list of
/// `(@id, name, optional_email)` tuples for Person entities in the
/// RO-Crate `@graph`. The `@id` is a deterministic anchor like
/// `#author-0`, `#author-1`, ... so multiple authors get distinct
/// referents.
fn parse_authors(authors: &str) -> Vec<(String, String, Option<String>)> {
    authors
        .split(':')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .enumerate()
        .map(|(i, raw)| {
            // Cargo's convention is `Name <email>`; either piece
            // may be missing.
            let (name, email) = match (raw.find('<'), raw.rfind('>')) {
                (Some(lt), Some(gt)) if lt < gt => (
                    raw[..lt].trim().to_string(),
                    Some(raw[lt + 1..gt].trim().to_string()),
                ),
                _ => (raw.to_string(), None),
            };
            // Prefer the resolvable ORCID URI as the Person's @id when
            // we know one for this email; fall back to the deterministic
            // `#author-N` anchor otherwise.
            let id = match orcid_for(email.as_deref()) {
                Some(orcid) => orcid.to_string(),
                None => format!("#author-{i}"),
            };
            (id, name, email)
        })
        .collect()
}

/// Just the @id of each author entity, in the same order as
/// `parse_authors`. Used to populate `Dataset.creator`.
fn authors_as_ids() -> Vec<String> {
    parse_authors(PKG.authors)
        .into_iter()
        .map(|(id, _, _)| id)
        .collect()
}

/// Parse Cargo's colon-separated `authors` string into a
/// comma-separated, email-stripped human display form. e.g.
/// `"Anton Pirogov <apirogov@users.noreply.github.com>:Other <x@y>"`
/// becomes `"Anton Pirogov, Other"`. Used in the README header
/// and the RO-Crate Person entities.
fn format_authors(authors: &str) -> String {
    authors
        .split(':')
        .map(|a| a.trim())
        .filter(|a| !a.is_empty())
        .map(|a| {
            // Strip any "<email>" suffix.
            if let Some(idx) = a.find('<') {
                a[..idx].trim().to_string()
            } else {
                a.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join(", ")
}

/// Year of the build timestamp, for the copyright notice. Uses
/// the same SOURCE_DATE_EPOCH-aware path as `build_endtime` so a
/// reproducible build sets a deterministic copyright year.
fn current_year() -> i64 {
    let secs: i64 = match std::env::var("SOURCE_DATE_EPOCH")
        .ok()
        .and_then(|s| s.parse().ok())
    {
        Some(s) => s,
        None => std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs() as i64)
            .unwrap_or(0),
    };
    let days = secs.div_euclid(86_400);
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = (yoe as i64) + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    if month <= 2 { y + 1 } else { y }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iso_utc_examples() {
        // 1970-01-01T00:00:00Z
        assert_eq!(format_iso_utc(0), "1970-01-01T00:00:00Z");
        // 2026-06-02T19:30:00Z -- a fixed reference point.
        // 56 years + leap days. Compute against
        // SOURCE_DATE_EPOCH-compatible values.
        let secs_2026_06_02 = 1_780_000_000; // approx; recheck below
        let s = format_iso_utc(secs_2026_06_02);
        // Sanity: must start with "2026-" and end with "Z".
        assert!(s.starts_with("2026-"), "got {s}");
        assert!(s.ends_with("Z"), "got {s}");
    }

    #[test]
    fn encoding_format_table() {
        assert_eq!(encoding_format_of("block_index.json"), "application/json");
        assert_eq!(
            encoding_format_of("blocks/abcdef0123456789.bin"),
            "application/gzip"
        );
        assert_eq!(
            encoding_format_of("schemas/blocks_schema.txt"),
            "text/plain"
        );
        assert_eq!(encoding_format_of("README.md"), "text/markdown");
        assert_eq!(encoding_format_of("xyzzy"), "application/octet-stream");
    }

    #[test]
    fn dataset_identifier_shape() {
        let p = AssetParams {
            ring: 12,
            max_steps: 10,
            step: 1,
            free: true,
            target_block_bytes: 512,
            n_sequences: 16_751,
            oeis_a_number: Some("A316192"),
            produced_via: ProducedVia::InMemory,
        };
        assert_eq!(dataset_identifier(&p), "tilezz-rat-zz12-n10-free");
        let p2 = AssetParams { free: false, ..p };
        assert_eq!(dataset_identifier(&p2), "tilezz-rat-zz12-n10-onesided");
        // Stepped datasets present as their EFFECTIVE ring (the step is
        // an implementation detail). ZZ12 step-3 is the ZZ4 subring...
        let p3 = AssetParams { step: 3, ..p2 };
        assert_eq!(dataset_identifier(&p3), "tilezz-rat-zz4-n10-onesided");
        // ...and the real deployed case: ZZ14 step-2 is the odd ring ZZ7.
        let p4 = AssetParams {
            ring: 14,
            step: 2,
            ..p
        };
        assert_eq!(dataset_identifier(&p4), "tilezz-rat-zz7-n10-free");
    }

    /// Pipeline-shape selector emits different reproduce.sh content.
    #[test]
    fn reproduce_commands_per_pipeline() {
        let base = AssetParams {
            ring: 12,
            max_steps: 10,
            step: 1,
            free: true,
            target_block_bytes: 512,
            n_sequences: 16_751,
            oeis_a_number: None,
            produced_via: ProducedVia::InMemory,
        };
        let in_mem = reproduce_commands(&base);
        assert_eq!(in_mem.len(), 1, "in-memory should be one command");
        assert!(in_mem[0].contains("--mode dafsa-blocks"), "{:?}", in_mem);

        let streamed = reproduce_commands(&AssetParams {
            produced_via: ProducedVia::StreamingPipeline,
            ..base
        });
        assert!(
            streamed.len() >= 3,
            "streaming pipeline needs stream + merge + build, got {:?}",
            streamed
        );
        assert!(streamed[0].contains("--mode stream"));
        assert!(streamed[1].contains("--mode merge"));
        assert!(streamed[2].contains("--mode build"));
    }

    /// End-to-end fixture used by all the heavy tests below: write
    /// a small RatDafsa to a temp dir + copy schemas + write
    /// ro-crate-metadata.json, return the temp dir path.
    fn build_test_asset() -> std::path::PathBuf {
        use crate::stringmatch::dafsa::RatDafsa;
        // Tiny enumeration: a handful of canonical sequences over an
        // arbitrary alphabet. RatDafsa::from_rats sorts internally.
        let rats: Vec<Vec<i8>> = vec![
            vec![1, 2, 3],
            vec![1, 2, 4],
            vec![1, 3, 1],
            vec![2, 1, 5],
            vec![3, 0, -1, 2],
        ];
        let dafsa = RatDafsa::from_rats(rats.iter().map(|r| r.as_slice()));
        let dir = std::env::temp_dir().join(format!(
            "tilezz_rocrate_test_{}",
            std::process::id() as u64 * 1000
                + std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .subsec_nanos() as u64
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dafsa.write_blocks(&dir, 2).expect("write_blocks");
        let params = AssetParams {
            ring: 12,
            max_steps: 4,
            step: 1,
            free: true,
            target_block_bytes: 2,
            n_sequences: rats.len() as u64,
            oeis_a_number: Some("A316192"),
            produced_via: ProducedVia::InMemory,
        };
        let counts = SequenceCounts::from_rats(dafsa.iter());
        write_archival_extras(&dir, &params).expect("write_archival_extras");
        write_ro_crate(&dir, &params, &counts, None, None).expect("write_ro_crate");
        dir
    }

    /// Pre-commit hash check: every `sha256` recorded in
    /// `ro-crate-metadata.json` matches the actual file's SHA-256.
    /// If the writer ever drifts (e.g. starts hashing pre-gzip data
    /// instead of the on-disk bytes), this catches it on the first
    /// test run.
    #[test]
    fn sha256_matches_actual_files() {
        let dir = build_test_asset();
        let meta: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("ro-crate-metadata.json")).unwrap(),
        )
        .unwrap();
        let graph = meta["@graph"].as_array().expect("@graph array");
        let mut checked = 0;
        for entity in graph {
            let Some(id) = entity["@id"].as_str() else {
                continue;
            };
            // Skip non-file entities (the Dataset, License, Software,
            // CreateAction, OEIS contextual entity, etc.). File-shaped
            // entities are exactly the ones with our `sha256` property.
            let Some(expected) = entity["sha256"].as_str() else {
                continue;
            };
            let on_disk = sha256_hex(&dir.join(id)).expect("hash file");
            assert_eq!(
                expected, on_disk,
                "ro-crate sha256 for {id} does not match on-disk",
            );
            checked += 1;
        }
        // Manifest + 3 schemas + at least one block file.
        assert!(checked >= 5, "expected >= 5 hashed files, got {checked}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// `SOURCE_DATE_EPOCH` controls the CreateAction.endTime so two
    /// builds with the same epoch produce bit-identical
    /// ro-crate-metadata.json. Standard reproducible-builds
    /// contract.
    #[test]
    fn source_date_epoch_controls_end_time() {
        // SAFETY: `std::env::set_var` is unsafe to call from
        // multi-threaded code (Rust 2024 lint). This test sets it
        // briefly, then restores; running it serially via the test
        // harness's single-thread-per-test execution makes it safe.
        unsafe { std::env::set_var("SOURCE_DATE_EPOCH", "1780000000") };
        let dir1 = build_test_asset();
        let dir2 = build_test_asset();
        let m1 = std::fs::read_to_string(dir1.join("ro-crate-metadata.json")).unwrap();
        let m2 = std::fs::read_to_string(dir2.join("ro-crate-metadata.json")).unwrap();

        // Same endTime in both manifests, and it matches the
        // SOURCE_DATE_EPOCH formatting.
        let want = "2026-05-28T20:26:40Z";
        assert!(m1.contains(want), "epoch-derived endTime missing in run 1");
        assert!(m2.contains(want), "epoch-derived endTime missing in run 2");

        unsafe { std::env::remove_var("SOURCE_DATE_EPOCH") };
        let _ = std::fs::remove_dir_all(&dir1);
        let _ = std::fs::remove_dir_all(&dir2);
    }

    /// The required RO-Crate 1.2 entities are all present and the
    /// `conformsTo` chain is well-formed (every referenced `@id`
    /// resolves to another entity in the @graph).
    #[test]
    fn ro_crate_structure_is_valid() {
        let dir = build_test_asset();
        let meta: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("ro-crate-metadata.json")).unwrap(),
        )
        .unwrap();

        // Required top-level wiring.
        assert_eq!(
            meta["@context"], "https://w3id.org/ro/crate/1.2/context",
            "wrong @context"
        );
        let graph = meta["@graph"].as_array().expect("@graph array");
        let by_id: std::collections::HashMap<&str, &Value> = graph
            .iter()
            .filter_map(|e| Some((e["@id"].as_str()?, e)))
            .collect();

        // Metadata file descriptor.
        let descriptor = by_id
            .get("ro-crate-metadata.json")
            .expect("metadata descriptor entity present");
        assert_eq!(
            descriptor["conformsTo"]["@id"], "https://w3id.org/ro/crate/1.2",
            "descriptor conformsTo wrong"
        );
        assert_eq!(descriptor["about"]["@id"], "./");

        // Root Data Entity.
        let root = by_id.get("./").expect("root Dataset present");
        assert_eq!(root["@type"], "Dataset");
        assert!(root["name"].as_str().is_some(), "Dataset.name missing");
        assert!(root["datePublished"].as_str().is_some());
        assert!(root["hasPart"].as_array().is_some());
        assert!(root["mainEntity"]["@id"].as_str().is_some());
        assert_eq!(
            root["license"]["@id"], "https://creativecommons.org/licenses/by-sa/4.0/",
            "license wrong"
        );
        // distribution: a DataDownload for the manifest. With no base
        // URL the contentUrl is the relative path inside the crate.
        let dist = &root["distribution"];
        assert_eq!(dist["@type"], "DataDownload", "distribution @type wrong");
        assert_eq!(
            dist["contentUrl"], "block_index.json",
            "distribution.contentUrl wrong for relative (no base_url) build"
        );
        // creator points at the author's resolvable ORCID @id.
        let orcid = "https://orcid.org/0000-0002-5077-7497";
        let creators = root["creator"].as_array().expect("creator array");
        assert!(
            creators.iter().any(|c| c["@id"] == orcid),
            "creator does not reference the author ORCID {orcid}"
        );
        // ...and the matching Person entity carries it as @id + identifier.
        let person = by_id.get(orcid).expect("Person entity at ORCID @id");
        assert_eq!(person["@type"], "Person");
        assert_eq!(person["identifier"], orcid, "Person.identifier wrong");

        // License entity, software, create-action, OEIS contextual.
        assert!(by_id.contains_key("https://creativecommons.org/licenses/by-sa/4.0/"));
        assert!(by_id.contains_key("#tilezz"));
        let build = by_id.get("#build").expect("CreateAction present");
        assert_eq!(build["@type"], "CreateAction");
        assert!(build["endTime"].as_str().is_some());
        assert!(by_id.contains_key("https://oeis.org/A316192"));

        // conformsTo references resolve. We walk every File entity
        // and assert each conformsTo @id is itself an entity in the
        // graph.
        for entity in graph {
            let Some(c) = entity.get("conformsTo") else {
                continue;
            };
            let refs: Vec<&str> = if let Some(arr) = c.as_array() {
                arr.iter().filter_map(|v| v["@id"].as_str()).collect()
            } else {
                c["@id"].as_str().into_iter().collect()
            };
            for r in refs {
                // External URLs (RO-Crate spec itself, schema.org)
                // can be unresolved; only check intra-crate refs.
                if r.starts_with("http://") || r.starts_with("https://") {
                    continue;
                }
                assert!(
                    by_id.contains_key(r),
                    "conformsTo points at unresolved @id {r}",
                );
            }
        }

        // hasPart points at File entities that actually exist.
        let parts = root["hasPart"].as_array().unwrap();
        assert!(parts.len() >= 5, "hasPart has too few entries");
        for part in parts {
            let pid = part["@id"].as_str().expect("hasPart entry has @id");
            assert!(
                by_id.contains_key(pid),
                "hasPart referent {pid} not in @graph"
            );
            assert!(
                dir.join(pid).exists(),
                "hasPart referent {pid} missing on disk"
            );
        }

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The writer emits `tools/count.py`, the RO-Crate describes it,
    /// and `--print` prints the per-perimeter family series (free read
    /// off the index and cross-checked). Needs `python3` on PATH for
    /// the end-to-end half; skipped when absent.
    #[test]
    fn count_py_prints_family_terms() {
        let dir = build_test_asset();

        let script = dir.join("tools/count.py");
        assert!(script.exists(), "tools/count.py not emitted");
        let meta: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("ro-crate-metadata.json")).unwrap(),
        )
        .unwrap();
        let graph = meta["@graph"].as_array().expect("@graph array");
        let entity = graph
            .iter()
            .find(|e| e["@id"] == "tools/count.py")
            .expect("File entity for tools/count.py");
        assert!(
            entity["name"].as_str().is_some_and(|s| !s.is_empty()),
            "count.py entity lacks a human label"
        );

        let have_python = std::process::Command::new("python3")
            .arg("--version")
            .output()
            .is_ok_and(|o| o.status.success());
        if !have_python {
            eprintln!("skipping end-to-end half: python3 not on PATH");
            let _ = std::fs::remove_dir_all(&dir);
            return;
        }

        // The test asset's rats (see build_test_asset) are four
        // 3-angle sequences and one 4-angle sequence -> the free series
        // over perimeters 1..4 is 0,0,4,1.
        let output = std::process::Command::new("python3")
            .arg(&script)
            .arg(&dir)
            .arg("--print")
            .output()
            .expect("run count.py");
        assert!(
            output.status.success(),
            "count.py --print failed:\nstdout: {}\nstderr: {}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
        let stdout = String::from_utf8_lossy(&output.stdout);
        assert!(
            stdout.lines().any(|l| l.starts_with("free: 0,0,4,1")),
            "expected free series 0,0,4,1 in:\n{stdout}"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// End-to-end for the sub-family sequences: the writer emits a
    /// `variableMeasured` block with the seven named sequences, and
    /// the shipped `tools/count.py --verify` re-derives them from the
    /// DAFSA and agrees. The Rust-side check pins the structure +
    /// the free series (must match n_sequences by perimeter); the
    /// Python half closes the loop independently (skipped if python3
    /// is absent).
    #[test]
    fn variable_measured_emitted_and_verified() {
        let dir = build_test_asset();
        let meta: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("ro-crate-metadata.json")).unwrap(),
        )
        .unwrap();
        let graph = meta["@graph"].as_array().expect("@graph");
        let root = graph
            .iter()
            .find(|e| e["@id"] == "./")
            .expect("root Dataset");
        let vm = root["variableMeasured"]
            .as_array()
            .expect("variableMeasured array");
        let names: std::collections::BTreeSet<&str> =
            vm.iter().filter_map(|p| p["name"].as_str()).collect();
        for want in [
            "free",
            "oneSided",
            "achiral",
            "rotationSymmetric",
            "symmetric",
            "subring",
            "coset",
        ] {
            assert!(names.contains(want), "variableMeasured missing `{want}`");
        }
        // The test asset has 4 rats of length 3 and 1 of length 4
        // (see build_test_asset); free series is perimeter-indexed
        // from 1 up to max_indexed_length (4).
        let free = vm
            .iter()
            .find(|p| p["name"] == "free")
            .and_then(|p| p["value"].as_str())
            .expect("free series");
        assert_eq!(free, "0,0,4,1", "free series wrong: {free}");

        // Independent re-derivation via the shipped Python verifier.
        let have_python = std::process::Command::new("python3")
            .arg("--version")
            .output()
            .is_ok_and(|o| o.status.success());
        if have_python {
            let out = std::process::Command::new("python3")
                .arg(dir.join("tools/count.py"))
                .arg(&dir)
                .arg("--verify")
                .output()
                .expect("run count.py --verify");
            assert!(
                out.status.success(),
                "count.py --verify failed:\nstdout: {}\nstderr: {}",
                String::from_utf8_lossy(&out.stdout),
                String::from_utf8_lossy(&out.stderr),
            );
        } else {
            eprintln!("skipping count.py --verify half: python3 not on PATH");
        }

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Regression: when the asset's `max_steps` exceeds the deepest
    /// perimeter that actually closes a polygon (as happens for odd
    /// perimeters on an even-perimeter-only ring), the emitted
    /// `variableMeasured` series and `maxIndexedLength` must still track
    /// the data (the block manifest), NOT `max_steps`. Sizing them to
    /// `max_steps` appended spurious trailing zeros and made
    /// `count.py --verify` -- which sizes by the manifest's
    /// `max_indexed_length` -- report a false mismatch.
    #[test]
    fn variable_measured_handles_max_steps_above_deepest_perimeter() {
        use crate::stringmatch::dafsa::RatDafsa;
        // Deepest rat is length 4; declare max_steps = 7 (a 3-perimeter
        // gap with no closing polygon at perimeters 5..7).
        let rats: Vec<Vec<i8>> = vec![vec![1, 2, 3], vec![1, 2, 4], vec![3, 0, -1, 2]];
        let dafsa = RatDafsa::from_rats(rats.iter().map(|r| r.as_slice()));
        let dir = std::env::temp_dir().join(format!(
            "tilezz_rocrate_gap_{}",
            std::process::id() as u64 * 1000
                + std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .subsec_nanos() as u64
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dafsa.write_blocks(&dir, 2).expect("write_blocks");
        let params = AssetParams {
            ring: 12,
            max_steps: 7,
            step: 1,
            free: true,
            target_block_bytes: 2,
            n_sequences: rats.len() as u64,
            oeis_a_number: Some("A316192"),
            produced_via: ProducedVia::InMemory,
        };
        let counts = SequenceCounts::from_rats(dafsa.iter());
        write_archival_extras(&dir, &params).expect("write_archival_extras");
        write_ro_crate(&dir, &params, &counts, None, None).expect("write_ro_crate");

        let meta: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("ro-crate-metadata.json")).unwrap(),
        )
        .unwrap();
        let graph = meta["@graph"].as_array().expect("@graph");
        let root = graph
            .iter()
            .find(|e| e["@id"] == "./")
            .expect("root Dataset");

        // Series sized to the deepest closing perimeter (4), not max_steps (7).
        let free = root["variableMeasured"]
            .as_array()
            .unwrap()
            .iter()
            .find(|p| p["name"] == "free")
            .and_then(|p| p["value"].as_str())
            .expect("free series");
        assert_eq!(
            free, "0,0,2,1",
            "free series should stop at perimeter 4: {free}"
        );

        // maxIndexedLength reflects the data (4), not max_steps (7).
        let max_indexed = root["additionalProperty"]
            .as_array()
            .unwrap()
            .iter()
            .find(|p| p["name"] == "maxIndexedLength")
            .and_then(|p| p["value"].as_u64())
            .expect("maxIndexedLength");
        assert_eq!(max_indexed, 4, "maxIndexedLength should track the data");

        // The shipped verifier must agree (this is what regressed).
        let have_python = std::process::Command::new("python3")
            .arg("--version")
            .output()
            .is_ok_and(|o| o.status.success());
        if have_python {
            let out = std::process::Command::new("python3")
                .arg(dir.join("tools/count.py"))
                .arg(&dir)
                .arg("--verify")
                .output()
                .expect("run count.py --verify");
            assert!(
                out.status.success(),
                "count.py --verify failed on the gap asset:\nstdout: {}\nstderr: {}",
                String::from_utf8_lossy(&out.stdout),
                String::from_utf8_lossy(&out.stderr),
            );
        } else {
            eprintln!("skipping count.py --verify half: python3 not on PATH");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The produced `block_index.json` validates against the formal
    /// JSON Schema we ship at `schemas/block_index.schema.json`.
    /// This catches drift between the Rust `BlockManifest` struct
    /// and the schema spec.
    #[test]
    fn block_index_json_validates_against_schema() {
        let dir = build_test_asset();
        let schema: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("schemas/block_index.schema.json")).unwrap(),
        )
        .unwrap();
        let instance: Value =
            serde_json::from_str(&std::fs::read_to_string(dir.join("block_index.json")).unwrap())
                .unwrap();
        let validator = jsonschema::validator_for(&schema).expect("build validator");
        let errors: Vec<String> = validator
            .iter_errors(&instance)
            .map(|e| format!("{} at {}", e, e.instance_path()))
            .collect();
        assert!(
            errors.is_empty(),
            "block_index.json failed validation:\n  {}",
            errors.join("\n  ")
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Collect every File entity's (`@id`, `sha256`) pair from a crate,
    /// sorted by @id. The fingerprint of the dataset payload as the
    /// crate records it: a re-host must leave this identical.
    fn file_sha_pairs(crate_json: &Value) -> Vec<(String, String)> {
        let mut out: Vec<(String, String)> = crate_json["@graph"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|e| {
                let id = e["@id"].as_str()?;
                let sha = e["sha256"].as_str()?;
                Some((id.to_string(), sha.to_string()))
            })
            .collect();
        out.sort();
        out
    }

    /// Re-hosting rewrites ONLY the host-coordinate fields and leaves
    /// the blocks, their recorded sha256s, and all provenance
    /// byte-for-byte intact -- then round-trips back to the relative
    /// form. This is the correctness pin for the metadata-only re-host
    /// principle: blocks are content-addressed and provenance is fixed
    /// history, so only `identifier` / `url` / `distribution` / `sameAs`
    /// may move.
    #[test]
    fn rehost_rewrites_only_host_coordinates() {
        let dir = build_test_asset();

        // Snapshot the pre-rehost metadata + the on-disk block bytes.
        let before: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("ro-crate-metadata.json")).unwrap(),
        )
        .unwrap();
        let before_root = before["@graph"]
            .as_array()
            .unwrap()
            .iter()
            .find(|e| e["@id"] == "./")
            .unwrap()
            .clone();
        let before_files = file_sha_pairs(&before);
        // Read a real block file's bytes to prove the payload is never
        // touched.
        let blocks_dir = dir.join("blocks");
        let block_path = std::fs::read_dir(&blocks_dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .find(|p| p.extension().and_then(|s| s.to_str()) == Some("bin"))
            .expect("at least one block file");
        let block_bytes_before = std::fs::read(&block_path).unwrap();

        // Re-host to a Zenodo record with a DOI.
        rehost_ro_crate(
            &dir,
            Some("https://zenodo.org/records/123/files"),
            Some("10.5281/zenodo.123"),
        )
        .expect("rehost to zenodo");

        let after: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("ro-crate-metadata.json")).unwrap(),
        )
        .unwrap();
        let after_root = after["@graph"]
            .as_array()
            .unwrap()
            .iter()
            .find(|e| e["@id"] == "./")
            .unwrap()
            .clone();

        // Host coordinates moved to the DOI / Zenodo forms.
        assert_eq!(after_root["identifier"]["@type"], "PropertyValue");
        assert_eq!(after_root["identifier"]["propertyID"], "DOI");
        assert_eq!(after_root["identifier"]["value"], "10.5281/zenodo.123");
        assert_eq!(
            after_root["distribution"]["contentUrl"],
            "https://zenodo.org/records/123/files/block_index.json"
        );
        assert_eq!(after_root["sameAs"], "https://doi.org/10.5281/zenodo.123");
        assert_eq!(after_root["url"], "https://zenodo.org/records/123/files");

        // Everything else on the root is byte-identical to before.
        for key in [
            "datePublished",
            "version",
            "variableMeasured",
            "additionalProperty",
            "name",
            "description",
            "keywords",
            "creator",
            "license",
            "hasPart",
            "mainEntity",
            "subjectOf",
        ] {
            assert_eq!(
                after_root.get(key),
                before_root.get(key),
                "re-host changed `{key}` on the root Dataset",
            );
        }

        // Provenance entities untouched: #build endTime, #tilezz
        // softwareVersion.
        let find = |c: &Value, id: &str| -> Value {
            c["@graph"]
                .as_array()
                .unwrap()
                .iter()
                .find(|e| e["@id"] == id)
                .unwrap()
                .clone()
        };
        assert_eq!(
            find(&after, "#build")["endTime"],
            find(&before, "#build")["endTime"],
            "re-host changed the CreateAction endTime"
        );
        assert_eq!(
            find(&after, "#build"),
            find(&before, "#build"),
            "re-host changed the #build CreateAction"
        );
        assert_eq!(
            find(&after, "#tilezz")["softwareVersion"],
            find(&before, "#tilezz")["softwareVersion"],
            "re-host changed the tool softwareVersion"
        );

        // Every File entity's sha256 is preserved, and the block bytes
        // on disk are unchanged.
        assert_eq!(
            file_sha_pairs(&after),
            before_files,
            "re-host changed a File entity's sha256"
        );
        assert_eq!(
            std::fs::read(&block_path).unwrap(),
            block_bytes_before,
            "re-host changed block file bytes on disk"
        );

        // Round-trip: re-host back to the location-independent form
        // drops url / sameAs and restores the PropertyValue identifier
        // + relative distribution.
        rehost_ro_crate(&dir, None, None).expect("rehost to relative");
        let reverted: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join("ro-crate-metadata.json")).unwrap(),
        )
        .unwrap();
        let reverted_root = reverted["@graph"]
            .as_array()
            .unwrap()
            .iter()
            .find(|e| e["@id"] == "./")
            .unwrap()
            .clone();
        assert_eq!(
            reverted_root["identifier"]["propertyID"],
            "tilezz-dataset-id"
        );
        // The reverted identifier must carry the ORIGINAL dataset slug,
        // recovered from immutable #build provenance -- NOT the DOI left
        // by the intermediate Zenodo re-host. (Regression: recovering the
        // slug from the mutated `identifier` leaked the DOI here.)
        assert_eq!(
            reverted_root["identifier"]["value"], before_root["identifier"]["value"],
            "reverted identifier slug must match the original, not the DOI"
        );
        assert_eq!(
            reverted_root["distribution"]["contentUrl"],
            "block_index.json"
        );
        assert!(
            reverted_root.get("url").is_none(),
            "reverting to relative form must remove url"
        );
        assert!(
            reverted_root.get("sameAs").is_none(),
            "reverting to relative form must remove sameAs"
        );
        // Provenance still intact after the round-trip.
        assert_eq!(
            find(&reverted, "#build")["endTime"],
            find(&before, "#build")["endTime"],
            "round-trip changed the CreateAction endTime"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }
}
