# Maintenance

Operational handbook for the moving parts that are not obvious from
the code alone. Currently one section: the RatDB dataset lifecycle.

## RatDB datasets

The big picture: enumerated rat datasets (tilezz-rat-dafsa-blocks
format) do not live in this repo. Each dataset is one orphan branch
of <https://github.com/apirogov/tilezz-ratdb>, carrying a complete,
self-describing RO-Crate at the branch root: `block_index.json`
(the wire manifest), `ro-crate-metadata.json`, `schemas/`, `tools/`,
`README.md`, and the content-addressed `blocks/<sha256>.bin` files.

This repo only pins them: `web/ratdb/datasets.json` lists
`{name, repo, ref}` per dataset, where `ref` is a tilezz-ratdb
commit sha. At deploy time, `.github/workflows/pages.yml` copies the
pinned ref's metadata (everything except `blocks/`) into
`web/ratdb/data/<name>/` and injects a `block_base_url` derived from the
same ref, so the explorer lazy-fetches blocks straight from
`raw.githubusercontent.com/apirogov/tilezz-ratdb/<ref>/blocks/`.

Why raw.githubusercontent.com: it is the only zero-cost GitHub
endpoint that serves with `Access-Control-Allow-Origin: *` (GitHub
Release downloads send no CORS headers at all and cannot be fetched
from a browser; Git LFS works but meters bandwidth). Pinning a
commit sha makes the served bytes immutable: updating a dataset
moves the branch, but old pins stay reachable as ancestors.

### 1. Getting a copy of a dataset

Each dataset branch is self-contained. Clone just that branch:

```sh
git clone --branch zz12-n14-free --single-branch --depth 1 \
    https://github.com/apirogov/tilezz-ratdb
cd tilezz-ratdb
python3 tools/verify_sha256.py    # must print: OK (<n> files verified)
python3 tools/decode.py > rats.txt  # plain-text export, no Rust needed
```

Alternatives, no git required:

```sh
# whole dataset as a tarball (any ref: branch name or commit sha)
curl -sL https://api.github.com/repos/apirogov/tilezz-ratdb/tarball/zz12-n14-free | tar -xz

# single file
curl -sLO https://raw.githubusercontent.com/apirogov/tilezz-ratdb/zz12-n14-free/block_index.json
```

The branch list (= dataset list) is indexed in the tilezz-ratdb
`main` README.

### 2. Computing a dataset with rat_enum

`rat_enum` is the producer for everything in tilezz-ratdb. Always
build it from a clean checkout of a commit on `origin/main` (see
section 5) and with `SOURCE_DATE_EPOCH` exported.

Scout first: `--mode bench` enumerates without writing anything,
so you can check counts and runtime before committing to a layout
(`--stats` adds a per-length breakdown). Expect each `+1` in `-n`
to multiply runtime by ~3-10x depending on ring.

Always pass for our published datasets:

- `--free` -- dihedral-canonical output (one representative per
  chiral pair). This is a SEMANTIC choice, not a tuning flag: it
  decides which counting sequence the dataset realises (e.g. ZZ12
  free matches OEIS A316192). All published datasets so far are
  free; do not mix conventions within a dataset.
- `--threads 0` -- use all cores.

Optional accelerators (pure pruning, never change the output;
combine freely): `--mod-prune`, `--closure-key-prune`. Use
`--paranoid` for an (expensive) self-checking validation run when
touching enumeration internals.

**Single-step path** (n small enough that the full rat set fits in
RAM -- fine for ZZ12 up to n~13 on a workstation):

```sh
./target/release/rat_enum --ring 12 -n 12 --free \
    --mode dafsa-blocks --threads 0 -o /tmp/zz12_n12_free
```

Enumerates in memory, builds the DAFSA, writes the complete asset
directory (manifest + blocks + schemas + RO-Crate + tools) in one
go.

**Multi-step streaming path** (bounded memory, for n >= 14-ish;
this is the n=16 plan). Three stages over one output directory,
with the SAME `--ring/-n/--step/--free` flags repeated each time
(they are recorded verbatim in `certificate.json` and
cross-checked):

```sh
D=/tmp/zz12_n16_free
# 1. parallel DFS, per-thread sorted runs (~16 MB/worker, no HashSet)
./target/release/rat_enum --ring 12 -n 16 --free --threads 0 \
    --mode stream -o $D
# 2. k-way merge runs -> unique.bin + certificate.json (BLAKE3 + counts)
./target/release/rat_enum --ring 12 -n 16 --free --mode merge -o $D
# 3. stream unique.bin through the DAFSA builder -> $D/dafsa/
./target/release/rat_enum --ring 12 -n 16 --free --mode build -o $D
```

The blocked asset lands in `$D/dafsa/` -- publish THAT directory,
not `$D` itself (`$D` also holds the intermediate `runs/`,
`unique.bin` and `certificate.json`, which are not part of the
asset). RO-Crate provenance and `--oeis-a-number` are honoured by
the asset-writing modes (`dafsa-blocks` and `build`).

**Deploying vs. local builds**: for a dataset you intend to PUBLISH,
pass `--base-url <public URL of the dataset directory>` to the
asset-writing mode (`dafsa-blocks` or `build`). That is the URL the
served directory will live at -- for our setup the
raw.githubusercontent.com tilezz-ratdb path, e.g.
`https://raw.githubusercontent.com/apirogov/tilezz-ratdb/<ref>` --
and it sets the RO-Crate `distribution.contentUrl`, `identifier`,
and `url` to the real location instead of relative paths. Omit
`--base-url` for local / test builds (the metadata then carries the
location-independent relative forms, which is what the archival
tilezz-ratdb branch keeps). If you already have a minted DOI at
compute time you can also pass `--doi <DOI>`; usually the DOI comes
later, so see section 6 to add it after the fact without recomputing.

**Multi-host splitting** (if one machine is not enough):
`--mode list-seeds` prints DFS prefixes at `--split-depth`; farm
each prefix out with `--seed a,b,c` on separate hosts, then feed
all resulting runs through the same merge + build stages.

**Block size**: `--target-block-bytes` counts UNCOMPRESSED bytes
per block (the `.bin` files are additionally gzipped, roughly 3x
smaller on disk). The 1 MiB default is the intended sweet spot for
HTTP-served assets; the n=10 probe used 8192 deliberately to force
a multi-block layout out of a tiny dataset. Don't do that for real
datasets.

**"Upgrading" to a higher n**: there is no incremental path -- a
larger perimeter bound is a fresh enumeration from scratch (the
DAFSA and the block set change globally anyway). What carries over
is everything in this section, plus section 4 for re-publishing.

### 3. Adding a new dataset (new ring / new parameters)

Compute the asset per section 2, then publish it as a new orphan
branch:

```sh
# 0. computed per section 2 into /tmp/zz24_n8_free; sanity-check it
python3 /tmp/zz24_n8_free/tools/verify_sha256.py /tmp/zz24_n8_free

# 1. publish as orphan branch of tilezz-ratdb
cd ../tilezz-ratdb
git checkout --orphan zz24-n8-free && git rm -rq --cached . && git clean -fdq
cp -r /tmp/zz24_n8_free/. .
git add -A && git commit -m "data: ZZ24 n=8 free (<count> rats)"
git push origin zz24-n8-free
git rev-parse zz24-n8-free          # note this sha: it is the pin

# 2. index it on tilezz-ratdb main (add a row to the README table)
git checkout main   # edit README.md, commit, push

# 3. pin it in tilezz: add an entry to web/ratdb/datasets.json
#    {"name": "zz24_n8_free", "repo": "apirogov/tilezz-ratdb",
#     "ref": "<sha from step 1>"}
#    then PR + merge as usual.
```

Explorer discovery is automatic from there: the deploy copies the
metadata into `web/ratdb/data/zz24_n8_free/`, `build_web_rocrate` scans
`web/ratdb/data/*/ro-crate-metadata.json` into the top-level RO-Crate,
and the explorer JS drives ring/dataset selection off that single
fetch. Nothing else to wire up.

Note: pass `--oeis-a-number` (at compute time) only when the count
sequence actually matches an OEIS entry; it lands in the RO-Crate
as a contextual entity.

### 4. Updating an existing dataset (pushing n upward)

Recompute per section 2 (a higher n is always a fresh enumeration),
then replace the branch content in place:

```sh
cd tilezz-ratdb
git checkout zz12-n14-free
find . -mindepth 1 -maxdepth 1 -not -name .git -exec rm -rf {} +
cp -r /tmp/zz12_n16_free/dafsa/. .
git add -A && git commit -m "data: extend to n=16 (<count> rats)"
git push origin zz12-n14-free        # rename branch first if the name encodes n
git rev-parse HEAD                   # the new pin
```

Then bump the `ref` (and `name`, if renamed) in `web/ratdb/datasets.json`
and merge. Two properties make this safe:

- Blocks are content-addressed: any block whose bytes did not
  change keeps its filename, and git stores it once across the
  branch history.
- The deployed site keeps serving the OLD pin until the
  datasets.json bump deploys -- old pins are ancestors of the new
  branch tip, so the sha-pinned raw URLs never break mid-rollout.

If the branch name encodes the perimeter bound (zz12-n14-free),
prefer creating the new branch name (zz12-n16-free) via the
section-3 flow and retiring the old one from the README index; the
old branch can be deleted once nothing pins it (check
web/ratdb/datasets.json history if unsure).

### 5. Reproducibility: keeping the hashes wired up

Every dataset records its provenance in three places (README.md,
ro-crate-metadata.json `CreateAction`, tools/reproduce.sh), all
derived from two inputs:

- The producing commit: `build.rs` bakes `git rev-parse HEAD` into
  the `rat_enum` binary at COMPILE time (`TILEZZ_GIT_COMMIT`).
- `SOURCE_DATE_EPOCH`: pins the `CreateAction.endTime` so reruns
  are bit-identical. We use `1780000000` by convention.

The reproducibility contract: `bash tools/reproduce.sh` in a
dataset must be able to `git clone` the tilezz repo, check out the
recorded commit, rebuild, and regenerate the dataset
file-for-file. That holds only if ALL of the following are true at
build time:

1. Working tree is clean (build.rs records bare HEAD and will
   silently overclaim if you build with uncommitted changes).
2. HEAD is a commit reachable from `origin/main` -- NOT an
   unmerged PR branch. GitHub's "Squash and merge" AND "Rebase and
   merge" both rewrite commits to new hashes (rebase-merge even
   when fast-forward is possible), which leaves the recorded
   commit dangling after the merge. This bit us in the n=10 probe.
   If a dataset must be cut from in-flight work, merge first
   (fast-forward from CLI or a UI merge commit preserves hashes),
   or regenerate the dataset after the merge -- blocks are
   deterministic, so only the three provenance files change.
3. `rat_enum` was rebuilt after the last commit (cargo re-runs
   build.rs when HEAD moves, so a plain `cargo build` suffices --
   just do not skip it).

Verification checklist after publishing (cheap, do all of it):

```sh
# recorded commit exists upstream (not just locally / in a PR ref)
git fetch origin && git merge-base --is-ancestor <recorded-sha> origin/main && echo ok

# dataset internally consistent
python3 tools/verify_sha256.py

# counts are sane: per-exact-length terms must match the OEIS
# sequence where published and any previously computed overlap.
# (Careful: OEIS A316192 counts perimeter EXACTLY n, while the
# manifest's n_sequences is cumulative over <= n.)
python3 tools/count_by_length.py

# live wiring: deployed manifest points at the right pin and the
# blocks serve with CORS
curl -s https://ratdb.app.pirogov.de/data/<name>/block_index.json | grep block_base_url
curl -sI -H "Origin: https://ratdb.app.pirogov.de" \
    "<block_base_url><first-sha256>.bin" | grep -i access-control

# full reproduce round-trip (slow but definitive)
SOURCE_DATE_EPOCH=1780000000 bash tools/reproduce.sh
```

### 6. Re-hosting a dataset to a new host / Zenodo (metadata only, no recompute)

A dataset's blocks are content-addressed and its provenance (the
producing commit, the `CreateAction.endTime`, the reproduce recipe,
every File `sha256`) is fixed history. Re-hosting the same dataset --
mirroring the raw.githubusercontent copy to a Zenodo deposit with a
DOI, or otherwise switching where the bytes are served -- therefore
changes ONLY the "host coordinates": `identifier`, `url`,
`distribution.contentUrl`, and (with a DOI) `sameAs`.

So we do NOT recompute and we do NOT re-run the emitter (its
`TILEZZ_GIT_COMMIT` and build time would differ and corrupt the
recorded provenance). Instead, surgically patch just those fields:

```sh
rat_enum --mode rehost -o <dataset-dir> \
    --base-url <new file base URL> [--doi <minted DOI>]
```

This rewrites only `identifier` / `url` / `distribution` / `sameAs`
in `<dataset-dir>/ro-crate-metadata.json` and leaves the blocks,
their recorded sha256s, the provenance, and `variableMeasured`
byte-for-byte untouched (a `git diff` after the patch shows only
those host fields). `--mode rehost` does no enumeration, so it does
not need `--ring` / `-n`. Re-run `verify_sha256.py` afterwards to
confirm the payload still matches the (unchanged) hashes.

Zenodo flow:

1. Create the Zenodo deposit and upload the dataset (the whole
   directory's files, or a single tarball of it). Zenodo mints a DOI
   like `10.5281/zenodo.123`.
2. Re-host with that DOI and the record's file base URL:
   - if Zenodo serves the files individually, use their base URL as
     `--base-url` (so `<base>/block_index.json` resolves), e.g.
     `--base-url https://zenodo.org/records/123/files`;
   - if you uploaded a single tarball instead, point
     `--base-url` at the record's file URL so `contentUrl` resolves
     to the artifact you actually published.
3. The identifier becomes a DOI PropertyValue, `sameAs` becomes
   `https://doi.org/<doi>`, and `url` / `distribution.contentUrl`
   point at the Zenodo location.

The same command covers any host switch, not just Zenodo: pass the
new `--base-url` (and `--doi` if one applies). Passing NEITHER flag
(`rat_enum --mode rehost -o <dataset-dir>`) reverts the crate to the
location-independent relative forms and drops `url` / `sameAs` --
useful to clear a stale absolute URL, e.g. before re-archiving the
canonical tilezz-ratdb branch copy, which deliberately stays
location-independent (see the note at the end of this document).

### Exactness / overflow verification gate

The geometry is exact by construction: `cell_floor` (the only input
to the self-intersection grid) is exact-by-default, so the f64 fast
path can no longer mis-bucket a boundary point and let a non-simple
polygon through. The exact-sign helpers run in i128; reachable
coordinate magnitudes (~hundreds even at n~20-30) sit ~36 orders of
magnitude below the i128 ceiling, and an overflow-checked sweep of
every ring (ZZ8/10/16/20/24/32/60 nested-sqrt + ZZ14/18 cubic) shows
no wrap.

For a count you intend to PUBLISH or SUBMIT (a new dataset, or an
OEIS extension/correction), run the producing enumeration once more
with integer overflow-checks on, as a belt-and-suspenders against a
silent i128 wrap at a magnitude beyond what's been measured:

```sh
RUSTFLAGS="-C overflow-checks=on" cargo build --release --bin rat_enum --features cli
# re-run the exact producing invocation; identical count, panics loudly on any wrap
RUSTFLAGS="-C overflow-checks=on" ./target/release/rat_enum --ring <r> -n <n> --free ...
```

This is cheap relative to the value (one extra run on the artifact
you stake a claim on). The deep nested-sqrt rings (ZZ32/60) have no
external count oracle yet, so for those a small-n independent brute
check (see `docs/oeis-A316200-correction/zz10_independent.py` as a
template) is also worth doing before trusting their counts.

The archival manifest (in tilezz-ratdb) deliberately has NO
`block_base_url` -- it stays location-independent, and
`tools/decode.py` works offline against the sibling `blocks/`
directory. Only the deploy-time copy in `web/ratdb/data/<name>/` gets
the URL injected; if the two ever disagree on anything else, the
tilezz-ratdb branch is the canonical one.
