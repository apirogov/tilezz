#!/usr/bin/env python3
"""Per-perimeter family counts for a tilezz-rat-dafsa-blocks asset --
display and validation in one tool, two modes.

Every rat is a closed walk's cyclic turn sequence; bucketing the stored
rats by perimeter gives the seven derived sequences this asset realises
(index = perimeter, from 1):

  free               all rats (up to rotation and reflection)
  oneSided           up to rotation only (= 2*free - achiral)
  achiral            equal to its own mirror image
  rotationSymmetric  repetition factor > 1
  symmetric          achiral OR rotationSymmetric
  subring            all turns even (the order-n/2 sub-ring)
  coset              all turns odd (no straight segments; even perimeter)

Both modes print the same seven series (one `name: v1,v2,...` line each,
ready to read off as OEIS-style terms). They differ in how much they
trust:

  --print  (default, FAST): the `free` series is read straight off the
           DAFSA rank index (no decode) and cross-checked against the
           RO-Crate `variableMeasured`; the other six are echoed from
           `variableMeasured` for display only (NOT independently
           re-derived). Exits nonzero only if the free cross-check
           fails. Use this to read the terms, and on huge assets.

  --verify (SLOW): stream and decode every rat, independently re-derive
           ALL seven families, print them, and check each against
           `variableMeasured`. Different language + its own bucketing,
           so a mismatch points at an emitter bug or a corrupted asset.
           Exits nonzero on any mismatch. This is the metadata-integrity
           gate.

    python3 tools/count.py [path/to/asset] [--print | --verify]

Exits 0 on success, 1 on a (free / full) mismatch, 2 on usage error.
Depends only on the stdlib and the sibling tools/decode.py.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
from collections import Counter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from decode import BlockedDafsa  # noqa: E402  (sibling tool, stdlib-only)

FAMILIES = [
    "free",
    "oneSided",
    "achiral",
    "rotationSymmetric",
    "symmetric",
    "subring",
    "coset",
]


def _rot_min(t: tuple) -> tuple:
    return min(t[i:] + t[:i] for i in range(len(t)))


def free_from_index(bd: BlockedDafsa) -> tuple[dict[int, int], int]:
    """Per-perimeter `free` counts read off the DAFSA rank index -- one
    state read per root edge, NO rat decoding. Returns (counts, root_count)."""
    _, root_count, _, root_edges = bd._state(0)
    counts: dict[int, int] = {}
    for label, target in sorted(root_edges):
        _, count, _, _ = bd._state(target)
        counts[label] = count
    return counts, root_count


def recompute_all(bd: BlockedDafsa) -> dict[str, Counter]:
    """Independently re-derive all seven families by decoding every rat."""
    free = Counter(); subring = Counter(); coset = Counter()
    achiral = Counter(); rotsym = Counter(); symm = Counter()
    for rat in bd.iter_rats():
        t = tuple(rat); L = len(t)
        free[L] += 1
        if all(a % 2 == 0 for a in t):
            subring[L] += 1
        if all(a % 2 != 0 for a in t):
            coset[L] += 1
        ach = _rot_min(t) == _rot_min(t[::-1])
        rot = any(t[d:] + t[:d] == t for d in range(1, L))
        if ach:
            achiral[L] += 1
        if rot:
            rotsym[L] += 1
        if ach or rot:
            symm[L] += 1
    one_sided = {n: 2 * free[n] - achiral[n] for n in free}
    return {
        "free": free, "oneSided": one_sided, "achiral": achiral,
        "rotationSymmetric": rotsym, "symmetric": symm,
        "subring": subring, "coset": coset,
    }


def load_emitted(asset_dir: pathlib.Path) -> dict[str, list[int]]:
    crate = json.loads((asset_dir / "ro-crate-metadata.json").read_text())
    out: dict[str, list[int]] = {}
    for e in crate.get("@graph", []):
        if e.get("@id") == "./" or e.get("@type") == "Dataset":
            for pv in (e.get("variableMeasured") or []):
                out[pv["name"]] = [int(x) for x in pv["value"].split(",")]
    return out


def to_series(counts, maxp: int) -> list[int]:
    return [counts.get(n, 0) for n in range(1, maxp + 1)]


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Per-perimeter family counts (print / verify).")
    ap.add_argument("asset_dir", nargs="?", default=".")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--print", dest="do_print", action="store_true",
                      help="fast: free off the index (cross-checked), others echoed from metadata (default)")
    mode.add_argument("--verify", action="store_true",
                      help="slow: decode all rats, re-derive ALL families, check vs metadata")
    args = ap.parse_args(argv[1:])

    d = pathlib.Path(args.asset_dir)
    if not (d / "block_index.json").exists() or not (d / "ro-crate-metadata.json").exists():
        sys.stderr.write(f"need block_index.json + ro-crate-metadata.json in {d}\n")
        return 2

    manifest = json.loads((d / "block_index.json").read_text())
    maxp = manifest["max_indexed_length"]
    n_sequences = manifest["n_sequences"]
    emitted = load_emitted(d)
    if not emitted:
        sys.stderr.write("no variableMeasured in ro-crate-metadata.json\n")
        return 2
    bd = BlockedDafsa(d)

    if args.verify:
        rec = recompute_all(bd)
        values = {f: to_series(rec[f], maxp) for f in FAMILIES}
        bad = 0
        for f in FAMILIES:
            ok = values[f] == emitted.get(f)
            bad += 0 if ok else 1
            print(f"{f}: " + ",".join(map(str, values[f])) + ("" if ok else "   *** MISMATCH vs metadata ***"))
        total = sum(values["free"])
        print(f"# n_sequences {n_sequences}; free total {total}; mode verify (all re-derived from rats)")
        if bad or total != n_sequences:
            print(f"# {bad} sequence(s) disagree with metadata"
                  + ("" if total == n_sequences else "; free total != n_sequences"))
            return 1
        print(f"# OK ({len(FAMILIES)} sequences re-derived and verified to perimeter {maxp})")
        return 0

    # default: --print (fast)
    free_idx, root_count = free_from_index(bd)
    free_series = to_series(free_idx, maxp)
    free_ok = free_series == emitted.get("free") and root_count == n_sequences
    print(f"free: " + ",".join(map(str, free_series))
          + ("   [index-derived, matches metadata]" if free_ok else "   *** free cross-check FAILED ***"))
    for f in FAMILIES[1:]:
        print(f"{f}: " + ",".join(map(str, emitted.get(f, []))) + "   [echoed from metadata]")
    print(f"# n_sequences {n_sequences}; free total {sum(free_series)}; "
          f"mode print (free independent off the index, others echoed -- use --verify for a full re-derivation)")
    if not free_ok:
        print("# free cross-check FAILED (index vs metadata)")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
