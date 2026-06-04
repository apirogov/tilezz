#!/usr/bin/env python3
"""Recompute the per-perimeter sub-family sequences directly from this
DAFSA and verify them against the `variableMeasured` block the writer
emitted into `ro-crate-metadata.json`.

This is an INDEPENDENT re-derivation (different language, its own
bucketing) of the same seven sequences the Rust writer computed, so a
disagreement points at a bug in the emitter (or a corrupted asset).

Sequences (index = perimeter, from 1; mirror the RO-Crate descriptions):
  free               -- all rats (up to rotation and reflection)
  oneSided           -- up to rotation only (= 2*free - achiral)
  achiral            -- equal to its own mirror image
  rotationSymmetric  -- repetition factor > 1
  symmetric          -- achiral OR rotationSymmetric
  subring            -- all turns even (order-n/2 sub-ring)
  coset              -- all turns odd (no straight segments; even perimeter)

Run from the asset directory root (or pass it):

    python3 tools/verify_counts.py [path/to/asset]

Exits 0 if every sequence matches, 1 on any mismatch, 2 on usage error.
Depends only on the stdlib and the sibling `tools/decode.py`.
"""

from __future__ import annotations

import json
import pathlib
import sys
from collections import Counter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from decode import BlockedDafsa  # noqa: E402  (sibling tool, stdlib-only)


def rot_min(t: tuple) -> tuple:
    return min(t[i:] + t[:i] for i in range(len(t)))


def recompute(bd: BlockedDafsa):
    free = Counter(); subring = Counter(); coset = Counter()
    achiral = Counter(); rotsym = Counter(); symm = Counter()
    for rat in bd.iter_rats():
        t = tuple(rat); L = len(t)
        free[L] += 1
        if all(a % 2 == 0 for a in t):
            subring[L] += 1
        if all(a % 2 != 0 for a in t):
            coset[L] += 1
        ach = rot_min(t) == rot_min(t[::-1])
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


def series(counts, max_perim: int) -> str:
    return ",".join(str(counts.get(n, 0)) for n in range(1, max_perim + 1))


def main(argv: list[str]) -> int:
    d = pathlib.Path(argv[1] if len(argv) >= 2 else ".")
    if not (d / "block_index.json").exists() or not (d / "ro-crate-metadata.json").exists():
        sys.stderr.write(f"need block_index.json + ro-crate-metadata.json in {d}\n")
        return 2
    max_perim = json.loads((d / "block_index.json").read_text())["max_indexed_length"]
    crate = json.loads((d / "ro-crate-metadata.json").read_text())
    emitted: dict[str, str] = {}
    for e in crate.get("@graph", []):
        if e.get("@id") == "./" or e.get("@type") == "Dataset":
            for pv in (e.get("variableMeasured") or []):
                emitted[pv["name"]] = pv["value"]
    if not emitted:
        sys.stderr.write("no variableMeasured in ro-crate-metadata.json\n")
        return 2

    recomputed = recompute(BlockedDafsa(d))
    bad = 0
    for name, counts in recomputed.items():
        want = series(counts, max_perim)
        got = emitted.get(name)
        if got != want:
            bad += 1
            print(f"{name}: MISMATCH\n  emitted:    {got}\n  recomputed: {want}")
    if bad:
        print(f"{bad} sequence(s) disagree")
        return 1
    print(f"OK ({len(recomputed)} sequences verified to perimeter {max_perim})")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
