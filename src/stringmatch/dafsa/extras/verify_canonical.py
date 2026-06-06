#!/usr/bin/env python3
"""Verify every stored sequence is the dihedral-canonical CCW
representative -- an INDEPENDENT re-derivation of the enumerator's
canonical form.

Picture each rat as a closed walk's cyclic turn sequence. Two walks are
the same *free* rat when one is a rotation or a reflection of the other,
and the dataset is meant to store exactly one representative per class:
the lexicographically smallest sequence over all cyclic rotations of the
sequence AND of its reverse, oriented counter-clockwise. This script
recomputes that representative from scratch (pure list arithmetic, no
ring geometry, no shared code with the Rust enumerator) and checks that
every stored sequence already equals it. Two things follow if all pass:

  * no NON-CANONICAL sequence is stored (the enumerator's
    canonicalization produced a true fixed point), and
  * since a DAFSA stores each string at most once, no two stored
    sequences can be the same class -- a correct canonical form maps a
    class to one string, so a second member would have collapsed onto
    it. This closes the "non-canonical duplicate" failure mode that a
    plain sha256 / count check cannot see.

Conventions (matching the enumerator, verified against the source):
  * CCW orientation  <=>  signed turn-sum > 0  (chirality = sign of the
    turn-sum); a closed CCW polygon turns through one positive full turn.
  * reflection within a CCW class is plain sequence REVERSAL (not
    reverse-and-negate) -- after CCW-normalization the two enantiomers
    of a shape are related by reversal.

Run from the asset directory root, or pass it as the first argument:

    python3 tools/verify_canonical.py [path/to/asset]

For very large assets, sample instead of a full sweep:

    python3 tools/verify_canonical.py --stride 1000   # every 1000th rat
    python3 tools/verify_canonical.py --max 5000000    # stop after N

Exits 0 if every checked sequence is canonical CCW, 1 on the first
batch of violations (a few are printed), 2 on usage error. Depends only
on the stdlib and its sibling tools/decode.py.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from decode import BlockedDafsa  # noqa: E402  (sibling tool, stdlib-only)


def canonical_form(seq: list[int]) -> list[int]:
    """Lex-min over all cyclic rotations of `seq` and of `reversed(seq)`."""
    n = len(seq)
    best: list[int] | None = None
    for base in (seq, seq[::-1]):
        doubled = base + base
        for k in range(n):
            rot = doubled[k : k + n]
            if best is None or rot < best:
                best = rot
    return best if best is not None else []


def is_canonical_ccw(seq: list[int]) -> tuple[bool, str]:
    """(ok, reason). A stored rat must be CCW (turn-sum > 0) and equal to
    its own dihedral-canonical form."""
    if sum(seq) <= 0:
        return False, f"not CCW (turn-sum {sum(seq)} <= 0)"
    if seq != canonical_form(seq):
        return False, "not the dihedral-canonical (lex-min rotation/reverse) form"
    return True, ""


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Verify canonical CCW form of stored rats.")
    ap.add_argument("asset_dir", nargs="?", default=".")
    ap.add_argument("--stride", type=int, default=1, help="check every Nth rat (sampling)")
    ap.add_argument("--max", type=int, default=0, help="stop after N checks (0 = no limit)")
    args = ap.parse_args(argv[1:])

    asset_dir = pathlib.Path(args.asset_dir)
    if not (asset_dir / "block_index.json").exists():
        sys.stderr.write(
            f"no block_index.json in {asset_dir}; pass the asset directory\n"
        )
        return 2

    bd = BlockedDafsa(asset_dir)
    checked = 0
    violations = 0
    for i, seq in enumerate(bd.iter_rats()):
        if args.stride > 1 and (i % args.stride) != 0:
            continue
        ok, reason = is_canonical_ccw(seq)
        checked += 1
        if not ok:
            violations += 1
            if violations <= 10:
                print(f"NON-CANONICAL: {seq} -- {reason}")
        if args.max and checked >= args.max:
            break

    mode = "full" if args.stride == 1 and not args.max else f"sampled (stride={args.stride}, max={args.max})"
    if violations:
        print(f"{violations} violation(s) in {checked} checked rat(s) [{mode}]")
        return 1
    print(f"OK ({checked} rat(s) verified canonical CCW [{mode}])")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
