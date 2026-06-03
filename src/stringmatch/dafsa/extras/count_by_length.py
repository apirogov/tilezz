#!/usr/bin/env python3
"""Print per-exact-length sequence counts (OEIS-style terms) for a
tilezz-rat-dafsa-blocks asset.

The asset stores each rat length-prefixed (`rat_schema.txt`:
stored_sequence = [len, rat...]), so every rat with exactly L
angles lives in the sub-DAFSA behind the root edge labeled L. And
every DAFSA state already records `count`, the number of accepted
sequences reachable from it (the rank/select index used for id
lookups). Combining the two, the full exact-perimeter histogram is
just the root's edge list paired with each target state's `count`
-- no decoding, no traversal beyond one state read per length.

Run from the asset directory root (or pass it as the argument):

    python3 tools/count_by_length.py [path/to/asset]

Output: one `<length> <count>` line per root edge, lengths
ascending, plus `#`-prefixed summary lines on stdout. Compare the
counts against the OEIS sequence named in the asset's RO-Crate
where published (note OEIS entries count perimeter EXACTLY n,
while the manifest's `n_sequences` is cumulative over <= n).

Exits 0 on success with consistent totals, 1 if the per-length
counts do not sum to `n_sequences`, 2 on usage error.

Depends only on the stdlib and its sibling `tools/decode.py`
(imported for the block reader).
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from decode import BlockedDafsa  # noqa: E402  (sibling tool, stdlib-only)


def main(argv: list[str]) -> int:
    asset_dir = pathlib.Path(argv[1] if len(argv) >= 2 else ".")
    if not (asset_dir / "block_index.json").exists():
        sys.stderr.write(
            f"no block_index.json in {asset_dir}; "
            "pass the asset directory as the first argument\n"
        )
        return 2
    bd = BlockedDafsa(asset_dir)
    _, root_count, _, root_edges = bd._state(0)
    total = 0
    for label, target in sorted(root_edges):
        _, count, _, _ = bd._state(target)
        print(f"{label} {count}")
        total += count
    n_sequences = bd.manifest["n_sequences"]
    print(f"# total {total}; manifest n_sequences {n_sequences}; root count {root_count}")
    if total != n_sequences or root_count != n_sequences:
        print("# MISMATCH: per-length counts do not sum to n_sequences")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
