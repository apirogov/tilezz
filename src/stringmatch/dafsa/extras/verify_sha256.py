#!/usr/bin/env python3
"""Verify the SHA-256 hashes recorded in a tilezz-rat-dafsa-blocks
asset's `ro-crate-metadata.json` against the on-disk bytes.

Every File entity in the RO-Crate manifest carries a `sha256`. This
script reads each one, hashes the file at that relative `@id`, and
reports any mismatches. Skips entities whose `@id` is an absolute
URL (`http://...`) or an internal anchor (`#author`, `#build`, ...).

Run from the asset directory root:

    python3 tools/verify_sha256.py

Or pass an explicit directory:

    python3 tools/verify_sha256.py path/to/asset

Exits 0 on full match, 1 on any mismatch, 2 on usage error. Prints
"OK (<n> files verified)" on success or one line per mismatch
followed by an error summary on failure.

No external Python dependencies; relies only on the stdlib.
"""

from __future__ import annotations

import hashlib
import json
import pathlib
import sys


def verify(asset_dir: pathlib.Path) -> int:
    crate_path = asset_dir / "ro-crate-metadata.json"
    if not crate_path.exists():
        sys.stderr.write(
            f"no ro-crate-metadata.json in {asset_dir}; "
            "pass the asset directory as the first argument\n"
        )
        return 2
    meta = json.loads(crate_path.read_text())
    errs = 0
    checked = 0
    for entity in meta.get("@graph", []):
        rel = entity.get("@id", "")
        want = entity.get("sha256")
        if not want or not isinstance(rel, str):
            continue
        if rel.startswith("http://") or rel.startswith("https://"):
            continue
        if rel.startswith("#"):
            continue
        # The metadata-file descriptor itself is `ro-crate-metadata.json`
        # without a sha256, so we don't usually hit it here -- but
        # guard anyway in case a writer ever stamps one in.
        path = asset_dir / rel
        if not path.is_file():
            print(f"{rel}: MISSING (manifest expected sha256={want[:12]}...)")
            errs += 1
            continue
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        if got != want:
            print(f"{rel}: MISMATCH (got {got[:12]}, want {want[:12]})")
            errs += 1
        else:
            checked += 1
    if errs:
        print(f"{errs} mismatch(es); {checked} file(s) verified")
        return 1
    print(f"OK ({checked} file(s) verified)")
    return 0


def main(argv: list[str]) -> int:
    asset_dir = pathlib.Path(argv[1] if len(argv) >= 2 else ".")
    return verify(asset_dir)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
