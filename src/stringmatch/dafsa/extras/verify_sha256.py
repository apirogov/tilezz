#!/usr/bin/env python3
"""Verify a tilezz-rat-dafsa-blocks asset against its
`ro-crate-metadata.json` -- paranoid mode.

Three independent checks, any of which fails the run:

1. Recorded -> disk: every File entity in the manifest carries a
   `sha256`; each is re-hashed from the file at its relative `@id`
   and compared. (Skips absolute-URL and internal-anchor `@id`s.)

2. Disk -> recorded: every file actually present in the asset
   directory must be accounted for in the manifest. An UNEXPECTED
   file (one on disk but not hashed in the manifest) fails the run --
   this catches injected/edited content and stray build artifacts
   that a recorded->disk pass alone would silently ignore.
   `ro-crate-metadata.json` itself (the manifest) and transient
   `__pycache__/` byproducts of running these tools are exempt.

3. Provenance sanity: the producing-tool `softwareVersion` (the build
   commit) is reported, and a loud WARNING is printed if it is
   `unknown` (built outside a git checkout) or `-dirty` (built from an
   unclean tree) -- i.e. when the recorded commit does not fully
   describe the binary that produced the asset.

Run from the asset directory root, or pass it as the first argument:

    python3 tools/verify_sha256.py [path/to/asset]

Exits 0 on full match, 1 on any mismatch / missing / unexpected file,
2 on usage error. No external dependencies; stdlib only.
"""

from __future__ import annotations

import hashlib
import json
import pathlib
import sys


def _recorded_rel_paths(meta: dict) -> dict[str, str]:
    """Map of relative-path -> recorded sha256 for local File entities."""
    out: dict[str, str] = {}
    for entity in meta.get("@graph", []):
        rel = entity.get("@id", "")
        want = entity.get("sha256")
        if not want or not isinstance(rel, str):
            continue
        if rel.startswith(("http://", "https://", "#")):
            continue
        out[rel] = want
    return out


def _provenance_warnings(meta: dict) -> list[str]:
    warns: list[str] = []
    for entity in meta.get("@graph", []):
        if entity.get("@id") == "#tilezz":
            ver = str(entity.get("softwareVersion", ""))
            if ver == "unknown":
                warns.append(
                    "provenance: producing tool reports softwareVersion "
                    "'unknown' (built outside a git checkout); the source "
                    "commit is NOT recorded."
                )
            elif ver.endswith("-dirty"):
                warns.append(
                    f"provenance: producing tool was built -dirty ({ver}); "
                    "the working tree had uncommitted changes, so the "
                    "recorded commit does not fully describe the binary."
                )
    return warns


def verify(asset_dir: pathlib.Path) -> int:
    crate_path = asset_dir / "ro-crate-metadata.json"
    if not crate_path.exists():
        sys.stderr.write(
            f"no ro-crate-metadata.json in {asset_dir}; "
            "pass the asset directory as the first argument\n"
        )
        return 2
    meta = json.loads(crate_path.read_text())
    recorded = _recorded_rel_paths(meta)

    errs = 0
    checked = 0

    # (1) recorded -> disk: hashes match.
    for rel, want in sorted(recorded.items()):
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

    # (2) disk -> recorded: no unaccounted files. The manifest file
    # itself has no self-hash, and __pycache__ is a transient byproduct
    # of running these very tools -- both are exempt.
    recorded_set = set(recorded)
    for path in sorted(asset_dir.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(asset_dir).as_posix()
        if rel == "ro-crate-metadata.json":
            continue
        if "__pycache__" in path.parts:
            continue
        if rel not in recorded_set:
            print(f"{rel}: UNEXPECTED (present on disk but not in the manifest)")
            errs += 1

    # (3) provenance sanity (warnings, not failures).
    for w in _provenance_warnings(meta):
        print(f"WARNING: {w}")

    if errs:
        print(f"{errs} problem(s); {checked} file(s) verified")
        return 1
    print(f"OK ({checked} file(s) verified; no unexpected files)")
    return 0


def main(argv: list[str]) -> int:
    asset_dir = pathlib.Path(argv[1] if len(argv) >= 2 else ".")
    return verify(asset_dir)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
