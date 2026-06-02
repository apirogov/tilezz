#!/usr/bin/env python3
"""Standalone decoder for a tilezz-rat-dafsa-blocks asset.

Reads a directory produced by `rat_enum --mode dafsa-blocks` and
streams every accepted sequence (one rat per line, space-separated
signed integers) to stdout.

Run from the asset directory root:

    python3 tools/decode.py > rats.txt

Or pass an explicit directory:

    python3 tools/decode.py path/to/asset > rats.txt

No external Python dependencies; relies only on the stdlib (gzip,
struct, json, pathlib). The full wire format is documented in
schemas/blocks_schema.txt and schemas/rat_schema.txt next to this
file -- the code here is a literal transcription of those specs.
"""

from __future__ import annotations

import gzip
import json
import struct
import sys
from pathlib import Path
from typing import Iterator


def decode_block(raw: bytes) -> tuple[int, list[tuple[int, int, bool]], list[tuple[int, int]]]:
    """Parse one gunzipped block file into (first_state_id, states, edges).

    states is a list of (edges_offset, count, is_accept); edges is a list
    of (label, target_state_id). Mirrors the layout in blocks_schema.txt.
    """
    if raw[:4] != b"TRB1":
        raise ValueError(f"bad block magic: {raw[:4]!r}")
    first_state, n_states, n_edges = struct.unpack_from("<III", raw, 4)
    cursor = 16  # 4 (magic) + 12 (three u32s)
    states = []
    for _ in range(n_states):
        edges_offset, count = struct.unpack_from("<IQ", raw, cursor)
        is_accept = raw[cursor + 12] != 0
        cursor += 16  # u32 + u64 + u8 + 3 bytes padding
        states.append((edges_offset, count, is_accept))
    edges = []
    for _ in range(n_edges):
        label = struct.unpack_from("<b", raw, cursor)[0]
        target = struct.unpack_from("<I", raw, cursor + 4)[0]
        cursor += 8  # i8 + 3 bytes padding + u32
        edges.append((label, target))
    return first_state, states, edges


def render_template(template: str, block_id: int) -> str:
    """Substitute {:0N} placeholders in the manifest's filename template."""
    open_brace = template.find("{")
    if open_brace < 0:
        return template
    close_brace = template.find("}", open_brace)
    if close_brace < 0:
        return template
    spec = template[open_brace + 1 : close_brace]
    width = 6
    if spec.startswith(":0"):
        try:
            width = int(spec[2:])
        except ValueError:
            pass
    return template[:open_brace] + f"{block_id:0{width}d}" + template[close_brace + 1 :]


class BlockedDafsa:
    """A minimal random-access reader for the blocked DAFSA format."""

    def __init__(self, asset_dir: Path) -> None:
        self.dir = asset_dir
        with (asset_dir / "block_index.json").open() as f:
            self.manifest = json.load(f)
        for required in ("format", "version", "block_size", "n_states", "root_state"):
            if required not in self.manifest:
                raise ValueError(f"manifest missing required field: {required}")
        if self.manifest["format"] != "tilezz-rat-dafsa-blocks":
            raise ValueError(f"unexpected format: {self.manifest['format']!r}")
        if self.manifest["version"] != 1:
            raise ValueError(f"unsupported manifest version: {self.manifest['version']}")
        if self.manifest["block_format"] != "tilezz-rat-block":
            raise ValueError(f"unsupported block format: {self.manifest['block_format']!r}")
        if self.manifest["block_version"] != 1:
            raise ValueError(f"unsupported block version: {self.manifest['block_version']}")
        self.block_size: int = self.manifest["block_size"]
        self.n_states: int = self.manifest["n_states"]
        self.root: int = self.manifest["root_state"]
        self.template: str = self.manifest["block_filename_template"]
        self._cache: dict[int, tuple[int, list, list]] = {}

    def _get_block(self, block_id: int):
        if block_id in self._cache:
            return self._cache[block_id]
        path = self.dir / render_template(self.template, block_id)
        with gzip.open(path, "rb") as f:
            raw = f.read()
        decoded = decode_block(raw)
        self._cache[block_id] = decoded
        return decoded

    def _state(self, state_id: int) -> tuple[int, int, bool, list[tuple[int, int]]]:
        """Return (edges_offset, count, is_accept, [(label, target), ...]) for state_id."""
        block_id = state_id // self.block_size
        local = state_id % self.block_size
        first_state, states, edges = self._get_block(block_id)
        edges_offset, count, is_accept = states[local]
        if local + 1 < len(states):
            next_offset = states[local + 1][0]
        else:
            next_offset = len(edges)
        return edges_offset, count, is_accept, edges[edges_offset:next_offset]

    def iter_sequences(self) -> Iterator[list[int]]:
        """Yield every accepted sequence in lex order under the DAFSA's
        own length-prefix encoding. The first byte of each emitted
        sequence is the length prefix described in rat_schema.txt;
        callers that want the raw rat can strip it (see iter_rats).
        """
        # Depth-first traversal of the DAFSA, emitting any accepted state.
        stack: list[tuple[int, list[int], int]] = [(self.root, [], 0)]
        # Each frame is (state_id, prefix_so_far, edge_cursor).
        while stack:
            state_id, prefix, cursor = stack[-1]
            edges_offset, count, is_accept, edges = self._state(state_id)
            if cursor == 0 and is_accept:
                yield list(prefix)
            if cursor < len(edges):
                stack[-1] = (state_id, prefix, cursor + 1)
                label, target = edges[cursor]
                stack.append((target, prefix + [label], 0))
            else:
                stack.pop()

    def iter_rats(self) -> Iterator[list[int]]:
        """Yield every stored rat (length-prefix already stripped)."""
        for seq in self.iter_sequences():
            # rat_schema.txt: stored_sequence(rat) = [len(rat), rat...]
            assert seq, "DAFSA emitted an empty accepted sequence"
            assert seq[0] == len(seq) - 1, f"length-prefix mismatch: {seq}"
            yield seq[1:]


def main(argv: list[str]) -> int:
    asset_dir = Path(argv[1] if len(argv) >= 2 else ".")
    if not (asset_dir / "block_index.json").exists():
        sys.stderr.write(
            f"no block_index.json in {asset_dir}; pass the asset directory as the first argument\n"
        )
        return 2
    bd = BlockedDafsa(asset_dir)
    sys.stderr.write(
        f"# tilezz-rat-dafsa-blocks: {bd.manifest.get('n_sequences', '?')} sequences "
        f"across {bd.manifest.get('n_blocks', '?')} block(s)\n"
    )
    out = sys.stdout
    for rat in bd.iter_rats():
        out.write(" ".join(str(x) for x in rat))
        out.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
