# Default target: build everything and start a local server.
default: dev

# Compile the Rust crate to wasm32-unknown-unknown and stage the
# JS bindings under web/ratdb/pkg/. Equivalent for the GitHub Pages
# workflow; tweak both together if the invocation changes.
wasm:
    wasm-pack build --target web --out-dir web/ratdb/pkg -- --no-default-features --features rat_explorer

# Generate the small local test datasets the web app loads at
# startup, plus the top-level RO-Crate that lists them. Mirrors the
# Pages workflow: ZZ12 free up to n=10, plus ZZ7 (the odd ring,
# computed as ZZ14 step-2; effectiveRing 7) as an end-to-end test of
# multi-dataset support and odd-ring presentation. Re-run after
# asset-format changes; otherwise once per checkout is enough.
data:
    cargo build --release --bin rat_enum --bin build_web_rocrate --features cli,rat_explorer
    rm -rf web/ratdb/data/zz12_n10_free web/ratdb/data/zz7_n10_free
    ./target/release/rat_enum --ring 12 -n 10 --free \
        --mode dafsa-blocks --threads 0 \
        --oeis-a-number A316192 \
        -o web/ratdb/data/zz12_n10_free
    ./target/release/rat_enum --ring 14 --step 2 -n 10 --free \
        --mode dafsa-blocks --threads 0 \
        -o web/ratdb/data/zz7_n10_free
    ./target/release/build_web_rocrate --web-dir web/ratdb/

# Serve the ratdb app (web/ratdb/) via Python's stdlib http server. Open the printed URL.
serve:
    @echo "Open http://localhost:8000/"
    cd web/ratdb && python3 -m http.server 8000 --bind 127.0.0.1

# Like `serve` but bound to 0.0.0.0 so phones / other devices on the
# same LAN can hit the page at http://<this-machine-IP>:8000/.
serve-lan:
    @echo "LAN: open http://$(hostname -I | awk '{print $1}'):8000/ on the other device"
    cd web/ratdb && python3 -m http.server 8000 --bind 0.0.0.0

# Build WASM + dataset + RO-Crate, then serve. One command from a
# fresh checkout to a working page in the browser.
dev: wasm data serve

# Optional: rebuild WASM on every Rust edit. Run alongside `just serve`.
# Requires `cargo install cargo-watch --locked`.
watch:
    cargo watch -s 'wasm-pack build --target web --out-dir web/ratdb/pkg -- --no-default-features --features rat_explorer'

# Remove generated WASM artifacts only.
clean:
    rm -rf web/ratdb/pkg

# Remove WASM artifacts + the generated datasets + top-level RO-Crate.
# Use after touching the asset format so stale on-disk data doesn't
# poison the next `just dev`.
clean-all: clean
    rm -rf web/ratdb/data web/ratdb/ro-crate-metadata.json
