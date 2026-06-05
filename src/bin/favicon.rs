//! Generate the Rat Explorer favicon: the spectre tile rendered through
//! the project's own `vis` pipeline, so the icon is an authentic ZZ12 rat
//! rather than a hand-drawn approximation. Emits `favicon.svg` plus 32px
//! and 180px PNG fallbacks into `web/ratdb/`.
//!
//! Run from the repo root:
//! `cargo run --bin favicon --features dev-tools,raster`
//! (`raster` pulls in resvg for the PNG fallbacks; `dev-tools` is the
//! marker that keeps this bin out of normal `cargo install` builds).

use std::path::Path;

use tilezz::cyclotomic::ZZ12;
use tilezz::geom::snake::{Snake, Turtle};
use tilezz::geom::tiles::spectre;
use tilezz::vis::draw::TileStyle;
use tilezz::vis::plotutils::P64;
use tilezz::vis::raster;
use tilezz::vis::scene::{Color, Fill, Scene, Stroke, Viewport};

fn main() {
    // The spectre as a closed polyline in math space (y-up).
    let snake: Snake<ZZ12> = spectre();
    let polyline: Vec<P64> = snake.to_polyline_f64(Turtle::default());

    // Tight square bbox centered on the shape's own center -- the spectre
    // sits off the origin, so centering on the origin would push it into a
    // corner. The 1.12 factor bakes in a margin so the outline never kisses
    // the canvas edge.
    let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
    let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for &(x, y) in &polyline {
        min_x = min_x.min(x);
        max_x = max_x.max(x);
        min_y = min_y.min(y);
        max_y = max_y.max(y);
    }
    let cx = (min_x + max_x) / 2.0;
    let cy = (min_y + max_y) / 2.0;
    let half = ((max_x - min_x).max(max_y - min_y) / 2.0) * 1.12;
    let bounds = ((cx - half, cy - half), (cx + half, cy + half));

    // Match the explorer's rendered tile: YELLOW@alpha96 fill composited on
    // a WHITE background (= the exact on-screen pale yellow) + a BLACK
    // outline. Two deliberate departures for an icon that has to survive a
    // 16px tab downscale: drop the explorer's blue vertex dots / mid-edge
    // arrows (they collapse into noise and fight the silhouette at small
    // sizes), and thicken the relative stroke (the explorer uses
    // 0.018 * half) so the outline stays visible.
    let stroke_w = 0.05 * half;
    let build_scene = || {
        let style = TileStyle::filled(
            Fill::solid(Color::YELLOW.with_alpha(96)),
            Stroke::solid(Color::BLACK, stroke_w),
        );
        let mut scene = Scene::new().with_background(Color::WHITE);
        scene.draw_tile(&polyline, &style);
        scene
    };

    let out = Path::new("web/ratdb");

    // SVG primary: `side_px` only sets the width/height attributes; the
    // vector itself scales to whatever the browser requests.
    let svg_vp = Viewport::square_for(64, bounds, 2);
    let svg = build_scene().to_svg(&svg_vp);
    std::fs::write(out.join("favicon.svg"), &svg).expect("write favicon.svg");

    // PNG fallbacks: 32px for standard/hi-DPI browser tabs, 180px for the
    // iOS home-screen apple-touch-icon. Padding scales with the side.
    for (side, pad, name) in [
        (32u32, 1u32, "favicon-32.png"),
        (180, 6, "apple-touch-icon.png"),
    ] {
        let vp = Viewport::square_for(side, bounds, pad);
        let png = raster::to_png(&build_scene(), &vp).expect("rasterize png");
        std::fs::write(out.join(name), &png).expect("write png");
    }

    println!("wrote favicon.svg, favicon-32.png, apple-touch-icon.png to web/ratdb/");
}
