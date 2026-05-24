//! PNG rasterization for `Scene` via `resvg` + `tiny-skia`.
//!
//! This module is feature-gated behind `raster` and pulls in
//! `resvg`/`usvg`/`tiny-skia`. All three are pure-Rust crates with
//! no native dependencies, so the raster path works on every target
//! `tiny-skia` supports (including WASM).
//!
//! Text rendering depends on at least one font being available to
//! `usvg`. By default we load system fonts; on platforms without a
//! system font store (e.g. WASM) text will render as empty boxes
//! unless callers preload fonts themselves and use the lower-level
//! [`render_with_options`] helper.

use super::scene::{Scene, Viewport};

#[derive(Debug)]
pub enum RasterError {
    /// `usvg` failed to parse the SVG string emitted by [`Scene::to_svg`].
    SvgParse(String),
    /// `tiny_skia::Pixmap::new` returned `None` (zero / overflowed
    /// dimensions).
    PixmapAlloc { width: u32, height: u32 },
    /// `tiny_skia::Pixmap::encode_png` failed.
    PngEncode(String),
}

impl std::fmt::Display for RasterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RasterError::SvgParse(msg) => write!(f, "SVG parse error: {msg}"),
            RasterError::PixmapAlloc { width, height } => {
                write!(f, "pixmap allocation failed for {width}×{height} pixels")
            }
            RasterError::PngEncode(msg) => write!(f, "PNG encode error: {msg}"),
        }
    }
}

impl std::error::Error for RasterError {}

/// Default `usvg::Options` for the tilezz scene pipeline.
///
/// Loads system fonts into fontdb and, importantly, maps CSS
/// generic family names (`sans-serif`, `serif`, `monospace`) to the
/// first matching real font face found among a portable set of
/// well-known fallbacks (DejaVu, Liberation, Noto, Arial,
/// Helvetica, …). Without this step, an SVG that says
/// `font-family="sans-serif"` would not render any text because
/// `sans-serif` is a CSS generic, not a real font family in fontdb.
fn default_usvg_options() -> resvg::usvg::Options<'static> {
    let mut opt = resvg::usvg::Options::default();
    {
        let db = opt.fontdb_mut();
        db.load_system_fonts();
        configure_generic_families(db);
    }
    opt
}

fn configure_generic_families(db: &mut resvg::usvg::fontdb::Database) {
    use resvg::usvg::fontdb::Family;

    // Find an available family whose name matches any of `candidates`.
    let pick = |db: &resvg::usvg::fontdb::Database, candidates: &[&str]| -> Option<String> {
        for c in candidates {
            if db.faces().any(|f| {
                f.families
                    .iter()
                    .any(|(fam, _)| fam.eq_ignore_ascii_case(c))
            }) {
                return Some((*c).to_string());
            }
        }
        // Last-ditch fallback: any face we have.
        db.faces()
            .next()
            .and_then(|f| f.families.first().map(|(fam, _)| fam.clone()))
    };

    if let Some(name) = pick(
        db,
        &[
            "DejaVu Sans",
            "Liberation Sans",
            "Noto Sans",
            "Arial",
            "Helvetica",
            "Segoe UI",
        ],
    ) {
        db.set_sans_serif_family(name);
    }
    if let Some(name) = pick(
        db,
        &[
            "DejaVu Serif",
            "Liberation Serif",
            "Noto Serif",
            "Times New Roman",
            "Times",
        ],
    ) {
        db.set_serif_family(name);
    }
    if let Some(name) = pick(
        db,
        &[
            "DejaVu Sans Mono",
            "Liberation Mono",
            "Noto Sans Mono",
            "Consolas",
            "Courier New",
            "Courier",
        ],
    ) {
        db.set_monospace_family(name);
    }

    // Silence "Family" import as it's used implicitly via fontdb internals
    // on some versions.
    let _ = std::marker::PhantomData::<Family>;
}

/// Render a scene to PNG bytes using the default options
/// (system-font fontdb). Convenience wrapper around
/// [`render_with_options`].
pub fn to_png(scene: &Scene, vp: &Viewport) -> Result<Vec<u8>, RasterError> {
    render_with_options(scene, vp, &default_usvg_options())
}

/// Render a scene to PNG bytes using caller-supplied `usvg::Options`.
/// Use this when you need to control font loading (e.g. embed a font
/// directly for WASM targets).
pub fn render_with_options(
    scene: &Scene,
    vp: &Viewport,
    opt: &resvg::usvg::Options<'_>,
) -> Result<Vec<u8>, RasterError> {
    let svg = scene.to_svg(vp);
    let tree =
        resvg::usvg::Tree::from_str(&svg, opt).map_err(|e| RasterError::SvgParse(e.to_string()))?;

    let mut pixmap =
        tiny_skia::Pixmap::new(vp.width_px, vp.height_px).ok_or(RasterError::PixmapAlloc {
            width: vp.width_px,
            height: vp.height_px,
        })?;

    resvg::render(
        &tree,
        tiny_skia::Transform::identity(),
        &mut pixmap.as_mut(),
    );

    pixmap
        .encode_png()
        .map_err(|e| RasterError::PngEncode(e.to_string()))
}

impl Scene {
    /// Render this scene to a PNG-encoded byte buffer at `vp`'s
    /// pixel resolution. Uses system fonts for text rendering.
    pub fn to_png(&self, vp: &Viewport) -> Result<Vec<u8>, RasterError> {
        to_png(self, vp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vis::draw::{MarkerStyle, TileStyle};
    use crate::vis::scene::{Color, Fill, Stroke, TextStyle};

    /// PNG magic bytes: 89 50 4E 47 0D 0A 1A 0A.
    const PNG_MAGIC: &[u8] = &[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a];

    /// IHDR chunk layout: 8 magic + 4 length + 4 "IHDR" +
    ///   4 width (big-endian) + 4 height (big-endian) + ...
    fn parse_png_dims(bytes: &[u8]) -> (u32, u32) {
        assert!(bytes.starts_with(PNG_MAGIC), "PNG magic missing");
        // After magic: 4-byte chunk length, 4-byte "IHDR".
        assert_eq!(&bytes[12..16], b"IHDR", "IHDR chunk missing");
        let width = u32::from_be_bytes(bytes[16..20].try_into().unwrap());
        let height = u32::from_be_bytes(bytes[20..24].try_into().unwrap());
        (width, height)
    }

    #[test]
    fn empty_scene_renders_blank_png() {
        let scene = Scene::new().with_background(Color::WHITE);
        let vp = Viewport::square_for(64, ((0.0, 0.0), (1.0, 1.0)), 0);
        let png = scene.to_png(&vp).expect("PNG render");
        assert_eq!(parse_png_dims(&png), (64, 64));
    }

    #[test]
    fn rectangular_viewport_produces_matching_dims() {
        let mut scene = Scene::new().with_background(Color::WHITE);
        scene.push(crate::vis::scene::Item::Marker {
            center: (0.5, 0.5),
            shape: crate::vis::scene::MarkerShape::Circle,
            size: 0.2,
            fill: Some(Fill::solid(Color::RED)),
            stroke: None,
        });
        let vp = Viewport::rect_for(200, 100, ((0.0, 0.0), (1.0, 1.0)), 0);
        let png = scene.to_png(&vp).expect("PNG render");
        assert_eq!(parse_png_dims(&png), (200, 100));
    }

    /// Render a representative scene using real tilezz tile data
    /// (hexagon / square / spectre snake polylines) + a rainbow
    /// point row + a dashed segment + a label, and dump both SVG
    /// and PNG to the project root for eyeball inspection.
    ///
    /// Run with: `cargo test --features raster -- --ignored dump_visual_smoke`.
    #[test]
    #[ignore = "writes files to project root; run on demand"]
    fn dump_visual_smoke() {
        use crate::cyclotomic::ZZ12;
        use crate::intgeom::rat::Rat;
        use crate::intgeom::snake::Turtle;
        use crate::intgeom::tiles;
        use crate::vis::draw::rainbow;
        use crate::vis::scene::{HAlign, Item, MarkerShape, VAlign};

        // Lay three real tiles in a row by walking each with a
        // shifted Turtle.
        // Three real tilezz tiles in a row, each demonstrating a
        // different arrow placement (none / per-edge-end / per-edge-mid).
        let positions = [(0.0_f64, 0.0_f64), (3.5, 0.0), (8.0, 0.0)];
        let snakes: [Vec<(f64, f64)>; 3] = [
            Rat::from_unchecked(&tiles::hexagon::<ZZ12>()).to_polyline_f64(Turtle::default()),
            Rat::from_unchecked(&tiles::square::<ZZ12>()).to_polyline_f64(Turtle::default()),
            Rat::from_unchecked(&tiles::spectre::<ZZ12>()).to_polyline_f64(Turtle::default()),
        ];
        let names = ["hex", "square↗", "spectre▶"];
        let palette = rainbow(3, 0.7, 0.55);

        let mut scene = Scene::new().with_background(Color::WHITE);

        for (i, ((shape, (dx, dy)), (label, fill_color))) in snakes
            .iter()
            .zip(positions.iter())
            .zip(names.iter().zip(palette.iter()))
            .enumerate()
        {
            let shifted: Vec<(f64, f64)> = shape.iter().map(|(x, y)| (x + dx, y + dy)).collect();
            let mut style = TileStyle::filled(
                Fill::solid(fill_color.with_alpha(96)),
                Stroke::solid(Color::BLACK, 0.06),
            )
            .with_vertex_marker(MarkerStyle::filled_circle(0.18, Color::RED))
            .with_vertex_labels(
                TextStyle::new(0.22, Color::WHITE)
                    .bold()
                    .align(HAlign::Center, VAlign::Middle),
            )
            .with_center_label(*label, TextStyle::new(0.4, Color::BLACK).bold());
            // i==0: plain. i==1: per-edge-end arrows. i==2: per-edge-mid.
            style = match i {
                1 => style.with_edge_arrows(0.25),
                2 => style.with_edge_arrows_mid(0.3),
                _ => style,
            };
            scene.draw_tile(&shifted, &style);
        }

        // Rainbow row below the tiles.
        for (i, c) in rainbow(10, 1.0, 0.55).iter().enumerate() {
            let x = -1.0 + (i as f64) * 1.2;
            scene.push(Item::Marker {
                center: (x, -2.0),
                shape: MarkerShape::Circle,
                size: 0.35,
                fill: Some(Fill::solid(*c)),
                stroke: Some(Stroke::solid(Color::BLACK, 0.04)),
            });
        }

        // A dashed underline segment.
        scene.draw_segment(
            (-1.5, -2.6),
            (11.5, -2.6),
            Stroke::dashed(Color::BLUE, 0.05, vec![0.3, 0.15]),
        );

        // A bold caption.
        scene.push(Item::Text {
            pos: (5.0, 2.5),
            text: "tilezz scene-graph demo".to_string(),
            style: TextStyle::new(0.6, Color::BLACK)
                .bold()
                .align(HAlign::Center, VAlign::Middle),
        });

        let bounds = scene.auto_bounds().unwrap();
        let vp = Viewport::rect_for(1200, 500, bounds, 24);

        let svg_path = concat!(env!("CARGO_MANIFEST_DIR"), "/vis_demo.svg");
        let png_path = concat!(env!("CARGO_MANIFEST_DIR"), "/vis_demo.png");
        std::fs::write(svg_path, scene.to_svg(&vp)).unwrap();
        std::fs::write(png_path, scene.to_png(&vp).unwrap()).unwrap();
        eprintln!("wrote {svg_path}");
        eprintln!("wrote {png_path}");
    }

    #[test]
    fn full_scene_renders_to_png() {
        // Exercise every primitive we have, end-to-end.
        let mut scene = Scene::new().with_background(Color::WHITE);
        scene.draw_segment((0.0, 0.0), (1.0, 1.0), Stroke::solid(Color::BLACK, 0.02));
        scene.draw_polyline(
            &[(0.0, 0.5), (0.3, 0.7), (0.6, 0.2)],
            Stroke::dashed(Color::BLUE, 0.02, vec![0.05, 0.05]),
        );
        scene.draw_tile(
            &[(0.5, 0.5), (1.0, 0.5), (0.75, 1.0)],
            &TileStyle::filled(
                Fill::solid(Color::YELLOW.with_alpha(96)),
                Stroke::solid(Color::BLACK, 0.02),
            )
            .with_vertex_marker(MarkerStyle::filled_square(0.05, Color::RED))
            .with_center_label("T", TextStyle::new(0.1, Color::BLACK).bold()),
        );
        let bounds = scene.auto_bounds().unwrap();
        let vp = Viewport::square_for(256, bounds, 8);
        let png = scene.to_png(&vp).expect("PNG render");
        assert_eq!(parse_png_dims(&png), (256, 256));
        // Sanity: non-blank canvas, more than a tiny header.
        assert!(
            png.len() > 200,
            "PNG suspiciously small ({} bytes)",
            png.len()
        );
    }
}
