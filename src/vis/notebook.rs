//! Inline display of [`Scene`]s in evcxr-powered Jupyter notebooks.
//!
//! The evcxr Rust kernel looks for an `evcxr_display(&self)` method
//! on whatever value the last expression of a cell evaluates to. If
//! the method exists, it's called and its `println!`-emitted
//! `EVCXR_BEGIN_CONTENT … EVCXR_END_CONTENT` blocks are captured and
//! rendered inline.
//!
//! This module supplies two thin wrappers around a [`Scene`] +
//! [`Viewport`] pair:
//!
//! - [`SceneSvgDisplay`] — always available, emits the scene as
//!   inline SVG (`image/svg+xml`). The browser renders it; scales
//!   crisply at any zoom.
//! - [`ScenePngDisplay`] — feature-gated behind `raster`, emits a
//!   base64-encoded PNG. Use this when the SVG would have so many
//!   primitives that browser rendering becomes sluggish.
//!
//! # Usage in a notebook
//!
//! ```ignore
//! :dep tilezz = { path = "../tilezz" }
//! use tilezz::vis::scene::*;
//! use tilezz::vis::draw::*;
//!
//! let mut scene = Scene::new().with_background(Color::WHITE);
//! scene.draw_tile(&pts, &TileStyle::outlined(Stroke::solid(Color::BLACK, 0.05)));
//! let vp = Viewport::square_for(400, scene.auto_bounds().unwrap(), 16);
//! scene.display(&vp)
//! ```
//!
//! The trailing `scene.display(&vp)` evaluates to a
//! [`SceneSvgDisplay`], on which evcxr calls `evcxr_display`.

use super::scene::{Scene, Viewport};

const EVCXR_BEGIN: &str = "EVCXR_BEGIN_CONTENT";
const EVCXR_END: &str = "EVCXR_END_CONTENT";

/// Borrowed `Scene` + `Viewport` pair that renders inline as SVG in
/// an evcxr notebook cell. Returned by [`Scene::display`].
pub struct SceneSvgDisplay<'a> {
    scene: &'a Scene,
    vp: &'a Viewport,
}

impl<'a> SceneSvgDisplay<'a> {
    /// Construct directly. Most callers use [`Scene::display`] instead.
    pub fn new(scene: &'a Scene, vp: &'a Viewport) -> Self {
        SceneSvgDisplay { scene, vp }
    }

    /// Print the EVCXR MIME envelope around the scene's SVG. Called
    /// automatically by the evcxr kernel when a cell ends with a
    /// `SceneSvgDisplay`.
    pub fn evcxr_display(&self) {
        let svg = self.scene.to_svg(self.vp);
        println!("{EVCXR_BEGIN} image/svg+xml\n{svg}\n{EVCXR_END}");
    }
}

impl Scene {
    /// Wrap this scene + viewport in a value that displays inline as
    /// SVG in an evcxr Jupyter notebook.
    pub fn display<'a>(&'a self, vp: &'a Viewport) -> SceneSvgDisplay<'a> {
        SceneSvgDisplay::new(self, vp)
    }
}

// ---- PNG display (raster) ----

#[cfg(feature = "raster")]
/// Borrowed `Scene` + `Viewport` pair that renders inline as a
/// base64-encoded PNG in an evcxr notebook cell. Returned by
/// [`Scene::display_png`].
pub struct ScenePngDisplay<'a> {
    scene: &'a Scene,
    vp: &'a Viewport,
}

#[cfg(feature = "raster")]
impl<'a> ScenePngDisplay<'a> {
    pub fn new(scene: &'a Scene, vp: &'a Viewport) -> Self {
        ScenePngDisplay { scene, vp }
    }

    /// Print the EVCXR MIME envelope around a base64-encoded PNG.
    /// Panics if rasterization fails — for notebook ergonomics; if
    /// you need an `Err` path use [`Scene::to_png`] directly.
    pub fn evcxr_display(&self) {
        let png = self.scene.to_png(self.vp).expect("PNG render");
        let b64 = encode_base64(&png);
        println!("{EVCXR_BEGIN} image/png\n{b64}\n{EVCXR_END}");
    }
}

#[cfg(feature = "raster")]
impl Scene {
    /// Wrap this scene + viewport in a value that displays inline as
    /// a base64-encoded PNG in an evcxr Jupyter notebook. Useful
    /// when the SVG variant has too many primitives for the browser
    /// to render smoothly.
    pub fn display_png<'a>(&'a self, vp: &'a Viewport) -> ScenePngDisplay<'a> {
        ScenePngDisplay::new(self, vp)
    }
}

/// Plain-old base64 encoder (RFC 4648 standard alphabet, with `=`
/// padding). Hand-rolled to avoid pulling in a dep just for the
/// notebook PNG-display payload.
#[cfg(feature = "raster")]
fn encode_base64(bytes: &[u8]) -> String {
    const ALPHA: &[u8; 64] =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let b0 = chunk[0] as u32;
        let b1 = chunk.get(1).copied().unwrap_or(0) as u32;
        let b2 = chunk.get(2).copied().unwrap_or(0) as u32;
        let v = (b0 << 16) | (b1 << 8) | b2;
        out.push(ALPHA[((v >> 18) & 0x3F) as usize] as char);
        out.push(ALPHA[((v >> 12) & 0x3F) as usize] as char);
        if chunk.len() > 1 {
            out.push(ALPHA[((v >> 6) & 0x3F) as usize] as char);
        } else {
            out.push('=');
        }
        if chunk.len() > 2 {
            out.push(ALPHA[(v & 0x3F) as usize] as char);
        } else {
            out.push('=');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vis::scene::{Color, Item, MarkerShape, Fill};

    fn simple_scene() -> Scene {
        let mut s = Scene::new().with_background(Color::WHITE);
        s.push(Item::Marker {
            center: (0.5, 0.5),
            shape: MarkerShape::Circle,
            size: 0.2,
            fill: Some(Fill::solid(Color::RED)),
            stroke: None,
        });
        s
    }

    /// `evcxr_display` writes only to stdout, so we can't capture
    /// it from a unit test without going through a child process.
    /// Instead, verify the wrapper builds and the SVG it would emit
    /// is well-formed.
    #[test]
    fn svg_display_wrapper_carries_scene_and_viewport() {
        let scene = simple_scene();
        let vp = Viewport::square_for(64, ((0.0, 0.0), (1.0, 1.0)), 4);
        let display = scene.display(&vp);
        let svg = display.scene.to_svg(display.vp);
        assert!(svg.starts_with("<svg "));
        assert!(svg.contains("<circle "));
    }

    #[cfg(feature = "raster")]
    #[test]
    fn png_display_wrapper_renders_png_under_the_hood() {
        let scene = simple_scene();
        let vp = Viewport::square_for(64, ((0.0, 0.0), (1.0, 1.0)), 4);
        let display = scene.display_png(&vp);
        let png = display.scene.to_png(display.vp).unwrap();
        // PNG magic.
        assert_eq!(&png[..8], &[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]);
    }

    #[cfg(feature = "raster")]
    #[test]
    fn base64_known_values() {
        assert_eq!(encode_base64(b""), "");
        assert_eq!(encode_base64(b"f"), "Zg==");
        assert_eq!(encode_base64(b"fo"), "Zm8=");
        assert_eq!(encode_base64(b"foo"), "Zm9v");
        assert_eq!(encode_base64(b"foob"), "Zm9vYg==");
        assert_eq!(encode_base64(b"fooba"), "Zm9vYmE=");
        assert_eq!(encode_base64(b"foobar"), "Zm9vYmFy");
    }
}
