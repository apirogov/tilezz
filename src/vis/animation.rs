//! Multi-frame GIF rendering from a sequence of `Scene`s.
//!
//! Feature-gated behind `animation`, which implies `raster` (the
//! per-frame rasterizer). Pure Rust: each `Scene` is rendered via
//! `resvg` + `tiny-skia` to an RGBA pixmap, then the `gif` crate
//! quantizes the RGBA buffer to a 256-colour palette per frame and
//! emits a GIF frame.
//!
//! All frames must use the same `Viewport` (pixel dimensions are
//! baked into the GIF header). Per-frame delays are configurable.

use super::raster::{default_options, render_to_pixmap, RasterError};
use super::scene::{Scene, Viewport};

#[derive(Debug)]
pub enum AnimationError {
    /// A scene failed to rasterise.
    Raster(RasterError),
    /// The `gif` crate's encoder reported an error.
    Encode(String),
}

impl std::fmt::Display for AnimationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AnimationError::Raster(e) => write!(f, "frame rasterisation failed: {e}"),
            AnimationError::Encode(msg) => write!(f, "GIF encode failed: {msg}"),
        }
    }
}

impl std::error::Error for AnimationError {}

impl From<RasterError> for AnimationError {
    fn from(e: RasterError) -> Self {
        AnimationError::Raster(e)
    }
}

/// Render a sequence of scenes into an animated GIF. Every scene
/// shares the same viewport; the frame delay applies to all frames.
///
/// `frame_delay_ms` is the wait between frames in milliseconds; GIF
/// stores this in centiseconds, so the value is divided by 10 (and
/// clamped to at least 1).
pub fn render_gif(
    scenes: &[Scene],
    vp: &Viewport,
    frame_delay_ms: u16,
) -> Result<Vec<u8>, AnimationError> {
    let opt = default_options();
    let mut buf: Vec<u8> = Vec::new();
    let mut encoder = gif::Encoder::new(&mut buf, vp.width_px as u16, vp.height_px as u16, &[])
        .map_err(|e| AnimationError::Encode(e.to_string()))?;
    encoder
        .set_repeat(gif::Repeat::Infinite)
        .map_err(|e| AnimationError::Encode(e.to_string()))?;

    let delay_cs = (frame_delay_ms / 10).max(1);

    for scene in scenes {
        let pixmap = render_to_pixmap(scene, vp, &opt)?;
        let mut rgba = pixmap.data().to_vec();
        let mut frame =
            gif::Frame::from_rgba(vp.width_px as u16, vp.height_px as u16, rgba.as_mut_slice());
        frame.delay = delay_cs;
        encoder
            .write_frame(&frame)
            .map_err(|e| AnimationError::Encode(e.to_string()))?;
    }
    drop(encoder);
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vis::scene::{Color, Fill, Item, MarkerShape};

    /// GIF magic: "GIF89a" or "GIF87a".
    fn parse_gif_header(bytes: &[u8]) -> Option<(u16, u16)> {
        if bytes.len() < 10 {
            return None;
        }
        if !bytes.starts_with(b"GIF89a") && !bytes.starts_with(b"GIF87a") {
            return None;
        }
        // After 6-byte signature: 2-byte width, 2-byte height (LE).
        let w = u16::from_le_bytes([bytes[6], bytes[7]]);
        let h = u16::from_le_bytes([bytes[8], bytes[9]]);
        Some((w, h))
    }

    fn marker_scene(color: Color) -> Scene {
        let mut s = Scene::new().with_background(Color::WHITE);
        s.push(Item::Marker {
            center: (0.5, 0.5),
            shape: MarkerShape::Circle,
            size: 0.3,
            fill: Some(Fill::solid(color)),
            stroke: None,
        });
        s
    }

    #[test]
    fn render_gif_writes_valid_header() {
        let scenes = vec![marker_scene(Color::RED), marker_scene(Color::BLUE)];
        let vp = Viewport::square_for(64, ((0.0, 0.0), (1.0, 1.0)), 4);
        let bytes = render_gif(&scenes, &vp, 200).expect("render gif");
        let (w, h) = parse_gif_header(&bytes).expect("valid GIF header");
        assert_eq!((w, h), (64, 64));
        // A 2-frame GIF should be larger than a 1-frame GIF; here
        // we just check it's not a header-only stub.
        assert!(
            bytes.len() > 200,
            "GIF suspiciously small ({} bytes)",
            bytes.len()
        );
    }

    #[test]
    fn render_gif_zero_delay_clamped_to_one_cs() {
        // Delay 5 ms → 0 centisecond after div 10, clamped to 1.
        let scenes = vec![marker_scene(Color::RED)];
        let vp = Viewport::square_for(32, ((0.0, 0.0), (1.0, 1.0)), 0);
        let bytes = render_gif(&scenes, &vp, 5).expect("render gif");
        assert!(parse_gif_header(&bytes).is_some());
    }

    #[test]
    fn render_gif_empty_scene_list_returns_header_only() {
        let vp = Viewport::square_for(16, ((0.0, 0.0), (1.0, 1.0)), 0);
        let bytes = render_gif(&[], &vp, 100).expect("render gif");
        // Empty animation still produces a valid header.
        assert!(parse_gif_header(&bytes).is_some());
    }
}
