//! High-level drawing helpers built on `Scene` primitives.
//!
//! [`MarkerStyle`] and [`TileStyle`] bundle the per-item style
//! choices needed for the common shapes (decorated polygons,
//! point clouds with labels), and the extension methods on
//! [`Scene`] turn one such bundle into the right sequence of
//! [`Item`]s.

use super::plotutils::{P64, tile_centroid};
use super::scene::{ArrowStyle, Color, Fill, Item, MarkerShape, Scene, Stroke, TextStyle};

// ---- Marker style ----

/// Visual style for a decorative point marker (vertex dot, etc.).
#[derive(Clone, Debug)]
pub struct MarkerStyle {
    pub shape: MarkerShape,
    /// Diameter (circle) or side length (square), in scene units.
    pub size: f64,
    pub fill: Option<Fill>,
    pub stroke: Option<Stroke>,
}

impl MarkerStyle {
    /// Solid filled circle.
    pub fn filled_circle(size: f64, color: Color) -> Self {
        MarkerStyle {
            shape: MarkerShape::Circle,
            size,
            fill: Some(Fill::solid(color)),
            stroke: None,
        }
    }

    /// Outlined circle (no fill).
    pub fn outlined_circle(size: f64, stroke: Stroke) -> Self {
        MarkerStyle {
            shape: MarkerShape::Circle,
            size,
            fill: None,
            stroke: Some(stroke),
        }
    }

    /// Solid filled square.
    pub fn filled_square(size: f64, color: Color) -> Self {
        MarkerStyle {
            shape: MarkerShape::Square,
            size,
            fill: Some(Fill::solid(color)),
            stroke: None,
        }
    }
}

// ---- Tile style ----

/// Visual style for a closed polygonal tile boundary plus optional
/// per-vertex decorations and a centroid label.
///
/// All fields are independent: any combination is valid. The empty
/// [`TileStyle::default`] draws nothing — use one of the
/// constructors ([`TileStyle::outlined`], [`TileStyle::filled`]) or
/// the builder methods to opt in to each visual element.
#[derive(Clone, Debug, Default)]
pub struct TileStyle {
    pub fill: Option<Fill>,
    pub border: Option<Stroke>,
    /// If set, the border draws a filled triangular arrowhead at
    /// the end of each edge (in vertex order). Requires `border` to
    /// also be set, since the arrowhead colour follows the border
    /// stroke's colour.
    pub edge_arrows: Option<ArrowStyle>,
    pub vertex_marker: Option<MarkerStyle>,
    /// If true, each vertex gets its 0-based index as a text label.
    pub vertex_labels: bool,
    /// Style for vertex labels. Required when `vertex_labels` is true.
    pub vertex_label_style: Option<TextStyle>,
    /// Optional text drawn at the tile's centroid.
    pub center_label: Option<(String, TextStyle)>,
}

impl TileStyle {
    pub fn outlined(border: Stroke) -> Self {
        TileStyle {
            border: Some(border),
            ..Self::default()
        }
    }

    pub fn filled(fill: Fill, border: Stroke) -> Self {
        TileStyle {
            fill: Some(fill),
            border: Some(border),
            ..Self::default()
        }
    }

    pub fn with_vertex_marker(mut self, marker: MarkerStyle) -> Self {
        self.vertex_marker = Some(marker);
        self
    }

    pub fn with_vertex_labels(mut self, style: TextStyle) -> Self {
        self.vertex_labels = true;
        self.vertex_label_style = Some(style);
        self
    }

    pub fn with_center_label(mut self, text: impl Into<String>, style: TextStyle) -> Self {
        self.center_label = Some((text.into(), style));
        self
    }

    /// Draw filled-triangle arrowheads at the **end** of each edge
    /// (uses [`ArrowStyle::per_edge`]). The tip aligns with the
    /// edge's terminal vertex. Only has an effect when `border` is
    /// also set, since the arrowhead colour follows it.
    pub fn with_edge_arrows(mut self, size: f64) -> Self {
        self.edge_arrows = Some(ArrowStyle::per_edge(size));
        self
    }

    /// Draw filled-triangle arrowheads at the **midpoint** of each
    /// edge (uses [`ArrowStyle::per_edge_mid`]). Useful as a direction
    /// indicator when vertices already carry markers / labels.
    pub fn with_edge_arrows_mid(mut self, size: f64) -> Self {
        self.edge_arrows = Some(ArrowStyle::per_edge_mid(size));
        self
    }
}

// ---- Scene extension methods ----

impl Scene {
    /// Draw an open polyline with the given stroke.
    pub fn draw_polyline(&mut self, points: &[P64], stroke: Stroke) -> &mut Self {
        self.push(Item::Polyline {
            points: points.to_vec(),
            stroke,
            arrow: None,
        })
    }

    /// Like [`Scene::draw_polyline`] but decorated with an arrowhead
    /// (or arrowheads, per [`ArrowStyle::place`]).
    pub fn draw_polyline_with_arrow(
        &mut self,
        points: &[P64],
        stroke: Stroke,
        arrow: ArrowStyle,
    ) -> &mut Self {
        self.push(Item::Polyline {
            points: points.to_vec(),
            stroke,
            arrow: Some(arrow),
        })
    }

    /// Draw a single line segment.
    pub fn draw_segment(&mut self, a: P64, b: P64, stroke: Stroke) -> &mut Self {
        self.push(Item::Segment {
            a,
            b,
            stroke,
            arrow: None,
        })
    }

    /// Draw a directed arrow from `a` to `b`: a stroked segment
    /// with a filled triangular arrowhead at `b`. Equivalent to a
    /// `Segment` with [`ArrowStyle::end`].
    pub fn draw_arrow(&mut self, a: P64, b: P64, stroke: Stroke, head_size: f64) -> &mut Self {
        self.push(Item::Segment {
            a,
            b,
            stroke,
            arrow: Some(ArrowStyle::end(head_size)),
        })
    }

    /// Draw a closed polygonal tile boundary with the given style.
    ///
    /// If `vertices` ends with a duplicate of its first entry
    /// (the format produced by `Snake::to_polyline_f64` /
    /// `Rat::to_polyline_f64`), the trailing duplicate is trimmed —
    /// SVG `<polygon>` closes implicitly.
    pub fn draw_tile(&mut self, vertices: &[P64], style: &TileStyle) -> &mut Self {
        if vertices.is_empty() {
            return self;
        }
        let mut pts: Vec<P64> = vertices.to_vec();
        if pts.len() > 1 && pts.first() == pts.last() {
            pts.pop();
        }

        // 1. Filled / stroked polygon (with optional edge arrows).
        if style.fill.is_some() || style.border.is_some() {
            self.push(Item::Polygon {
                points: pts.clone(),
                fill: style.fill,
                stroke: style.border.clone(),
                arrow: style.edge_arrows.clone(),
            });
        }

        // 2. Per-vertex markers.
        if let Some(marker) = &style.vertex_marker {
            for &p in &pts {
                self.push(Item::Marker {
                    center: p,
                    shape: marker.shape,
                    size: marker.size,
                    fill: marker.fill,
                    stroke: marker.stroke.clone(),
                });
            }
        }

        // 3. Per-vertex index labels.
        if style.vertex_labels {
            let lbl_style = style.vertex_label_style.clone().unwrap_or_else(|| {
                let size = style
                    .vertex_marker
                    .as_ref()
                    .map(|m| m.size * 0.7)
                    .unwrap_or(0.5);
                TextStyle::new(size, Color::BLACK)
            });
            for (i, &p) in pts.iter().enumerate() {
                self.push(Item::Text {
                    pos: p,
                    text: i.to_string(),
                    style: lbl_style.clone(),
                });
            }
        }

        // 4. Tile-center label.
        if let Some((text, style)) = &style.center_label {
            let center = tile_centroid(pts.iter());
            self.push(Item::Text {
                pos: center,
                text: text.clone(),
                style: style.clone(),
            });
        }

        self
    }

    /// Draw plain vertex markers, no labels.
    pub fn draw_points(&mut self, pts: &[P64], marker: &MarkerStyle) -> &mut Self {
        for &p in pts {
            self.push(Item::Marker {
                center: p,
                shape: marker.shape,
                size: marker.size,
                fill: marker.fill,
                stroke: marker.stroke.clone(),
            });
        }
        self
    }

    /// Draw markers with per-point text labels. `label_fn(i, p)`
    /// receives the 0-based index and the math-coord point; empty
    /// labels are skipped.
    pub fn draw_labeled_points<F>(
        &mut self,
        pts: &[P64],
        marker: &MarkerStyle,
        label_style: &TextStyle,
        label_fn: F,
    ) -> &mut Self
    where
        F: Fn(usize, P64) -> String,
    {
        for (i, &p) in pts.iter().enumerate() {
            self.push(Item::Marker {
                center: p,
                shape: marker.shape,
                size: marker.size,
                fill: marker.fill,
                stroke: marker.stroke.clone(),
            });
            let s = label_fn(i, p);
            if !s.is_empty() {
                self.push(Item::Text {
                    pos: p,
                    text: s,
                    style: label_style.clone(),
                });
            }
        }
        self
    }
}

// ---- Palette ----

/// HSL → RGB. `h` is wrapped to `[0, 1)`; `s`, `l` are clamped to
/// `[0, 1]`.
fn hsl_to_rgb(h: f32, s: f32, l: f32) -> Color {
    let h = h.rem_euclid(1.0);
    let s = s.clamp(0.0, 1.0);
    let l = l.clamp(0.0, 1.0);
    if s == 0.0 {
        let v = (l * 255.0).round().clamp(0.0, 255.0) as u8;
        return Color::rgb(v, v, v);
    }
    let q = if l < 0.5 {
        l * (1.0 + s)
    } else {
        l + s - l * s
    };
    let p = 2.0 * l - q;
    let r = hue_chunk(p, q, h + 1.0 / 3.0);
    let g = hue_chunk(p, q, h);
    let b = hue_chunk(p, q, h - 1.0 / 3.0);
    Color::rgb(
        (r * 255.0).round().clamp(0.0, 255.0) as u8,
        (g * 255.0).round().clamp(0.0, 255.0) as u8,
        (b * 255.0).round().clamp(0.0, 255.0) as u8,
    )
}

fn hue_chunk(p: f32, q: f32, mut t: f32) -> f32 {
    if t < 0.0 {
        t += 1.0;
    }
    if t > 1.0 {
        t -= 1.0;
    }
    if t < 1.0 / 6.0 {
        return p + (q - p) * 6.0 * t;
    }
    if t < 1.0 / 2.0 {
        return q;
    }
    if t < 2.0 / 3.0 {
        return p + (q - p) * (2.0 / 3.0 - t) * 6.0;
    }
    p
}

/// Evenly-spaced rainbow palette of `n` colors. `saturation` and
/// `lightness` are in `[0, 1]`. Returns an empty vec for `n = 0`.
pub fn rainbow(n: usize, saturation: f32, lightness: f32) -> Vec<Color> {
    (0..n)
        .map(|i| hsl_to_rgb(i as f32 / n as f32, saturation, lightness))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vis::scene::Viewport;

    #[test]
    fn rainbow_returns_n_colors() {
        let p = rainbow(7, 1.0, 0.5);
        assert_eq!(p.len(), 7);
        // Should not all be the same.
        let first = p[0];
        assert!(p.iter().any(|c| *c != first));
    }

    #[test]
    fn rainbow_zero_is_empty() {
        assert!(rainbow(0, 1.0, 0.5).is_empty());
    }

    #[test]
    fn hsl_pure_colors() {
        // h=0 sat=1 lit=0.5 → pure red.
        assert_eq!(hsl_to_rgb(0.0, 1.0, 0.5), Color::rgb(255, 0, 0));
        // h=1/3 → pure green.
        assert_eq!(hsl_to_rgb(1.0 / 3.0, 1.0, 0.5), Color::rgb(0, 255, 0));
        // h=2/3 → pure blue.
        assert_eq!(hsl_to_rgb(2.0 / 3.0, 1.0, 0.5), Color::rgb(0, 0, 255));
    }

    #[test]
    fn hsl_sat_zero_is_gray() {
        let g = hsl_to_rgb(0.5, 0.0, 0.5);
        assert_eq!(g.r, g.g);
        assert_eq!(g.g, g.b);
    }

    #[test]
    fn draw_tile_trims_duplicate_close_point() {
        let mut scene = Scene::new();
        let pts: Vec<P64> = vec![(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (0.0, 0.0)];
        scene.draw_tile(
            &pts,
            &TileStyle::outlined(Stroke::solid(Color::BLACK, 0.05)),
        );
        // Only one polygon item, with 3 (not 4) points.
        let Item::Polygon { points, .. } = &scene.items[0] else {
            panic!("expected polygon, got {:?}", scene.items[0]);
        };
        assert_eq!(points.len(), 3);
    }

    #[test]
    fn draw_tile_emits_expected_item_counts() {
        let pts: Vec<P64> = vec![(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)];
        let style = TileStyle::filled(
            Fill::solid(Color::YELLOW.with_alpha(64)),
            Stroke::solid(Color::BLACK, 0.05),
        )
        .with_vertex_marker(MarkerStyle::filled_circle(0.1, Color::RED))
        .with_vertex_labels(TextStyle::new(0.1, Color::BLACK))
        .with_center_label("foo", TextStyle::new(0.2, Color::BLACK));
        let mut scene = Scene::new();
        scene.draw_tile(&pts, &style);
        // Expected items: 1 polygon + 3 markers + 3 vertex labels + 1 center label = 8.
        assert_eq!(scene.items.len(), 8);
        assert!(matches!(scene.items[0], Item::Polygon { .. }));
        for i in 1..=3 {
            assert!(matches!(scene.items[i], Item::Marker { .. }));
        }
        for i in 4..=6 {
            assert!(matches!(scene.items[i], Item::Text { .. }));
        }
        assert!(matches!(scene.items[7], Item::Text { .. }));
    }

    #[test]
    fn draw_labeled_points_skips_empty_labels() {
        let mut scene = Scene::new();
        let pts: Vec<P64> = vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)];
        let marker = MarkerStyle::filled_circle(0.1, Color::RED);
        let lbl_style = TextStyle::new(0.1, Color::BLACK);
        // Only label the middle point.
        scene.draw_labeled_points(&pts, &marker, &lbl_style, |i, _| {
            if i == 1 {
                "mid".to_string()
            } else {
                String::new()
            }
        });
        // 3 markers + 1 label = 4 items. Interleaving: marker(0),
        // marker(1), label("mid"), marker(2) — the label for point
        // 1 is pushed immediately after its marker.
        assert_eq!(scene.items.len(), 4);
        let Item::Text { text, .. } = &scene.items[2] else {
            panic!("expected text at index 2, got {:?}", scene.items[2]);
        };
        assert_eq!(text, "mid");
    }

    #[test]
    fn integration_roundtrip_to_svg_smoke() {
        // Build a Scene exercising every Item kind and confirm
        // to_svg produces a string starting with <svg ...> and
        // containing one tag per item kind.
        let mut scene = Scene::new().with_background(Color::WHITE);
        scene.draw_segment((0.0, 0.0), (1.0, 1.0), Stroke::solid(Color::BLACK, 0.02));
        scene.draw_polyline(
            &[(0.0, 0.5), (0.3, 0.7), (0.6, 0.2)],
            Stroke::dashed(Color::BLUE, 0.02, vec![0.05, 0.05]),
        );
        scene.draw_tile(
            &[(0.5, 0.5), (1.0, 0.5), (0.75, 1.0)],
            &TileStyle::filled(
                Fill::solid(Color::YELLOW.with_alpha(80)),
                Stroke::solid(Color::BLACK, 0.02),
            )
            .with_vertex_marker(MarkerStyle::filled_square(0.05, Color::RED))
            .with_vertex_labels(TextStyle::new(0.06, Color::BLACK))
            .with_center_label("T", TextStyle::new(0.1, Color::BLACK).bold()),
        );

        let bounds = scene.auto_bounds().unwrap();
        let svg = scene.to_svg(&Viewport::square_for(200, bounds, 10));
        assert!(svg.starts_with("<svg "));
        assert!(svg.contains("<line "));
        assert!(svg.contains("<polyline "));
        assert!(svg.contains("<polygon "));
        assert!(svg.contains("<rect ")); // background + square markers
        assert!(svg.contains("<text "));
        assert!(svg.contains("stroke-dasharray="));
    }
}
