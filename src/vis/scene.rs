//! Declarative 2D scene model + SVG emitter.
//!
//! A `Scene` is a flat list of geometry primitives (lines, polygons,
//! markers, text) in math coordinates (y-up). Rendering is a
//! separate step: [`Scene::to_svg`] emits a well-formed SVG string
//! with zero external dependencies. The raster path
//! ([`Scene::to_png`], feature `raster`) lives in `super::raster`.
//!
//! The scene model is intentionally narrow — it covers exactly the
//! primitives the tilezz binaries / notebooks actually use:
//!
//! - line segments and open polylines,
//! - filled / stroked polygons,
//! - circle / square markers (for vertex decorations),
//! - aligned text labels.
//!
//! No charts, no axes, no grids, no auto-legends.
//!
//! # Units
//!
//! Everything in a [`Scene`] is in math coordinates: stroke width,
//! marker size, and font size are in the same units as point
//! positions. The [`Viewport`] supplied to `to_svg` / `to_png`
//! converts to pixel space (with uniform scale, equal aspect, and
//! configurable padding).

use super::plotutils::{P64, R64};

// ---- Colors / styles ----

/// 8-bit RGBA color. Use [`Color::rgb`] or [`Color::rgba`] to
/// construct, or one of the convenience constants.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

impl Color {
    pub const BLACK: Color = Color::rgb(0, 0, 0);
    pub const WHITE: Color = Color::rgb(255, 255, 255);
    pub const RED: Color = Color::rgb(255, 0, 0);
    pub const GREEN: Color = Color::rgb(0, 200, 0);
    pub const BLUE: Color = Color::rgb(0, 0, 255);
    pub const YELLOW: Color = Color::rgb(255, 220, 0);
    pub const TRANSPARENT: Color = Color::rgba(0, 0, 0, 0);

    pub const fn rgb(r: u8, g: u8, b: u8) -> Self {
        Color { r, g, b, a: 255 }
    }
    pub const fn rgba(r: u8, g: u8, b: u8, a: u8) -> Self {
        Color { r, g, b, a }
    }
    pub const fn with_alpha(self, a: u8) -> Self {
        Color {
            r: self.r,
            g: self.g,
            b: self.b,
            a,
        }
    }
}

/// Stroke style: color, width (in scene units), optional dash array.
#[derive(Clone, Debug)]
pub struct Stroke {
    pub color: Color,
    pub width: f64,
    pub dash: Option<Vec<f64>>,
}

impl Stroke {
    pub fn solid(color: Color, width: f64) -> Self {
        Stroke {
            color,
            width,
            dash: None,
        }
    }

    pub fn dashed(color: Color, width: f64, dash: Vec<f64>) -> Self {
        Stroke {
            color,
            width,
            dash: Some(dash),
        }
    }
}

/// Fill style. Alpha lives on the color.
#[derive(Clone, Copy, Debug)]
pub struct Fill {
    pub color: Color,
}

impl Fill {
    pub fn solid(color: Color) -> Self {
        Fill { color }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MarkerShape {
    Circle,
    Square,
}

/// Where arrowheads sit on a line / polyline / polygon edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ArrowPlace {
    /// One arrowhead at the line's terminal endpoint (`b` for a
    /// `Segment`, last point for a `Polyline`). Not meaningful for
    /// `Polygon` (use [`Self::PerEdge`]).
    End,
    /// One arrowhead at the line's start (`a` for a `Segment`,
    /// first point for a `Polyline`).
    Start,
    /// Arrowheads at both endpoints.
    Both,
    /// One arrowhead at the **end of every edge** — for a `Polyline`
    /// of `N` points that's `N - 1` arrowheads (at `points[1..]`);
    /// for a `Polygon` of `N` vertices that's `N` arrowheads (one
    /// per edge, in cyclic order).
    PerEdge,
    /// One arrowhead at the **midpoint of every edge**, pointing
    /// along the edge's traversal direction. Triangle is centered
    /// on the midpoint rather than ending there. For a `Segment`
    /// this puts a single direction triangle at its center; for a
    /// `Polyline` / `Polygon` it places one per edge.
    PerEdgeMid,
}

/// Arrowhead decoration on a stroked shape. The triangular tip is
/// filled in the stroke's color; size is the arrowhead's length
/// along the line, in scene units.
#[derive(Clone, Debug)]
pub struct ArrowStyle {
    pub size: f64,
    pub place: ArrowPlace,
}

impl ArrowStyle {
    /// Single arrowhead at the line's end.
    pub fn end(size: f64) -> Self {
        ArrowStyle {
            size,
            place: ArrowPlace::End,
        }
    }
    /// Single arrowhead at the line's start.
    pub fn start(size: f64) -> Self {
        ArrowStyle {
            size,
            place: ArrowPlace::Start,
        }
    }
    /// Arrowheads at both endpoints.
    pub fn both(size: f64) -> Self {
        ArrowStyle {
            size,
            place: ArrowPlace::Both,
        }
    }
    /// One arrowhead at the end of each edge.
    pub fn per_edge(size: f64) -> Self {
        ArrowStyle {
            size,
            place: ArrowPlace::PerEdge,
        }
    }
    /// One arrowhead at the midpoint of each edge (a direction
    /// indicator centered on the edge).
    pub fn per_edge_mid(size: f64) -> Self {
        ArrowStyle {
            size,
            place: ArrowPlace::PerEdgeMid,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HAlign {
    Left,
    Center,
    Right,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VAlign {
    Top,
    Middle,
    Bottom,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Weight {
    Normal,
    Bold,
}

/// Text style. `size` is in scene units.
#[derive(Clone, Debug)]
pub struct TextStyle {
    pub size: f64,
    pub color: Color,
    pub family: String,
    pub weight: Weight,
    pub italic: bool,
    pub align: HAlign,
    pub baseline: VAlign,
}

impl TextStyle {
    pub fn new(size: f64, color: Color) -> Self {
        TextStyle {
            size,
            color,
            family: "sans-serif".into(),
            weight: Weight::Normal,
            italic: false,
            align: HAlign::Center,
            baseline: VAlign::Middle,
        }
    }

    pub fn bold(mut self) -> Self {
        self.weight = Weight::Bold;
        self
    }
    pub fn italic(mut self) -> Self {
        self.italic = true;
        self
    }
    pub fn align(mut self, h: HAlign, v: VAlign) -> Self {
        self.align = h;
        self.baseline = v;
        self
    }
    pub fn family(mut self, f: impl Into<String>) -> Self {
        self.family = f.into();
        self
    }
}

// ---- Items ----

/// One drawable element in a [`Scene`]. All coordinates are math-
/// space (y-up). Sizes / widths are also in math units.
#[derive(Clone, Debug)]
pub enum Item {
    /// Single line segment from `a` to `b`. Optional `arrow`
    /// decorates the endpoints with filled triangular tips.
    Segment {
        a: P64,
        b: P64,
        stroke: Stroke,
        arrow: Option<ArrowStyle>,
    },
    /// Connected open polyline. Optional `arrow` decorates the
    /// polyline's endpoints (or per-edge — see [`ArrowPlace`]).
    Polyline {
        points: Vec<P64>,
        stroke: Stroke,
        arrow: Option<ArrowStyle>,
    },
    /// Closed polygon. Either fill or stroke (or both) must be
    /// `Some` to be visible. Optional `arrow` decorates polygon
    /// edges (per-edge is the only meaningful place for a closed
    /// shape).
    Polygon {
        points: Vec<P64>,
        fill: Option<Fill>,
        stroke: Option<Stroke>,
        arrow: Option<ArrowStyle>,
    },
    /// Decorative marker at a point. `size` is the full diameter
    /// (circle) or side length (square) in scene units.
    Marker {
        center: P64,
        shape: MarkerShape,
        size: f64,
        fill: Option<Fill>,
        stroke: Option<Stroke>,
    },
    /// Aligned text label at a math-coord position. Font size is in
    /// scene units; the renderer converts to pixels via the
    /// [`Viewport`] transform.
    Text {
        pos: P64,
        text: String,
        style: TextStyle,
    },
}

// ---- Scene ----

#[derive(Clone, Debug, Default)]
pub struct Scene {
    pub items: Vec<Item>,
    pub background: Option<Color>,
}

impl Scene {
    pub fn new() -> Self {
        Scene::default()
    }

    pub fn with_background(mut self, c: Color) -> Self {
        self.background = Some(c);
        self
    }

    pub fn push(&mut self, item: Item) -> &mut Self {
        self.items.push(item);
        self
    }

    /// Math-coord bounds tight around all items. `None` for an empty
    /// scene; in that case callers should supply an explicit
    /// [`Viewport::bounds`].
    pub fn auto_bounds(&self) -> Option<R64> {
        let mut iter = self.items.iter().flat_map(item_points);
        let first = iter.next()?;
        let (mut min_x, mut min_y) = first;
        let (mut max_x, mut max_y) = first;
        for (x, y) in iter {
            min_x = min_x.min(x);
            min_y = min_y.min(y);
            max_x = max_x.max(x);
            max_y = max_y.max(y);
        }
        Some(((min_x, min_y), (max_x, max_y)))
    }
}

fn item_points(item: &Item) -> Box<dyn Iterator<Item = P64> + '_> {
    match item {
        Item::Segment { a, b, .. } => Box::new([*a, *b].into_iter()),
        Item::Polyline { points, .. } => Box::new(points.iter().copied()),
        Item::Polygon { points, .. } => Box::new(points.iter().copied()),
        Item::Marker { center, size, .. } => {
            // Expand around the marker so its full extent is in
            // bounds, not just the center.
            let r = *size * 0.5;
            let (cx, cy) = *center;
            Box::new(
                [
                    (cx - r, cy - r),
                    (cx + r, cy - r),
                    (cx - r, cy + r),
                    (cx + r, cy + r),
                ]
                .into_iter(),
            )
        }
        // Text bounds are not meaningfully computable without font
        // metrics — exclude from auto-bounds. Caller should pad
        // viewport if labels run off-edge.
        Item::Text { .. } => Box::new(std::iter::empty()),
    }
}

// ---- Viewport ----

/// Pixel canvas + math-coord bounds the scene gets fitted into.
/// The fit is always uniform (equal aspect); the math-space window
/// is centered inside the pixel canvas with `padding_px` of empty
/// space on each side.
#[derive(Clone, Copy, Debug)]
pub struct Viewport {
    pub width_px: u32,
    pub height_px: u32,
    pub bounds: R64,
    pub padding_px: u32,
}

impl Viewport {
    /// Square pixel canvas (`side_px` × `side_px`) around `bounds`.
    pub fn square_for(side_px: u32, bounds: R64, padding_px: u32) -> Self {
        Viewport {
            width_px: side_px,
            height_px: side_px,
            bounds,
            padding_px,
        }
    }

    /// Rectangular pixel canvas around `bounds`.
    pub fn rect_for(w_px: u32, h_px: u32, bounds: R64, padding_px: u32) -> Self {
        Viewport {
            width_px: w_px,
            height_px: h_px,
            bounds,
            padding_px,
        }
    }
}

// ---- math → pixel transform ----

/// Internal: converts math-space (y-up) to pixel-space (y-down).
/// Uniform scale, equal aspect, centered inside the padded canvas.
struct PixelTransform {
    scale: f64,
    offset_x: f64,
    offset_y: f64,
    height_px: f64,
    min_x: f64,
    min_y: f64,
}

impl PixelTransform {
    fn from_viewport(vp: &Viewport) -> Self {
        let ((min_x, min_y), (max_x, max_y)) = vp.bounds;
        let math_w = (max_x - min_x).max(f64::EPSILON);
        let math_h = (max_y - min_y).max(f64::EPSILON);
        let pad = vp.padding_px as f64;
        let avail_w = (vp.width_px as f64 - 2.0 * pad).max(1.0);
        let avail_h = (vp.height_px as f64 - 2.0 * pad).max(1.0);
        let scale = (avail_w / math_w).min(avail_h / math_h);
        let used_w = math_w * scale;
        let used_h = math_h * scale;
        let offset_x = (vp.width_px as f64 - used_w) * 0.5;
        let offset_y = (vp.height_px as f64 - used_h) * 0.5;
        PixelTransform {
            scale,
            offset_x,
            offset_y,
            height_px: vp.height_px as f64,
            min_x,
            min_y,
        }
    }

    fn point(&self, (x, y): P64) -> (f64, f64) {
        let px = self.offset_x + (x - self.min_x) * self.scale;
        // y-flip: math max_y maps to the top of the pixel canvas.
        let py = self.height_px - (self.offset_y + (y - self.min_y) * self.scale);
        (px, py)
    }

    fn distance(&self, d: f64) -> f64 {
        d * self.scale
    }
}

// ---- SVG emission ----

impl Scene {
    /// Render the scene as a self-contained SVG document string.
    /// No external dependencies, no font loading; the output is
    /// trivially embeddable in HTML, Jupyter, or saved as `.svg`.
    pub fn to_svg(&self, vp: &Viewport) -> String {
        let t = PixelTransform::from_viewport(vp);
        let mut out = String::new();
        out.push_str(&format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\">",
            vp.width_px, vp.height_px, vp.width_px, vp.height_px,
        ));
        if let Some(bg) = self.background {
            out.push_str(&format!(
                "<rect width=\"{}\" height=\"{}\" fill=\"{}\" />",
                vp.width_px,
                vp.height_px,
                svg_color(bg),
            ));
        }
        for item in &self.items {
            emit_item(&mut out, item, &t);
        }
        out.push_str("</svg>");
        out
    }
}

fn emit_item(out: &mut String, item: &Item, t: &PixelTransform) {
    match item {
        Item::Segment { a, b, stroke, arrow } => {
            let (x1, y1) = t.point(*a);
            let (x2, y2) = t.point(*b);
            out.push_str(&format!(
                "<line x1=\"{:.3}\" y1=\"{:.3}\" x2=\"{:.3}\" y2=\"{:.3}\" {} />",
                x1,
                y1,
                x2,
                y2,
                stroke_attrs(stroke, t),
            ));
            if let Some(arrow) = arrow {
                emit_segment_arrows(out, *a, *b, stroke, arrow, t);
            }
        }
        Item::Polyline { points, stroke, arrow } => {
            if points.is_empty() {
                return;
            }
            let pts = points_str(points, t);
            out.push_str(&format!(
                "<polyline points=\"{}\" fill=\"none\" {} />",
                pts,
                stroke_attrs(stroke, t),
            ));
            if let Some(arrow) = arrow {
                emit_polyline_arrows(out, points, stroke, arrow, t, false);
            }
        }
        Item::Polygon {
            points,
            fill,
            stroke,
            arrow,
        } => {
            if points.is_empty() {
                return;
            }
            let pts = points_str(points, t);
            out.push_str(&format!(
                "<polygon points=\"{}\" {} {} />",
                pts,
                fill_attrs(fill.as_ref()),
                stroke_attrs_opt(stroke.as_ref(), t),
            ));
            if let (Some(stroke), Some(arrow)) = (stroke.as_ref(), arrow.as_ref()) {
                emit_polyline_arrows(out, points, stroke, arrow, t, true);
            }
        }
        Item::Marker {
            center,
            shape,
            size,
            fill,
            stroke,
        } => {
            let (cx, cy) = t.point(*center);
            let r = t.distance(*size) * 0.5;
            match shape {
                MarkerShape::Circle => {
                    out.push_str(&format!(
                        "<circle cx=\"{:.3}\" cy=\"{:.3}\" r=\"{:.3}\" {} {} />",
                        cx,
                        cy,
                        r,
                        fill_attrs(fill.as_ref()),
                        stroke_attrs_opt(stroke.as_ref(), t),
                    ));
                }
                MarkerShape::Square => {
                    out.push_str(&format!(
                        "<rect x=\"{:.3}\" y=\"{:.3}\" width=\"{:.3}\" height=\"{:.3}\" {} {} />",
                        cx - r,
                        cy - r,
                        2.0 * r,
                        2.0 * r,
                        fill_attrs(fill.as_ref()),
                        stroke_attrs_opt(stroke.as_ref(), t),
                    ));
                }
            }
        }
        Item::Text { pos, text, style } => {
            let (x, y) = t.point(*pos);
            let font_px = t.distance(style.size);
            let anchor = match style.align {
                HAlign::Left => "start",
                HAlign::Center => "middle",
                HAlign::Right => "end",
            };
            // SVG `dominant-baseline` mapping. We use math y-up but
            // text rasterizes top-down in SVG; the baseline names
            // refer to where on the glyph the anchor sits.
            let baseline = match style.baseline {
                VAlign::Top => "hanging",
                VAlign::Middle => "central",
                VAlign::Bottom => "alphabetic",
            };
            let weight = match style.weight {
                Weight::Normal => "normal",
                Weight::Bold => "bold",
            };
            let italic = if style.italic { "italic" } else { "normal" };
            out.push_str(&format!(
                "<text x=\"{:.3}\" y=\"{:.3}\" font-family=\"{}\" font-size=\"{:.3}\" \
                 font-weight=\"{}\" font-style=\"{}\" text-anchor=\"{}\" \
                 dominant-baseline=\"{}\" fill=\"{}\"{}>{}</text>",
                x,
                y,
                escape_attr(&style.family),
                font_px,
                weight,
                italic,
                anchor,
                baseline,
                svg_color(style.color),
                opacity_attr("fill-opacity", style.color.a),
                escape_text(text),
            ));
        }
    }
}

/// Emit arrowhead triangle(s) at one or both ends of a single
/// segment, depending on `arrow.place`.
fn emit_segment_arrows(
    out: &mut String,
    a: P64,
    b: P64,
    stroke: &Stroke,
    arrow: &ArrowStyle,
    t: &PixelTransform,
) {
    match arrow.place {
        ArrowPlace::End | ArrowPlace::PerEdge => {
            emit_edge_arrow(out, a, b, false, stroke, arrow.size, t);
        }
        ArrowPlace::Start => {
            emit_edge_arrow(out, b, a, false, stroke, arrow.size, t);
        }
        ArrowPlace::Both => {
            emit_edge_arrow(out, a, b, false, stroke, arrow.size, t);
            emit_edge_arrow(out, b, a, false, stroke, arrow.size, t);
        }
        ArrowPlace::PerEdgeMid => {
            emit_edge_arrow(out, a, b, true, stroke, arrow.size, t);
        }
    }
}

/// Emit arrowhead triangle(s) along an open polyline (or, if
/// `closed`, a polygon's cyclic edges).
fn emit_polyline_arrows(
    out: &mut String,
    points: &[P64],
    stroke: &Stroke,
    arrow: &ArrowStyle,
    t: &PixelTransform,
    closed: bool,
) {
    let n = points.len();
    if n < 2 {
        return;
    }
    let last_edge = if closed { n } else { n - 1 };
    match arrow.place {
        ArrowPlace::Start => {
            // First edge, reversed → tip at points[0].
            emit_edge_arrow(out, points[1], points[0], false, stroke, arrow.size, t);
        }
        ArrowPlace::End => {
            // For closed shapes "end" means the wrap-edge's end =
            // points[0]. For open shapes it's the last point.
            let (from, to) = if closed {
                (points[n - 1], points[0])
            } else {
                (points[n - 2], points[n - 1])
            };
            emit_edge_arrow(out, from, to, false, stroke, arrow.size, t);
        }
        ArrowPlace::Both => {
            emit_edge_arrow(out, points[1], points[0], false, stroke, arrow.size, t);
            let (from, to) = if closed {
                (points[n - 1], points[0])
            } else {
                (points[n - 2], points[n - 1])
            };
            emit_edge_arrow(out, from, to, false, stroke, arrow.size, t);
        }
        ArrowPlace::PerEdge => {
            for i in 0..last_edge {
                emit_edge_arrow(
                    out,
                    points[i],
                    points[(i + 1) % n],
                    false,
                    stroke,
                    arrow.size,
                    t,
                );
            }
        }
        ArrowPlace::PerEdgeMid => {
            for i in 0..last_edge {
                emit_edge_arrow(
                    out,
                    points[i],
                    points[(i + 1) % n],
                    true,
                    stroke,
                    arrow.size,
                    t,
                );
            }
        }
    }
}

/// Emit a filled triangular arrowhead for the edge `from → to`.
///
/// - If `mid == false`: tip at `to`, base sits behind it on the
///   line. The arrowhead "points at" the edge's endpoint.
/// - If `mid == true`: triangle is centered on the edge's midpoint,
///   still pointing along the traversal direction. Useful as a
///   direction indicator for an edge whose endpoints are themselves
///   decorated (vertex markers etc.).
///
/// `size` is the length of the triangle from tip to base center,
/// in scene units. Base width is fixed at 0.7 × size — gives a
/// classic isoceles arrowhead shape.
fn emit_edge_arrow(
    out: &mut String,
    from: P64,
    to: P64,
    mid: bool,
    stroke: &Stroke,
    size: f64,
    t: &PixelTransform,
) {
    let (fx, fy) = from;
    let (tx, ty) = to;
    let dx = tx - fx;
    let dy = ty - fy;
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-12 {
        return;
    }
    let ux = dx / len;
    let uy = dy / len;
    let half_base = size * 0.35;

    let (tip_x, tip_y) = if mid {
        // Centered on midpoint, tip half a size forward.
        ((fx + tx) * 0.5 + ux * size * 0.5, (fy + ty) * 0.5 + uy * size * 0.5)
    } else {
        // Tip at the endpoint.
        (tx, ty)
    };

    let base_cx = tip_x - ux * size;
    let base_cy = tip_y - uy * size;
    let nx = -uy;
    let ny = ux;
    let tip = (tip_x, tip_y);
    let left = (base_cx + nx * half_base, base_cy + ny * half_base);
    let right = (base_cx - nx * half_base, base_cy - ny * half_base);
    let pts = points_str(&[tip, left, right], t);
    out.push_str(&format!(
        "<polygon points=\"{}\" fill=\"{}\"{} stroke=\"none\" />",
        pts,
        svg_color(stroke.color),
        opacity_attr("fill-opacity", stroke.color.a),
    ));
}

fn points_str(points: &[P64], t: &PixelTransform) -> String {
    points
        .iter()
        .map(|p| {
            let (x, y) = t.point(*p);
            format!("{:.3},{:.3}", x, y)
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn stroke_attrs(s: &Stroke, t: &PixelTransform) -> String {
    let mut out = format!(
        "stroke=\"{}\" stroke-width=\"{:.3}\"{} stroke-linejoin=\"round\" stroke-linecap=\"round\"",
        svg_color(s.color),
        t.distance(s.width),
        opacity_attr("stroke-opacity", s.color.a),
    );
    if let Some(dash) = &s.dash {
        let segs = dash
            .iter()
            .map(|d| format!("{:.3}", t.distance(*d)))
            .collect::<Vec<_>>()
            .join(",");
        out.push_str(&format!(" stroke-dasharray=\"{}\"", segs));
    }
    out
}

fn stroke_attrs_opt(s: Option<&Stroke>, t: &PixelTransform) -> String {
    match s {
        Some(s) => stroke_attrs(s, t),
        None => "stroke=\"none\"".to_string(),
    }
}

fn fill_attrs(f: Option<&Fill>) -> String {
    match f {
        Some(f) => format!(
            "fill=\"{}\"{}",
            svg_color(f.color),
            opacity_attr("fill-opacity", f.color.a),
        ),
        None => "fill=\"none\"".to_string(),
    }
}

fn opacity_attr(attr: &str, alpha: u8) -> String {
    if alpha == 255 {
        String::new()
    } else {
        format!(" {}=\"{:.3}\"", attr, alpha as f64 / 255.0)
    }
}

fn svg_color(c: Color) -> String {
    format!("#{:02x}{:02x}{:02x}", c.r, c.g, c.b)
}

fn escape_text(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn escape_attr(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('"', "&quot;")
        .replace('<', "&lt;")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_viewport() -> Viewport {
        Viewport::square_for(100, ((0.0, 0.0), (1.0, 1.0)), 10)
    }

    #[test]
    fn color_constants_have_full_alpha() {
        assert_eq!(Color::BLACK.a, 255);
        assert_eq!(Color::WHITE.a, 255);
    }

    #[test]
    fn color_with_alpha_overrides() {
        let c = Color::RED.with_alpha(128);
        assert_eq!((c.r, c.g, c.b, c.a), (255, 0, 0, 128));
    }

    #[test]
    fn empty_scene_has_no_bounds() {
        let scene = Scene::new();
        assert!(scene.auto_bounds().is_none());
    }

    #[test]
    fn auto_bounds_covers_all_items() {
        let mut scene = Scene::new();
        scene.push(Item::Segment {
            a: (0.0, 0.0),
            b: (3.0, 4.0),
            stroke: Stroke::solid(Color::BLACK, 0.1),
            arrow: None,
        });
        scene.push(Item::Marker {
            center: (-1.0, -2.0),
            shape: MarkerShape::Circle,
            size: 0.5,
            fill: Some(Fill::solid(Color::RED)),
            stroke: None,
        });
        let ((min_x, min_y), (max_x, max_y)) = scene.auto_bounds().unwrap();
        assert!((min_x - -1.25).abs() < 1e-9, "min_x={min_x}");
        assert!((min_y - -2.25).abs() < 1e-9, "min_y={min_y}");
        assert!((max_x - 3.0).abs() < 1e-9, "max_x={max_x}");
        assert!((max_y - 4.0).abs() < 1e-9, "max_y={max_y}");
    }

    #[test]
    fn pixel_transform_y_up_to_y_down() {
        let vp = Viewport::square_for(100, ((0.0, 0.0), (1.0, 1.0)), 0);
        let t = PixelTransform::from_viewport(&vp);
        // (0,0) math = bottom-left pixel.
        let (x, y) = t.point((0.0, 0.0));
        assert!((x - 0.0).abs() < 1e-6);
        assert!((y - 100.0).abs() < 1e-6);
        // (1,1) math = top-right pixel.
        let (x, y) = t.point((1.0, 1.0));
        assert!((x - 100.0).abs() < 1e-6);
        assert!((y - 0.0).abs() < 1e-6);
    }

    #[test]
    fn pixel_transform_padding_centers() {
        let vp = Viewport::square_for(120, ((0.0, 0.0), (1.0, 1.0)), 10);
        let t = PixelTransform::from_viewport(&vp);
        // (0,0) math = inner bottom-left = padding offset.
        let (x, y) = t.point((0.0, 0.0));
        assert!((x - 10.0).abs() < 1e-6);
        assert!((y - 110.0).abs() < 1e-6);
    }

    #[test]
    fn pixel_transform_non_square_keeps_aspect() {
        // 200×100 px viewport with 1×1 math bounds → uniform scale
        // 100 (limited by height), content centered horizontally.
        let vp = Viewport::rect_for(200, 100, ((0.0, 0.0), (1.0, 1.0)), 0);
        let t = PixelTransform::from_viewport(&vp);
        // math (0,0) = bottom-left of centered 100×100 region → x=50.
        let (x, y) = t.point((0.0, 0.0));
        assert!((x - 50.0).abs() < 1e-6, "x={x}");
        assert!((y - 100.0).abs() < 1e-6, "y={y}");
        let (x, y) = t.point((1.0, 1.0));
        assert!((x - 150.0).abs() < 1e-6, "x={x}");
        assert!((y - 0.0).abs() < 1e-6, "y={y}");
    }

    #[test]
    fn svg_well_formed_for_segment() {
        let mut scene = Scene::new();
        scene.push(Item::Segment {
            a: (0.0, 0.0),
            b: (1.0, 1.0),
            stroke: Stroke::solid(Color::BLACK, 0.05),
            arrow: None,
        });
        let svg = scene.to_svg(&unit_viewport());
        assert!(svg.starts_with("<svg "));
        assert!(svg.ends_with("</svg>"));
        assert!(svg.contains("<line "));
        assert!(svg.contains("stroke=\"#000000\""));
    }

    #[test]
    fn svg_polygon_emits_points_attribute() {
        let mut scene = Scene::new();
        scene.push(Item::Polygon {
            points: vec![(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)],
            fill: Some(Fill::solid(Color::YELLOW.with_alpha(128))),
            stroke: Some(Stroke::solid(Color::BLACK, 0.02)),
            arrow: None,
        });
        let svg = scene.to_svg(&unit_viewport());
        assert!(svg.contains("<polygon "));
        assert!(svg.contains("points=\""));
        assert!(svg.contains("fill-opacity"));
    }

    #[test]
    fn svg_marker_circle_and_square() {
        let mut scene = Scene::new();
        scene.push(Item::Marker {
            center: (0.5, 0.5),
            shape: MarkerShape::Circle,
            size: 0.1,
            fill: Some(Fill::solid(Color::RED)),
            stroke: None,
        });
        scene.push(Item::Marker {
            center: (0.25, 0.25),
            shape: MarkerShape::Square,
            size: 0.1,
            fill: Some(Fill::solid(Color::BLUE)),
            stroke: None,
        });
        let svg = scene.to_svg(&unit_viewport());
        assert!(svg.contains("<circle "));
        assert!(svg.contains("<rect "));
    }

    #[test]
    fn svg_text_uses_anchor_and_baseline() {
        let mut scene = Scene::new();
        scene.push(Item::Text {
            pos: (0.5, 0.5),
            text: "hello".into(),
            style: TextStyle::new(0.1, Color::BLACK).align(HAlign::Center, VAlign::Middle),
        });
        let svg = scene.to_svg(&unit_viewport());
        assert!(svg.contains("text-anchor=\"middle\""));
        assert!(svg.contains("dominant-baseline=\"central\""));
        assert!(svg.contains(">hello</text>"));
    }

    #[test]
    fn svg_escapes_text_and_attrs() {
        let mut scene = Scene::new();
        scene.push(Item::Text {
            pos: (0.0, 0.0),
            text: "a < b & c > d".into(),
            style: TextStyle::new(0.1, Color::BLACK).family("\"weird\" <font>"),
        });
        let svg = scene.to_svg(&unit_viewport());
        assert!(svg.contains("a &lt; b &amp; c &gt; d"));
        assert!(svg.contains("&quot;weird&quot;"));
    }

    #[test]
    fn arrow_end_emits_triangle_at_b() {
        let mut scene = Scene::new();
        scene.push(Item::Segment {
            a: (0.0, 0.0),
            b: (1.0, 0.0),
            stroke: Stroke::solid(Color::BLACK, 0.05),
            arrow: Some(ArrowStyle::end(0.2)),
        });
        let svg = scene.to_svg(&unit_viewport());
        // Line + one arrowhead polygon = 2 main shape tags.
        assert_eq!(svg.matches("<line ").count(), 1);
        assert_eq!(svg.matches("<polygon ").count(), 1);
    }

    #[test]
    fn arrow_per_edge_emits_one_triangle_per_edge() {
        // Open polyline with 4 points → 3 edges → 3 arrowheads.
        let mut scene = Scene::new();
        scene.push(Item::Polyline {
            points: vec![(0.0, 0.0), (0.5, 0.5), (1.0, 0.0), (0.5, -0.5)],
            stroke: Stroke::solid(Color::BLACK, 0.02),
            arrow: Some(ArrowStyle::per_edge(0.1)),
        });
        let svg = scene.to_svg(&unit_viewport());
        assert_eq!(svg.matches("<polyline ").count(), 1);
        assert_eq!(svg.matches("<polygon ").count(), 3);
    }

    #[test]
    fn arrow_per_edge_on_polygon_wraps() {
        // Closed polygon with 3 vertices → 3 edges (cyclic) →
        // 3 arrowheads.
        let mut scene = Scene::new();
        scene.push(Item::Polygon {
            points: vec![(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)],
            fill: None,
            stroke: Some(Stroke::solid(Color::BLACK, 0.02)),
            arrow: Some(ArrowStyle::per_edge(0.1)),
        });
        let svg = scene.to_svg(&unit_viewport());
        // One polygon (the shape) + 3 arrowhead polygons.
        assert_eq!(svg.matches("<polygon ").count(), 4);
    }

    #[test]
    fn arrow_per_edge_mid_segment_centers_triangle() {
        // For a segment a=(0,0), b=(2,0) with size=0.4, the
        // mid-edge triangle should be centered on x=1. Tip is at
        // (1 + 0.2, 0) = (1.2, 0); base center at (0.8, 0).
        // Viewport is 100×100 over (0,0)-(2,0)... let's use a
        // wider viewport.
        let mut scene = Scene::new();
        scene.push(Item::Segment {
            a: (0.0, 0.0),
            b: (2.0, 0.0),
            stroke: Stroke::solid(Color::BLACK, 0.05),
            arrow: Some(ArrowStyle::per_edge_mid(0.4)),
        });
        let vp = Viewport::square_for(200, ((0.0, -1.0), (2.0, 1.0)), 0);
        let svg = scene.to_svg(&vp);
        // Just check we emitted exactly one triangle.
        assert_eq!(svg.matches("<polygon ").count(), 1);
        // Coordinates should land near the center of the canvas
        // (math x=1.2 → pixel x ≈ 120 for a 200px canvas spanning
        // x∈[0,2]).
        assert!(svg.contains("120"), "expected mid-tip near px 120, svg: {svg}");
    }

    #[test]
    fn svg_background_renders_first_rect() {
        let scene = Scene::new().with_background(Color::WHITE);
        let svg = scene.to_svg(&unit_viewport());
        // Background must be emitted before items.
        let bg_idx = svg.find("<rect ").expect("background rect");
        let body_idx = svg.find("</svg>").unwrap();
        assert!(bg_idx < body_idx);
        assert!(svg[bg_idx..body_idx].contains("fill=\"#ffffff\""));
    }
}
