//! Helper fuctions specifically for rendering with plotters

use num_traits::Zero;

use plotters::prelude::*;
use plotters::{
    coord::types::RangedCoordf64,
    coord::Shift,
    style::text_anchor::{HPos, Pos, VPos},
};

use crate::plotutils::{tile_bounds, tile_centroid, tile_viewport};

/// Metadata and style to be applied to a rendered tile.
// TODO: could allow customizing border rendering with a pair of funcs (lines, segments)
// and provide a outline style builder to simplify describing a style.
#[derive(Clone)]
pub struct TileStyle<'a> {
    pub fill_style: ShapeStyle,
    pub border_style: ShapeStyle,

    pub node_style: ShapeStyle,
    pub node_size: i32,

    pub node_labels: bool,
    pub node_zero_only: bool,
    pub node_font: TextStyle<'a>,

    pub label: Option<String>,
    pub label_font: TextStyle<'a>,
}

impl<'a> TileStyle<'a> {
    pub fn with_label(mut self, lbl: &str) -> Self {
        self.label = Some(lbl.to_string());
        self
    }

    pub fn with_fill(mut self, style: ShapeStyle) -> Self {
        self.fill_style = style;
        self
    }

    pub fn with_border(mut self, style: ShapeStyle) -> Self {
        self.border_style = style;
        self
    }

    // TODO: add more builder-pattern methods for convenient tile styling?
}

impl<'a> Default for TileStyle<'a> {
    fn default() -> Self {
        Self {
            fill_style: YELLOW.mix(0.2).into(),
            border_style: BLACK.into(),

            node_style: RED.filled().into(),
            node_size: 10,

            node_labels: true,
            node_zero_only: false,
            node_font: ("sans-serif", 16, FontStyle::Bold).into_font().color(&BLUE),

            label: None,
            label_font: ("sans-serif", 32, FontStyle::Normal)
                .into_font()
                .color(&BLACK),
        }
    }
}

/// Plot a sequence of points of a given size using the given style.
pub fn plot_points<'a, DB: DrawingBackend, F: Fn((f64, f64)) -> String>(
    chart: &mut ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    pts: &[(f64, f64)],
    label_func: F,
    size: i32,
    style: ShapeStyle,
) {
    let node_func = |c: (f64, f64), s: i32, st: ShapeStyle| {
        EmptyElement::at(c) + Circle::<(i32, i32), i32>::new((0.into(), 0.into()), s.into(), st)
    };
    let nodes = PointSeries::of_element(pts.to_vec(), size, style, &node_func);
    chart.draw_series(nodes).unwrap();

    let centered = Pos::new(HPos::Center, VPos::Center);
    let mut txtstyle: TextStyle = TileStyle::default().node_font;
    txtstyle.font = txtstyle.font.resize(size as f64 * 2.0);
    txtstyle = txtstyle.color(&WHITE);
    txtstyle = txtstyle.pos(centered);
    let node_lbl_func = |c: (f64, f64), _: i32, _: ShapeStyle| {
        EmptyElement::at(c) + Text::new(format!("{}", label_func(c)), (0, 0), &txtstyle)
    };
    let node_lbls = PointSeries::of_element(pts.to_vec(), size, style, &node_lbl_func);
    chart.draw_series(node_lbls).unwrap();
}

pub fn plot_seg<'a, DB: DrawingBackend, S: Into<ShapeStyle>>(
    chart: &mut ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    (x, y): ((f64, f64), (f64, f64)),
    style: S,
) {
    let line = LineSeries::new(vec![x, y], style);
    chart.draw_series(line).unwrap();
}

/// Render a tile into a prepared chart.
fn plot_tile_into<'a, DB: DrawingBackend>(
    chart: &mut ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    tile: &[(f64, f64)],
    style: &TileStyle,
) {
    // we want to center all rendered texts
    let centered = Pos::new(HPos::Center, VPos::Center);

    // draw segments (border and/or fill)
    let area =
        AreaSeries::new(tile.to_vec(), 0., style.fill_style).border_style(style.border_style);
    chart.draw_series(area).unwrap();

    // draw vertex marker circles
    if !style.node_size.is_zero() {
        let mut pts = tile.to_vec();
        if style.node_zero_only {
            pts.truncate(1);
        }

        plot_points(
            chart,
            pts.as_slice(),
            |_| "".to_string(),
            style.node_size,
            style.node_style,
        );
    }

    // draw vertex labels (i.e. show indices of vertices)
    if style.node_labels {
        let mut ixd_pts: Vec<(usize, (f64, f64))> =
            tile.iter().enumerate().map(|(i, c)| (i, *c)).collect();
        ixd_pts.pop(); // drop last point (same as first)
        if style.node_zero_only {
            ixd_pts.truncate(1);
        }

        let node_lbl_style: TextStyle = style.node_font.pos(centered);
        let node_lbl_func = |(idx, c): (usize, (f64, f64)), _: i32, _: ShapeStyle| {
            EmptyElement::at(c) + Text::new(format!("{idx}"), (0, 0), &node_lbl_style)
        };
        let node_lbls = PointSeries::of_element(
            ixd_pts,
            style.node_size.into(),
            style.node_style,
            &node_lbl_func,
        );
        chart.draw_series(node_lbls).unwrap();
    }

    // render tile label
    if let Some(tile_lbl) = style.label.as_ref() {
        let tile_lbl_style: TextStyle = style.label_font.pos(centered);
        let tile_lbl_func = |c, _, _| {
            return EmptyElement::at(c) + Text::new(format!("{tile_lbl}"), (0, 0), &tile_lbl_style);
        };
        let tile_lbl_series =
            PointSeries::of_element(vec![tile_centroid(tile)], 20, BLACK, &tile_lbl_func);
        chart.draw_series(tile_lbl_series).unwrap();
    }
}

/// Prepare a chart with a coordinate system of sufficient size for the given tiles.
/// Returns the chart and a callback function to render the tiles.
/// This can be used to e.g. customize how the mesh is rendered.
pub fn tile_chart<'a, DB: DrawingBackend>(
    cb: &mut ChartBuilder<DB>,
    tiles: &'a [(&[(f64, f64)], &TileStyle)],
) -> (
    ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    impl FnOnce(&mut ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>) -> (),
) {
    // compute square view centered on tile(s)
    let all_points = tiles.iter().map(|(tile, _)| *tile).flatten();
    let ((min_x, min_y), (max_x, max_y)) = tile_viewport(tile_bounds(all_points));

    let chart = cb.build_cartesian_2d(min_x..max_x, min_y..max_y).unwrap();
    let t = tiles.to_vec();

    // FIXME: as I did not figure out how to make a zero-argument closure,
    // I have to force the user to give me the chart again.
    // To ensure that the chart is at least compatible, the ranges are checked.
    let x_range = chart.x_range();
    let y_range = chart.y_range();

    // Plotting function to be called when the mesh and chart are set up as desired.
    let plot_func = move |c: &mut plotters::chart::ChartContext<
        'a,
        DB,
        plotters::prelude::Cartesian2d<RangedCoordf64, RangedCoordf64>,
    >| {
        assert_eq!(c.x_range(), x_range);
        assert_eq!(c.y_range(), y_range);
        for (tile, style) in t {
            plot_tile_into(c, tile, style);
        }
    };

    (chart, plot_func)
}

/// Plot a sequence of tiles into a chart that is constructed using the given chart builder.
pub fn plot_tiles<'a, DB: DrawingBackend>(
    cb: &mut ChartBuilder<DB>,
    tiles: &'a [(&[(f64, f64)], &TileStyle)],
) -> ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>> {
    let (mut chart, plot_func) = tile_chart(cb, tiles);

    chart
        .configure_mesh()
        .disable_mesh()
        .disable_axes()
        .draw()
        .unwrap();

    plot_func(&mut chart);
    chart
}

/// Plot a single tile into a bare-bones chart. Check the other plotting functions for more flexibility.
pub fn plot_tile<'a, DB: DrawingBackend>(
    da: &DrawingArea<DB, Shift>,
    tile: &[(f64, f64)],
    style: &TileStyle,
) {
    plot_tiles(&mut ChartBuilder::on(&da), &[(tile, style)]);
}
