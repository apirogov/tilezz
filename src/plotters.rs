use crate::plotutils::{tile_mass_center, tile_viewport};

use num_traits::Zero;
use plotters::prelude::*;
use plotters::{
    coord::{types::RangedCoordf64, Shift},
    style::text_anchor::{HPos, Pos, VPos},
};

/// Metadata and style to be applied to a rendered tile.
pub struct TileStyle<'a> {
    pub fill_style: ShapeStyle,
    pub border_style: ShapeStyle,

    pub node_style: ShapeStyle,
    pub node_size: i32,

    pub node_labels: bool,
    pub node_font: TextStyle<'a>,

    pub label: Option<String>,
    pub label_font: TextStyle<'a>,
}

impl<'a> TileStyle<'a> {
    pub fn with_label(mut self, lbl: &str) -> Self {
        self.label = Some(lbl.to_string());
        self
    }
}

impl<'a> Default for TileStyle<'a> {
    fn default() -> Self {
        Self {
            fill_style: YELLOW.mix(0.2).into(),
            border_style: BLACK.into(),

            node_style: RED.filled().into(),
            node_size: 10,

            node_labels: true,
            node_font: ("sans-serif", 16, FontStyle::Bold).into_font().color(&BLUE),

            label: None,
            label_font: ("sans-serif", 32, FontStyle::Normal)
                .into_font()
                .color(&BLACK),
        }
    }
}

/// Plot a tile into a chart based on the pre-configured custom
/// chart builder (which can be used to e.g. add axes and a caption).
pub fn plot_tile_with<'a, 'b, DB: DrawingBackend>(
    cb: &mut ChartBuilder<DB>,
    tile: &[(f64, f64)],
    style: &TileStyle,
) -> ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>> {
    // we want to center all rendered texts
    let centered = Pos::new(HPos::Center, VPos::Center);
    // square view centered on tile
    let ((min_x, min_y), (max_x, max_y)) = tile_viewport(tile);

    // prepare coordinate system
    let mut chart = cb.build_cartesian_2d(min_x..max_x, min_y..max_y).unwrap();
    chart.configure_mesh().draw().unwrap();

    // draw segments (border and/or fill)
    let area =
        AreaSeries::new(tile.to_vec(), 0., style.fill_style).border_style(style.border_style);
    chart.draw_series(area).unwrap();

    // draw vertex marker circles
    if !style.node_size.is_zero() {
        let node_func = |c: (f64, f64), s: i32, st: ShapeStyle| {
            EmptyElement::at(c) + Circle::<(i32, i32), i32>::new((0.into(), 0.into()), s.into(), st)
        };
        let nodes =
            PointSeries::of_element(tile.to_vec(), style.node_size, style.node_style, &node_func);
        chart.draw_series(nodes).unwrap();
    }

    // draw vertex labels (i.e. show indices of vertices)
    if style.node_labels {
        let node_lbl_style: TextStyle = style.node_font.pos(centered);
        let mut ixd_pts: Vec<(usize, (f64, f64))> =
            tile.iter().enumerate().map(|(i, c)| (i, *c)).collect();
        ixd_pts.pop(); // drop last point (same as first)
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
            PointSeries::of_element(vec![tile_mass_center(&tile)], 20, BLACK, &tile_lbl_func);
        chart.draw_series(tile_lbl_series).unwrap();
    } else {
    };

    chart
}

/// Plot a chart depicting a tile using default chart settings.
pub fn plot_tile<'a, 'b, DB: DrawingBackend>(
    da: &DrawingArea<DB, Shift>,
    tile: &[(f64, f64)],
    style: &TileStyle,
) -> ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>> {
    plot_tile_with(&mut ChartBuilder::on(da), &tile, &style)
}
