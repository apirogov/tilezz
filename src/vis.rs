use plotters::{
    chart::ChartContext,
    coord::types::RangedCoordf64,
    prelude::{
        Cartesian2d, ChartBuilder, Circle, DrawingBackend, EmptyElement, LineSeries, PointSeries,
        Text, BLACK, RED,
    },
    series::AreaSeries,
    style::{
        text_anchor::{HPos, Pos, VPos},
        Color, IntoTextStyle, RGBColor, ShapeStyle, TextStyle, YELLOW,
    },
};

/// Given a list of (x,y) coordinates, return the bounds
/// ((min_x, min_y), (max_x, max_y)).
pub fn bounds(pts: &[(f64, f64)]) -> ((f64, f64), (f64, f64)) {
    let (x, y) = pts[0];
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (x, y, x, y);
    for (x, y) in pts {
        (min_x, min_y) = (min_x.min(*x), min_y.min(*y));
        (max_x, max_y) = (max_x.max(*x), max_y.max(*y));
    }
    ((min_x, min_y), (max_x, max_y))
}

/// Given a list of (x,y) coordinates, returns the bounds
/// ((min_x, min_y), (max_x, max_y)) of a square
/// centered on and including all the points.
pub fn viewport(pts: &[(f64, f64)]) -> ((f64, f64), (f64, f64)) {
    let ((min_x, min_y), (max_x, max_y)) = bounds(pts);
    let (w, h) = (max_x - min_x, max_y - min_y);
    let half_d = 0.5 * if w >= h { w - h } else { h - w };
    let pad = 1_f64;
    let pad_x = pad + if w < h { half_d } else { 0_f64 };
    let pad_y = pad + if h < w { half_d } else { 0_f64 };

    (
        ((min_x - pad_x), (min_y - pad_y)),
        ((max_x + pad_x), (max_y + pad_y)),
    )
}

pub fn plot_tile<'a, 'b, DB: DrawingBackend>(
    cb: &mut ChartBuilder<DB>,
    tile: &[(f64, f64)],
) -> ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>> {
    // prepare coordinate system

    // square view centered on tile
    let ((min_x, min_y), (max_x, max_y)) = viewport(tile);

    let mut chart = cb.build_cartesian_2d(min_x..max_x, min_y..max_y).unwrap();
    chart.configure_mesh().draw().unwrap();

    // draw segments
    let lines: LineSeries<DB, (f64, f64)> = LineSeries::new(tile.to_vec(), &BLACK);
    let area = AreaSeries::new(tile.to_vec(), 0., &YELLOW.mix(0.2));
    chart.draw_series(lines).unwrap();
    chart.draw_series(area).unwrap();

    // draw vertex labels (i.e. show indices of vertices)
    let mut ixd_pts: Vec<(usize, (f64, f64))> =
        tile.iter().enumerate().map(|(i, c)| (i, *c)).collect();
    ixd_pts.pop(); // drop last point (same as first)

    let mk_font = |size, clr| {
        ("sans-serif", size)
            .with_color(clr)
            .with_anchor::<RGBColor>(Pos::new(HPos::Center, VPos::Center))
            .into_text_style(chart.plotting_area())
    };

    let tile_lbl = "hello";
    let tile_lbl_style: TextStyle = mk_font(32, BLACK);
    let node_lbl_style: TextStyle = mk_font(16, BLACK);
    let node_shape_style: ShapeStyle = RED.filled();

    let node_func = |(idx, c): (usize, (f64, f64)), s, st: ShapeStyle| {
        return EmptyElement::at(c)
            + Circle::new((0, 0), s, st)
            + Text::new(format!("{idx}"), (0, 0), &node_lbl_style);
    };
    let nodes = PointSeries::of_element(ixd_pts.clone(), 10, node_shape_style, &node_func);
    chart.draw_series(nodes).unwrap();

    let tile_lbl_func = |c, _, _| {
        return EmptyElement::at(c) + Text::new(format!("{tile_lbl}"), (0, 0), &tile_lbl_style);
    };
    let tile_center = (min_x + (max_x - min_x) / 2., min_y + (max_y - min_y) / 2.);
    let tile_lbl_series = PointSeries::of_element(vec![tile_center], 20, BLACK, &tile_lbl_func);
    chart.draw_series(tile_lbl_series).unwrap();

    chart
}
