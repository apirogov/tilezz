use tilezz::snake::{Snake, Turtle};
use tilezz::zz::ZZ12;

use plotters::prelude::*;

/// Given a list of (x,y) coordinates, return the bounds
/// ((min_x, min_y), (max_x, max_y)).
fn bounds(pts: &[(f64, f64)]) -> ((f64, f64), (f64, f64)) {
    let (x, y) = pts[0];
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (x, y, x, y);
    for (x, y) in pts {
        (min_x, min_y) = (min_x.min(*x), min_y.min(*y));
        (max_x, max_y) = (max_x.max(*x), max_y.max(*y));
    }
    ((min_x, min_y), (max_x, max_y))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s: Snake<ZZ12> = Snake::from(&[0, 5, 3, 4, 1, 4, 1, 1, 3, 5]);
    let l = s.to_polyline_f64(&Turtle::default());
    let ((min_x, min_y), (max_x, max_y)) = bounds(&l);
    let (w, h) = (max_x - min_x, max_y - min_y);
    let half_d = 0.5 * if w >= h { w - h } else { h - w };
    let pad = 1_f64;
    let pad_x = pad + if w < h { half_d } else { 0_f64 };
    let pad_y = pad + if h < w { half_d } else { 0_f64 };

    // for p in l

    let root = BitMapBackend::new("test.png", (1000, 1000)).into_drawing_area();
    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);
    // After this point, we should be able to construct a chart context
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        .caption("My first snake", ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(40)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(
            (min_x - pad_x)..(max_x + pad_x),
            (min_y - pad_y)..(max_y + pad_y),
        )?;

    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        // .x_labels(5)
        // .y_labels(5)
        .draw()?;

    // And we can draw something in the drawing area
    chart.draw_series(LineSeries::new(l.clone(), &RED))?;
    // Similarly, we can draw point series
    chart.draw_series(PointSeries::of_element(l.clone(), 2, &RED, &|c, s, st| {
        return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
            + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
            + Text::new(format!("{:.3},{:.3}", c.0, c.1), (10, 0), ("sans-serif", 10).into_font());
    }))?;
    root.present()?;
    Ok(())
}
