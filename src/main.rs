use plotters::prelude::*;
use tilezz::plotters::plot_tile;
use tilezz::snake::constants::spectre;
use tilezz::snake::{Snake, Turtle};
use tilezz::zz::ZZ12;

fn main() {
    let s: Snake<ZZ12> = spectre();

    let tile = s.to_polyline_f64(&Turtle::default());

    let root = BitMapBackend::new("test.png", (1000, 1000)).into_drawing_area();
    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);

    plot_tile(
        &mut ChartBuilder::on(&root)
            .caption("Spectre tile", ("sans-serif", 40).into_font())
            .x_label_area_size(20)
            .y_label_area_size(40),
        &tile,
    );

    root.present().unwrap();
}
