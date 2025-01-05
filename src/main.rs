use plotters::prelude::*;
use tilezz::plotters::{plot_tile_with, TileStyle};
use tilezz::rat::Rat;
use tilezz::snake::constants::spectre;
use tilezz::snake::Turtle;
use tilezz::zz::ZZ12;

// fn my_custom_style<'a>() -> TileStyle<'a> {
//     // FIXME: add builder pattern for convenient tile styling?
//     let mut st = TileStyle::default();
//     st.label_font = st.label_font.color(&BLUE);
//     st.border_style = st.border_style.stroke_width(5);
//     st.fill_style.color = TRANSPARENT;
//     st.node_size = 0;
//     st.node_font = st.node_font.color(&RED);
//     st.node_font.font = st.node_font.font.resize(40.);
//     return st;
// }

fn main() {
    let root = BitMapBackend::new("test.png", (1000, 1000)).into_drawing_area();
    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);

    // get a spectre tile
    let spectre: Rat<ZZ12> = Rat::new(&spectre());
    // glue a copy of it to itself to get a mystic
    let mystic = spectre.glue((2, 0), &spectre);

    // TODO: implement tile patch structure that keeps the original tiles accessible
    // let pts = mystic.to_polyline_f64(Turtle::new(0.into(), -2));
    // plot_tile(&root, &pts, &TileStyle::default().with_label("The Spectre"));

    // TODO: refactor snake to allow representing sequences self-intersections (optionally)
    // TODO: provide a way to initialize snake from possibly invalid slice without panic (use Result)

    let pts2 = mystic.to_polyline_f64(Turtle::new(0.into(), -3));
    plot_tile_with(
        &mut ChartBuilder::on(&root)
            .caption("Custom Tile Plot", ("sans-serif", 40).into_font())
            .x_label_area_size(20)
            .y_label_area_size(40),
        &pts2,
        &TileStyle::default().with_label("The Mystic"),
    );

    // let pts3 = spectre
    //     .reflect()
    //     .cycle(5)
    //     .to_polyline_f64(Turtle::new(ZZ12::from(3) + ZZ12::one_i().scale(3), 0));
    // plot_tile(&root, &pts3, &my_custom_style());

    // FIXME: render multiple tiles in one plot (take a sequence of pairs of points and styles)
    // FIXME: figure out how to align multiple charts in one drawing area

    root.present().unwrap();
}
