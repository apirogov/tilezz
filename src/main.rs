use plotters::prelude::*;
use tilezz::plotters::{plot_tile, TileStyle};
use tilezz::rat::Rat;
use tilezz::snake::constants::spectre;
use tilezz::snake::Turtle;
use tilezz::zz::ZZ12;
use tilezz::zzbase::ZZBase;

fn my_custom_style<'a>() -> TileStyle<'a> {
    // FIXME: add builder pattern for convenient tile styling?
    let mut st = TileStyle::default();
    st.label_font = st.label_font.color(&BLUE);
    st.border_style = st.border_style.stroke_width(5);
    st.fill_style.color = TRANSPARENT;
    st.node_size = 0;
    st.node_font = st.node_font.color(&RED);
    st.node_font.font = st.node_font.font.resize(40.);
    return st;
}

fn main() {
    let root = BitMapBackend::new("test.png", (1000, 1000)).into_drawing_area();
    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);

    let tile: Rat<ZZ12> = Rat::new(&spectre());
    let pts = tile.to_polyline_f64(Turtle::new(0.into(), -2));
    plot_tile(&root, &pts, &TileStyle::default().with_label("The Spectre"));

    let _pts2 = tile
        .reflect()
        .cycle(5)
        .to_polyline_f64(Turtle::new(ZZ12::from(3) + ZZ12::one_i().scale(3), 0));

    // FIXME: render multiple tiles in one plot (take a sequence of pairs of points and styles)
    // plot_tile(&root, &pts2, &my_custom_style());

    root.present().unwrap();
}
