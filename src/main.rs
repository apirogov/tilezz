use plotters::prelude::*;
use tilezz::plotters::{plot_tile, tile_chart, TileStyle};
use tilezz::rat::Rat;
use tilezz::snake::constants::spectre;
use tilezz::snake::Snake;
use tilezz::snake::Turtle;
use tilezz::zz::ZZ12;

// Return a customized tile style.
fn my_custom_style<'a>() -> TileStyle<'a> {
    let mut st = TileStyle::default();
    st.label_font = st.label_font.color(&RED);
    st.border_style = st.border_style.stroke_width(5);
    st.fill_style.color = TRANSPARENT;
    st.node_zero_only = true;
    st.node_size = 10;
    st.node_style = BLUE.filled();
    st.node_labels = false;
    return st;
}

fn main() {
    let root = BitMapBackend::new("test.png", (2000, 1000)).into_drawing_area();
    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);

    // get a spectre tile
    let spectre: Rat<ZZ12> = Rat::try_from(&spectre()).unwrap();
    // glue a copy of the spectre to itself to get a mystic
    let mystic = spectre.glue((2, 0), &spectre);
    // take the mirror image and shift the starting point
    let spectre_mod = spectre.reflect().cycle(5);

    // get concrete points and define styles for the tiles
    let t1 = spectre.to_polyline_f64(Turtle::new(0.into(), -2));
    let s1 = TileStyle::default().with_label("The Spectre");

    let t2 = mystic.to_polyline_f64(Turtle::new(ZZ12::from((2, 3)), -3));
    let s2 = TileStyle::default().with_label("The Mystic");

    let t3 = spectre_mod.to_polyline_f64(Turtle::new(0.into(), -1));
    let s3 = my_custom_style().with_label("Customized Spectre");

    // TODO: implement tile patch structure that keeps the original tiles accessible/renderable individually

    // split the drawing area
    let (left, right) = root.split_horizontally(1000);

    // plot two tiles in the chart on the left
    let tiles = vec![(t1.as_slice(), &s1), (t2.as_slice(), &s2)];
    let (mut c1, plot_tiles) = tile_chart(
        &mut ChartBuilder::on(&left)
            .caption("Custom Tile Plot", ("sans-serif", 40).into_font())
            .x_label_area_size(20)
            .y_label_area_size(40),
        tiles.as_slice(),
    );
    c1.configure_mesh().draw().unwrap();
    plot_tiles(&mut c1);

    // plot a tile in the chart on the right
    right.fill(&GREEN.mix(0.1)).unwrap();
    plot_tile(&right, t3.as_slice(), &s3);

    // test unchecked snake rendering
    let unchecked = Snake::<ZZ12>::from_slice_unchecked(&[5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]);
    let unchecked_pts = unchecked.to_polyline_f64(&Turtle::default());
    right.fill(&WHITE).unwrap();
    plot_tile(
        &right,
        unchecked_pts.as_slice(),
        &TileStyle::default().with_fill(TRANSPARENT.into()),
    );

    // make sure the image is actually updated
    root.present().unwrap();
}
