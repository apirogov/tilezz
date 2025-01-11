use plotters::prelude::*;
use tilezz::plotters::{plot_tile, tile_chart, TileStyle};
use tilezz::rat::Rat;
use tilezz::snake::constants::spectre;
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
    let spectre_mod = spectre.reflected().cycle(5);

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

    let r1 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 0, 3, 3, -3, -3, 3, 3, 0, 3]);
    let upts = r1.to_polyline_f64(Turtle::default());
    let r2 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 3, 0, 0, 3, 3]);
    let upts2 = r2.to_polyline_f64(Turtle::new(ZZ12::from((0, 3)), 0));

    let r3 = r1.glue_unchecked((0, 3), &r2);
    // let r3 = r1.glue_unchecked((1, 4), &r2).cycle(-3);
    let upts3 = r3.to_polyline_f64(Turtle::new(ZZ12::from((5, 3)), 3));

    let s = TileStyle::default().with_fill(TRANSPARENT.into());
    let tiles2 = vec![
        (upts.as_slice(), &s),
        (upts2.as_slice(), &s),
        (upts3.as_slice(), &s),
    ];
    let (mut c2, plot_tiles2) = tile_chart(&mut ChartBuilder::on(&right), tiles2.as_slice());
    c2.configure_mesh().draw().unwrap();
    right.fill(&WHITE).unwrap();
    plot_tiles2(&mut c2);

    // test unchecked snake rendering (self-intersecting star)
    // let unchecked = Snake::<ZZ12>::from_slice_unchecked(&[5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]);
    // plot_tile(
    //     &right,
    //     upts2.as_slice(),
    //     &TileStyle::default().with_fill(TRANSPARENT.into()),
    // );

    // make sure the image is actually updated
    root.present().unwrap();
}
