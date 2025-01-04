use plotters::prelude::*;
use tilezz::plotters::plot_tile;
use tilezz::snake::constants::spectre;
use tilezz::snake::{Snake, Turtle};
use tilezz::traits::Ccw;
use tilezz::zz::{ZZ12, ZZ60};
use tilezz::zzbase::{GInt, ZZBase};

fn main() {
    let p = ZZ60::ccw();
    println!("{}", p * p);

    //                               x  x  x <- incorrect result
    // let v1 = ZZ30::new(&[o, o, o, z, z, z, o, o]);
    for i in 0..8 {
        println!("------------");
        for j in 0..8 {
            println!("----");
            println!("l={i} r={j}");
            let mut vec1: Vec<GInt> = vec![0.into(); 8];
            let mut vec2: Vec<GInt> = vec![0.into(); 8];
            vec1[i] = 1.into();
            vec2[j] = 1.into();
            let v1 = ZZ60::new(vec1.as_slice());
            let v2 = ZZ60::new(vec2.as_slice());

            println!("{} {}", (v1.complex() * v2.complex()), (v1 * v2).complex());
            println!("{}\n{}\n{}", v1, v2, v1 * v2);
        }
    }

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
