use num_complex::ComplexFloat;
use std::collections::HashSet;

use plotters::prelude::*;

use tilezz::angles::to_abs_seq;
use tilezz::cyclotomic::*;
use tilezz::plotters::plot_points;
use tilezz::plotutils::tile_bounds;
use tilezz::snake::constants::spectre;
use tilezz::snake::Snake;

fn main() {
    // let root = BitMapBackend::new("test.png", (1000, 1000))
    let w: u32 = 3000;
    let root = BitMapBackend::new("scratch/test.png", (w, w))
        // let root = BitMapBackend::new("test.png", (w, (w as f64 * ((3.1415 / 3.0).sin())) as u32))
        .into_drawing_area();

    let _ = root.fill(&WHITE);
    let root = root.margin(10, 10, 10, 10);

    let num_rounds = 1 + 96;

    type ZZ = ZZ12;
    type Z = <ZZ as ZZBase>::Real;

    let sn: Snake<ZZ> = spectre();
    println!("{:?} -> {:?}", sn.angles(), to_abs_seq::<ZZ>(sn.angles()));

    let origin: ZZ = 0.into();
    let one: ZZ = 1.into();
    let one_i: ZZ = ZZ::one_i();
    let one_i_p_one: ZZ = one + one_i;

    let mut edges: HashSet<(ZZ, ZZ)> = HashSet::new();
    let mut visited: HashSet<ZZ> = HashSet::new();
    visited.insert(origin);
    visited.insert(one);
    visited.insert(one_i);
    visited.insert(one_i_p_one);

    let mut round_pts: Vec<Vec<ZZ>> = Vec::new();
    // let mut pts_ang: Vec<Vec<i8>> = Vec::new();

    round_pts.push(vec![origin, one, one_i_p_one, one_i]);
    // pts_ang.push(vec![0, 3, 6, 9]);
    for i in 1..=num_rounds {
        println!("Computing points at distance {i}");

        let last = round_pts.last().unwrap();
        // let last_ang = pts_ang.last().unwrap();

        let mut curr: Vec<ZZ> = Vec::new();
        let mut curr_ang: Vec<i8> = Vec::new();

        for (/* idx */ _, p) in last.iter().enumerate() {
            for i in 0..ZZ::turn() {
                let new_ang = (/* last_ang[idx] + */i) % ZZ::turn();

                let mut p_dest = *p + ZZ::unit(new_ang);
                let (mut p_x, mut p_y) = p_dest.re_im();
                if (p_x - Z::zero()).is_negative() {
                    p_x = p_x + Z::one();
                }
                if (Z::one() - p_x).is_negative() {
                    p_x = p_x - Z::one();
                }
                if (p_y - Z::zero()).is_negative() {
                    p_y = p_y + Z::one();
                }
                if (Z::one() - p_y).is_negative() {
                    p_y = p_y - Z::one();
                }
                p_dest = ZZ::from(p_x) + ZZ::one_i() * ZZ::from(p_y);

                // if norm(&p_dest).complex().abs() > 4. {
                //     continue;
                // }
                // if !angle_between(&p_dest, (&ZZ::one(), &ZZ::ccw())) {
                //     continue; //outside of slice
                // }

                if !visited.contains(&p_dest) {
                    // let mut already_have_angle = false;
                    // for p in visited.iter() {
                    //     if p.is_zero() {
                    //         continue;
                    //     }
                    //     let ang_sgn = angle_signum(&p_dest, &p);
                    //     // println!("{ang_sgn:?}");
                    //     if ang_sgn.is_zero() {
                    //         // already_have_angle = true;
                    //         break;
                    //     }
                    // }
                    // if already_have_angle {
                    //     continue;
                    // }

                    visited.insert(p_dest);

                    curr.push(p_dest);
                    curr_ang.push(new_ang);

                    edges.insert((*p, p_dest));
                    edges.insert((p_dest, *p));
                }
            }
        }
        println!("Number of new points: {}", curr.len());
        round_pts.push(curr);
        // pts_ang.push(curr_ang);
    }
    let num_pts: usize = round_pts.iter().map(|v| v.len()).sum();
    println!("Total points: {num_pts}");

    let render_pts: Vec<Vec<(f64, f64)>> = round_pts
        .iter()
        .map(|v| {
            v.iter()
                .map(|p| {
                    let c = p.complex64();
                    (c.re, c.im)
                })
                .collect()
        })
        .collect();

    println!("Computing bounds");
    let mut bound_pts: Vec<(f64, f64)> = Vec::new();
    for v in render_pts.iter() {
        let (x, y) = tile_bounds(v.iter());
        bound_pts.push(x);
        bound_pts.push(y);
    }
    let ((x_min, y_min), (x_max, y_max)) = tile_bounds(bound_pts.as_slice());

    let mut styles: Vec<(i32, ShapeStyle)> = Vec::new();
    for i in 0..=num_rounds {
        let step = 255 / num_rounds;
        let c: RGBColor = RGBColor(step / 2 * (num_rounds - i) + 128, step / 8 * i, step * i);
        // styles.push((2, c.filled().into()));
        styles.push(((30. * (0.99.powi(i as i32))) as i32, c.filled().into()));
    }

    println!("Preparing chart {x_min} {x_max} {y_min} {y_max}");
    // TODO: preserve aspect ratio helper (create sub draw area to avoid streching of coorinate system)
    let mut chart = ChartBuilder::on(&root)
        .x_label_area_size(20)
        .y_label_area_size(40)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .unwrap();
    chart.configure_mesh().draw().unwrap();

    // println!("plotting {} edges", edges.len());
    // for (x, y) in edges {
    //     let (pt_x, pt_y) = (x.complex(), y.complex());
    //     plot_seg(&mut chart, ((pt_x.re, pt_x.im), (pt_y.re, pt_y.im)), &BLACK);
    // }
    // for i in 0..=num_rounds {
    for i in num_rounds..=num_rounds {
        // for i in [8_usize, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96] {
        let (sz, st) = styles[i as usize];
        println!("plotting level {i}");
        plot_points(
            &mut chart,
            render_pts[i as usize].as_slice(),
            |_| format!("{i}"),
            sz,
            st,
        );
    }

    // make sure the image is actually updated
    root.present().unwrap();
}
