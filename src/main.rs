use std::collections::HashSet;
use std::io::{stdout, Write};

use plotters::coord::Shift;
use plotters::prelude::*;

use tilezz::cyclotomic::geometry::point_mod_rect;
use tilezz::cyclotomic::*;
use tilezz::plotters::{plot_points, rainbow};
use tilezz::plotutils::{tile_bounds, P64, R64};

/// Compute the levels of points reachable within the unit square from any Gaussian integer in n steps.
fn explore<ZZ: ZZType + HasZZ4>(n: usize, mod_unit_square: bool) -> Vec<Vec<ZZ>>
where
    <ZZ as IsComplex>::Field: From<(<ZZ as IsRingOrField>::Real, <ZZ as IsRingOrField>::Real)>,
{
    // we start at the corners of the unit square
    let start_pts: &[ZZ] = if mod_unit_square {
        &[ZZ::zero(), ZZ::one(), ZZ::one_i(), ZZ::one() + ZZ::one_i()]
    } else {
        &[ZZ::zero()]
    };

    let mut visited: HashSet<ZZ> = HashSet::from_iter(start_pts.to_vec().into_iter());
    let mut round_pts: Vec<Vec<ZZ>> = Vec::new();
    round_pts.push(start_pts.to_vec());

    let unit_square: (ZZ, ZZ) = ((0, 0).into(), (1, 1).into());
    // in each round, we go one unit step in every possible direction
    for i in 1..=n {
        let last = round_pts.last().unwrap();
        let mut curr: Vec<ZZ> = Vec::new();
        for p in last.iter() {
            for i in 0..ZZ::turn() {
                let mut p_dest: ZZ = *p + ZZ::unit(i);
                if mod_unit_square {
                    // normalize back into unit square (if enabled)
                    p_dest = point_mod_rect(&p_dest, &unit_square).coerce_ring();
                }
                let p_dest = p_dest;

                if !visited.contains(&p_dest) {
                    visited.insert(p_dest);
                    curr.push(p_dest);
                }
            }
        }
        print!("{}{}", curr.len(), if i == n { "" } else { "+" }); // print number of new points
        stdout().flush().unwrap();
        if i % 10 == 0 {
            println!("")
        }
        round_pts.push(curr);
    }
    let num_pts: usize = round_pts.iter().map(|v| v.len()).sum();
    println!("\n= {num_pts}");
    return round_pts;
}

/// Helper function to plot the points with given settings into a drawing area.
pub fn render<'a, DB: DrawingBackend>(
    da: &DrawingArea<DB, Shift>,

    ((x_min, y_min), (x_max, y_max)): R64,
    point_levels: &[Vec<P64>],
    level_styles: &[(i32, ShapeStyle)],
    offset: usize,
    stride: u32,
) {
    // prepare coordinate system
    let mut chart = ChartBuilder::on(&da)
        .x_label_area_size(20)
        .y_label_area_size(40)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .unwrap();
    chart.configure_mesh().draw().unwrap();
    // plot points level by level with the correct style
    for i in (0..point_levels.len()).step_by(stride as usize).rev() {
        let (sz, st) = level_styles[i as usize];
        plot_points(
            &mut chart,
            point_levels[i as usize].as_slice(),
            |_| format!("{}", i + offset),
            sz,
            st,
        );
    }
}

/// Generate combinations of point sizes and colors
fn get_styles(n: usize, sz: Option<i32>) -> Vec<(i32, ShapeStyle)> {
    // get a rainbow color gradient
    let colors = rainbow(n as u32 + 1, 1.);
    (0..=n)
        .collect::<Vec<_>>()
        .iter()
        .map(|i| {
            (
                match sz {
                    Some(size) => size,
                    None => (40. * (0.99_f64.powi(*i as i32))) as i32,
                },
                colors[*i as usize].filled().into(),
            )
        })
        .collect()
}

fn prepare_render<ZZ: ZZType + HasZZ4>(
    num_rounds: usize,
    mod_unit_square: bool,
    pt_size: Option<i32>,
) -> (Vec<Vec<P64>>, R64, Vec<(i32, ShapeStyle)>)
where
    <ZZ as IsComplex>::Field: From<(<ZZ as IsRingOrField>::Real, <ZZ as IsRingOrField>::Real)>,
{
    let points: Vec<Vec<P64>> = explore::<ZZ>(num_rounds, mod_unit_square)
        .iter()
        .map(|v| v.iter().map(|p| p.xy()).collect())
        .collect();

    let bounds = tile_bounds(
        points
            .iter()
            .map(|p| <(P64, P64) as Into<[P64; 2]>>::into(tile_bounds(p.iter())).into_iter())
            .flatten()
            .collect::<Vec<P64>>()
            .as_slice(),
    );
    let styles = get_styles(points.len(), pt_size);

    (points, bounds, styles)
}

fn main() {
    let w: u32 = 2000;
    let root = BitMapBackend::gif("scratch/test.gif", (w, w), 250)
        .unwrap()
        .into_drawing_area();
    // let root = BitMapBackend::new("test.png", (w, (w as f64 * ((3.1415 / 3.0).sin())) as u32))

    let _ = root.fill(&WHITE);

    let (points, bounds, styles) = prepare_render::<ZZ12>(16, false, Some(2));

    for i in 0..points.len() {
        println!("plotting level {i}");
        let _ = root.fill(&WHITE);
        let root = root.margin(40, 40, 40, 40);

        render(&root, bounds, &points[i..=i], &styles[i..=i], i, 1);

        root.present().unwrap();
    }

    // TODO: move tilezz_examples notebook to examples and update it, test all 3 notebooks
    // TODO: preserve aspect ratio helper (create sub draw area to avoid streching of coorinate system)

    // TODO: 2D point cloud helper (put into plotters util module) - PointSetSequence of points, bounds and shared styles
    // TODO: png/gif creation helper ?
    // TODO: make use of 3D point cloud? how to render?
}
