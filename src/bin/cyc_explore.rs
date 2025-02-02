use std::collections::HashSet;
use std::io::{stdout, Write};

use clap::Parser;

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
fn get_styles(n: usize, pt_sz: i32, pt_sz_const: bool) -> Vec<(i32, ShapeStyle)> {
    // get a rainbow color gradient
    let colors = rainbow(n as u32 + 1, 1.);
    (0..=n)
        .collect::<Vec<_>>()
        .iter()
        .map(|i| {
            (
                ((pt_sz as f64)
                    * (if pt_sz_const {
                        1.
                    } else {
                        0.95_f64.powi(*i as i32)
                    })) as i32,
                colors[*i as usize].filled().into(),
            )
        })
        .collect()
}

fn prepare_render<ZZ: ZZType + HasZZ4>(
    num_rounds: usize,
    mod_unit_square: bool,
    chart_width: u32,
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

    let pt_sz_const = !mod_unit_square;
    let pt_sz = if pt_sz_const {
        2
    } else {
        chart_width as i32 * 3 / 100
    };
    let styles = get_styles(points.len(), pt_sz, pt_sz_const);

    (points, bounds, styles)
}

fn prepare_render_for(
    ring: u8,
    num_rounds: usize,
    mod_unit_square: bool,
    chart_width: u32,
) -> (Vec<Vec<P64>>, R64, Vec<(i32, ShapeStyle)>) {
    match ring {
        4 => prepare_render::<ZZ4>(num_rounds, mod_unit_square, chart_width),
        8 => prepare_render::<ZZ8>(num_rounds, mod_unit_square, chart_width),
        12 => prepare_render::<ZZ12>(num_rounds, mod_unit_square, chart_width),
        16 => prepare_render::<ZZ16>(num_rounds, mod_unit_square, chart_width),
        20 => prepare_render::<ZZ20>(num_rounds, mod_unit_square, chart_width),
        24 => prepare_render::<ZZ24>(num_rounds, mod_unit_square, chart_width),
        60 => prepare_render::<ZZ60>(num_rounds, mod_unit_square, chart_width),
        _ => panic!("invalid ring selected"),
    }
}

// ------------------------------------------------------------------------

#[derive(Clone, Copy)]
enum OutputFormat {
    Png,
    Gif,
}

#[derive(Parser, Debug)]
#[command(version, about = "Explore cyclotomic rings and render the discovered points", long_about = None)]
struct Cli {
    #[arg(short, long)]
    ring: u8,

    #[arg(
        short,
        long,
        help = "Number of BFS exploration rounds (distance from the starting point(s))"
    )]
    num_rounds: usize,

    #[arg(
        short,
        long,
        help = "Run exploration modulo a unit square, starting in its corners"
    )]
    unit_square: bool,

    #[arg(
        short = 'o',
        long,
        help = "Filename (with .png or .gif extension), if missing => dry run"
    )]
    filename: Option<String>,

    #[arg(short, long, default_value_t = 1000, help = "Image width (in px)")]
    width: u32,

    #[arg(short, long, default_value_t = 500, help = "GIF frame delay (in ms)")]
    delay: u32,

    #[arg(short = 'p', long, default_value_t = 4, help = "PNG plots per row")]
    row: usize,

    #[arg(
        short,
        long,
        help = "Number of threads (= # of available cores if unset)"
    )]
    threads: Option<usize>,
}

#[cfg(feature = "examples")]
fn main() {
    let cli = Cli::parse();
    if cli.ring % 4 != 0 {
        panic!("ZZ{} not supported for unit square exploration!", cli.ring);
    }

    let filename = cli.filename.unwrap_or_default();
    let output_format = if filename.is_empty() {
        None
    } else if filename.ends_with(".gif") {
        Some(OutputFormat::Gif)
    } else if filename.ends_with(".png") {
        Some(OutputFormat::Png)
    } else {
        panic!("Unknown image format!")
    };

    println!("Computing points...");
    let chart_width = match output_format.unwrap_or(OutputFormat::Png) {
        OutputFormat::Gif => cli.width,
        OutputFormat::Png => cli.width / (cli.row as u32),
    };
    let (points, bounds, styles) =
        prepare_render_for(cli.ring, cli.num_rounds, cli.unit_square, chart_width);

    if output_format.is_none() {
        return; // dry run, no image
    }
    let output_format = output_format.unwrap();

    let img_dims = (cli.width, cli.width);
    let root = match output_format {
        OutputFormat::Gif => BitMapBackend::gif(&filename, img_dims, cli.delay).unwrap(),
        OutputFormat::Png => BitMapBackend::new(&filename, img_dims),
    }
    .into_drawing_area();

    let grid = match output_format {
        OutputFormat::Gif => (1, 1),
        OutputFormat::Png => (
            cli.row as usize,
            (points.len() / cli.row + points.len() % cli.row) as usize,
        ),
    };
    let areas: Vec<_> = root.split_evenly(grid);

    println!("Plotting points...");
    for i in 0..points.len() {
        let area = &areas[match output_format {
            OutputFormat::Gif => 0,
            OutputFormat::Png => i,
        }];

        // plot into gif frame or png chart matrix grid cell
        let _ = area.fill(&WHITE);
        let area = area.margin(40, 40, 40, 40);
        render(&area, bounds, &points[i..=i], &styles[i..=i], i, 1);

        area.present().unwrap();
    }
}
