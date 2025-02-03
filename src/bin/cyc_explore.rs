use std::collections::HashSet;
use std::io::{stdout, Write};
use std::marker::{Send, Sync};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use clap::Parser;

use crossbeam;
use plotters::coord::Shift;
use plotters::prelude::*;

use tilezz::cyclotomic::geometry::point_mod_rect;
use tilezz::cyclotomic::*;
use tilezz::vis::plotters::{plot_points, rainbow};
use tilezz::vis::plotutils::{chart_padding, tiles_bounds, P64, R64};

static VERBOSE: Mutex<bool> = Mutex::new(false);

/// Compute the levels of points reachable within the unit square from any Gaussian integer in n steps.
fn explore<ZZ: ZZType + HasZZ4 + Send + Sync>(
    n: usize,
    mod_unit_square: bool,
    num_threads: usize,
) -> Vec<Vec<ZZ>>
where
    <ZZ as IsComplex>::Field: From<(<ZZ as IsRingOrField>::Real, <ZZ as IsRingOrField>::Real)>,
{
    // we start at the corners of the unit square
    let start_pts: &[ZZ] = if mod_unit_square {
        &[ZZ::zero(), ZZ::one(), ZZ::one_i(), ZZ::one() + ZZ::one_i()]
    } else {
        &[ZZ::zero()]
    };

    // FIXME: as distributing chunks alone does not make it as fast as expected,
    // it seems like access to the visited map is the bottleneck!
    let protected_visited: Arc<Mutex<HashSet<ZZ>>> = {
        let visited: HashSet<ZZ> = HashSet::from_iter(start_pts.to_vec().into_iter());
        Arc::new(Mutex::new(visited))
    };
    let mut round_pts: Vec<Vec<ZZ>> = Vec::new();
    round_pts.push(start_pts.to_vec());

    let unit_square: (ZZ, ZZ) = ((0, 0).into(), (1, 1).into());
    // in each round, we go one unit step in every possible direction
    for i in 1..=n {
        let last = round_pts.last().unwrap();
        let protected_curr: Arc<Mutex<Vec<ZZ>>> = Arc::new(Mutex::new(Vec::new()));

        let per_thread = num_threads.max(last.len() / num_threads);

        let thread_id_counter = Arc::new(AtomicUsize::new(0));
        crossbeam::scope(|s| {
            last.as_slice().chunks(per_thread).for_each(|chunk| {
                let _ = thread_id_counter.fetch_add(1, Ordering::SeqCst);
                let curr_arc = Arc::clone(&protected_curr);
                let visited_arc = Arc::clone(&protected_visited);

                s.spawn(move |_| {
                    for p in chunk.iter() {
                        for i in 0..ZZ::turn() {
                            let mut p_dest: ZZ = *p + ZZ::unit(i);
                            if mod_unit_square {
                                // normalize back into unit square (if enabled)
                                p_dest = point_mod_rect(&p_dest, &unit_square).coerce_ring();
                            }
                            let p_dest = p_dest;

                            if !visited_arc.lock().unwrap().contains(&p_dest) {
                                visited_arc.lock().unwrap().insert(p_dest);
                                curr_arc.lock().unwrap().push(p_dest);
                            }
                        }
                    }
                });
            });
        })
        .unwrap();

        let used_threads = Arc::into_inner(thread_id_counter).unwrap().into_inner();

        let mutex = Arc::into_inner(protected_curr).unwrap();
        let curr: Vec<ZZ> = mutex.into_inner().unwrap();

        print!("{}{}", curr.len(), if i == n { "" } else { " + " }); // print number of new points
        stdout().flush().unwrap();
        if *VERBOSE.lock().unwrap() {
            println!("\n(used {used_threads} threads, each on around {per_thread} numbers)");
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
    bounds @ ((x_min, y_min), (x_max, y_max)): R64,
    point_levels: &[Vec<P64>],
    level_styles: &[(i32, ShapeStyle)],
    offset: usize,
    stride: u32,
) {
    let (pad_x, pad_y) = chart_padding(da.dim_in_pixel(), bounds);
    println!("{pad_x} {pad_y}");
    let da = da.margin(pad_y / 2, pad_y / 2, pad_x / 2, pad_x / 2);

    // prepare coordinate system
    let mut chart = ChartBuilder::on(&da)
        .x_label_area_size(20)
        .y_label_area_size(20)
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
                        0.98_f64.powi(*i as i32)
                    })) as i32,
                colors[*i as usize].filled().into(),
            )
        })
        .collect()
}

fn prepare_render<ZZ: ZZType + HasZZ4 + Send + Sync>(
    num_rounds: usize,
    mod_unit_square: bool,
    chart_width: u32,
    num_threads: usize,
) -> (Vec<Vec<P64>>, R64, Vec<(i32, ShapeStyle)>)
where
    <ZZ as IsComplex>::Field: From<(<ZZ as IsRingOrField>::Real, <ZZ as IsRingOrField>::Real)>,
{
    let points: Vec<Vec<P64>> = explore::<ZZ>(num_rounds, mod_unit_square, num_threads)
        .iter()
        .map(|v| v.iter().map(|p| p.xy()).collect())
        .collect();

    let bounds = tiles_bounds(&points);

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
    num_threads: usize,
) -> (Vec<Vec<P64>>, R64, Vec<(i32, ShapeStyle)>) {
    match ring {
        4 => prepare_render::<ZZ4>(num_rounds, mod_unit_square, chart_width, num_threads),
        8 => prepare_render::<ZZ8>(num_rounds, mod_unit_square, chart_width, num_threads),
        12 => prepare_render::<ZZ12>(num_rounds, mod_unit_square, chart_width, num_threads),
        16 => prepare_render::<ZZ16>(num_rounds, mod_unit_square, chart_width, num_threads),
        20 => prepare_render::<ZZ20>(num_rounds, mod_unit_square, chart_width, num_threads),
        24 => prepare_render::<ZZ24>(num_rounds, mod_unit_square, chart_width, num_threads),
        60 => prepare_render::<ZZ60>(num_rounds, mod_unit_square, chart_width, num_threads),
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

    #[arg(short, long)]
    verbose: bool,
}

#[cfg(feature = "examples")]
fn main() {
    let cli = Cli::parse();
    if cli.verbose {
        let mut verbose = VERBOSE.lock().unwrap();
        *verbose = true;
    }
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

    let img_dims = (cli.width, cli.width);
    let chart_width = match output_format.unwrap_or(OutputFormat::Png) {
        OutputFormat::Gif => cli.width,
        OutputFormat::Png => cli.width / (cli.row as u32),
    };

    // -------- Compute --------

    let num_threads = cli.threads.unwrap_or(num_cpus::get());
    if *VERBOSE.lock().unwrap() {
        println!("Computing points using {num_threads} threads...");
    }

    let (points, bounds, styles) = prepare_render_for(
        cli.ring,
        cli.num_rounds,
        cli.unit_square,
        chart_width,
        num_threads,
    );

    // -------- Render --------

    if output_format.is_none() {
        return; // dry run -> computation with no rendering
    }
    let output_format = output_format.unwrap();

    // get image of desired size for chosen backend depending of format
    let root = match output_format {
        OutputFormat::Gif => BitMapBackend::gif(&filename, img_dims, cli.delay).unwrap(),
        OutputFormat::Png => BitMapBackend::new(&filename, img_dims),
    }
    .into_drawing_area();

    if *VERBOSE.lock().unwrap() {
        println!("Plotting points...");
    }

    // frame chart of GIF <-> cell of a grid of plots in PNG
    let grid = match output_format {
        OutputFormat::Gif => (1, 1),
        OutputFormat::Png => (
            (points.len() / cli.row + points.len() % cli.row) as usize,
            cli.row as usize,
        ),
    };
    let areas: Vec<_> = root.split_evenly(grid);

    for i in 0..points.len() {
        let area = &areas[match output_format {
            OutputFormat::Gif => 0,
            OutputFormat::Png => i,
        }];

        let _ = area.fill(&WHITE);

        let pad = cli.width * 2 / 100;
        let area = area.margin(pad, pad, pad, pad);

        // plot into gif frame or png chart matrix grid cell
        render(&area, bounds, &points[i..=i], &styles[i..=i], i, 1);
        area.present().unwrap();
    }
}
