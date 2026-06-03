use std::collections::HashSet;
use std::io::{Write, stdout};
use std::marker::{Send, Sync};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use clap::Parser;

use tilezz::cyclotomic::Units;
use tilezz::cyclotomic::geometry::point_mod_rect;
use tilezz::cyclotomic::*;
use tilezz::vis::animation::render_gif;
use tilezz::vis::draw::rainbow;
use tilezz::vis::plotutils::{P64, R64, points_bounds};
use tilezz::vis::scene::{Color, Fill, Item, MarkerShape, Scene, Stroke, Viewport};

static VERBOSE: Mutex<bool> = Mutex::new(false);

/// Compute the levels of points reachable within the unit square from any Gaussian integer in n steps.
fn explore<ZZ>(n: usize, mod_unit_square: bool, num_threads: usize) -> Vec<Vec<ZZ>>
where
    ZZ: HasZZ4 + Send + Sync,
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
        let visited: HashSet<ZZ> = HashSet::from_iter(start_pts.to_vec());
        Arc::new(Mutex::new(visited))
    };
    let mut round_pts: Vec<Vec<ZZ>> = Vec::new();
    round_pts.push(start_pts.to_vec());

    let anchor: ZZ = ZZ::zero();
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
                            let mut p_dest: ZZ = *p + <ZZ as Units>::unit(i);
                            if mod_unit_square {
                                // normalize back into unit square (if enabled)
                                p_dest = point_mod_rect(&p_dest, &anchor, (1, 1));
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
    round_pts
}

fn prepare_render<ZZ>(
    num_rounds: usize,
    mod_unit_square: bool,
    num_threads: usize,
) -> (Vec<Vec<P64>>, R64)
where
    ZZ: HasZZ4 + Send + Sync,
{
    let points: Vec<Vec<P64>> = explore::<ZZ>(num_rounds, mod_unit_square, num_threads)
        .iter()
        .map(|v| v.iter().map(|p| p.xy()).collect())
        .collect();
    let bounds = points_bounds(points.iter()).unwrap_or(((-0.5, -0.5), (0.5, 0.5)));
    (points, bounds)
}

fn prepare_render_for(
    ring: u8,
    num_rounds: usize,
    mod_unit_square: bool,
    num_threads: usize,
) -> (Vec<Vec<P64>>, R64) {
    match ring {
        4 => prepare_render::<ZZ4>(num_rounds, mod_unit_square, num_threads),
        8 => prepare_render::<ZZ8>(num_rounds, mod_unit_square, num_threads),
        12 => prepare_render::<ZZ12>(num_rounds, mod_unit_square, num_threads),
        16 => prepare_render::<ZZ16>(num_rounds, mod_unit_square, num_threads),
        20 => prepare_render::<ZZ20>(num_rounds, mod_unit_square, num_threads),
        24 => prepare_render::<ZZ24>(num_rounds, mod_unit_square, num_threads),
        60 => prepare_render::<ZZ60>(num_rounds, mod_unit_square, num_threads),
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

#[cfg(feature = "cli")]
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

    // -------- Compute --------

    let num_threads = cli.threads.unwrap_or(num_cpus::get());
    if *VERBOSE.lock().unwrap() {
        println!("Computing points using {num_threads} threads...");
    }

    let (points, bounds) =
        prepare_render_for(cli.ring, cli.num_rounds, cli.unit_square, num_threads);

    // -------- Render --------

    let Some(output_format) = output_format else {
        return; // dry run -> computation with no rendering
    };

    render_vis(
        &filename,
        output_format,
        img_dims,
        &points,
        bounds,
        cli.num_rounds,
        cli.row,
        cli.delay,
        cli.unit_square,
    );
}

// ------------------------------------------------------------------------
// Scene-graph rendering (vis::scene / vis::raster / vis::animation).
// ------------------------------------------------------------------------

/// Build a single Scene containing the round-`level` point cloud,
/// shifted by `(dx, dy)` math units and rendered with the given
/// `marker_size` (scene units) and color.
fn add_cell(
    scene: &mut Scene,
    points: &[P64],
    bounds: R64,
    (dx, dy): (f64, f64),
    marker_size: f64,
    color: Color,
) {
    let ((mn_x, mn_y), (mx_x, mx_y)) = bounds;
    // Light cell border.
    scene.push(Item::Polygon {
        points: vec![
            (dx + mn_x, dy + mn_y),
            (dx + mx_x, dy + mn_y),
            (dx + mx_x, dy + mx_y),
            (dx + mn_x, dy + mx_y),
        ],
        fill: None,
        stroke: Some(Stroke::solid(
            Color::rgb(180, 180, 180),
            0.005 * (mx_x - mn_x),
        )),
        arrow: None,
    });
    // Cross-hairs at origin: helps see scale + symmetry.
    scene.push(Item::Segment {
        a: (dx + mn_x, dy),
        b: (dx + mx_x, dy),
        stroke: Stroke::solid(Color::rgb(225, 225, 225), 0.003 * (mx_x - mn_x)),
        arrow: None,
    });
    scene.push(Item::Segment {
        a: (dx, dy + mn_y),
        b: (dx, dy + mx_y),
        stroke: Stroke::solid(Color::rgb(225, 225, 225), 0.003 * (mx_x - mn_x)),
        arrow: None,
    });
    // Points. Thin black outline so light-coloured palette entries
    // (yellow, cyan, …) stand out against the white background;
    // ratio 0.10 keeps the outline visible without swallowing the
    // fill colour at 6-pixel marker diameters.
    let outline_width = marker_size * 0.10;
    for &(x, y) in points {
        scene.push(Item::Marker {
            center: (x + dx, y + dy),
            shape: MarkerShape::Circle,
            size: marker_size,
            fill: Some(Fill::solid(color)),
            stroke: Some(Stroke::solid(Color::BLACK, outline_width)),
        });
    }
}

/// Render via the new scene-graph backend. `output_format` decides
/// PNG-grid vs animated GIF.
#[allow(clippy::too_many_arguments)]
fn render_vis(
    filename: &str,
    output_format: OutputFormat,
    img_dims: (u32, u32),
    points: &[Vec<P64>],
    bounds: R64,
    num_rounds: usize,
    cols: usize,
    delay_ms: u32,
    unit_square: bool,
) {
    let n = points.len();
    let palette: Vec<Color> = rainbow(n.max(1), 1.0, 0.5);
    let ((mn_x, mn_y), (mx_x, mx_y)) = bounds;
    let cell_w = mx_x - mn_x;
    let cell_h = mx_y - mn_y;

    // Compute a marker size that renders to roughly 6 pixels of
    // diameter at the actual per-cell pixel resolution of each
    // output format. For the PNG grid the cell is `img_w / cols`
    // pixels wide; for the GIF each frame fills the full image.
    let marker_size_for = |cell_px_w: f64| -> f64 {
        if unit_square {
            0.04 * cell_w
        } else {
            let pixel_per_unit = cell_px_w / cell_w;
            6.0 / pixel_per_unit
        }
    };

    match output_format {
        OutputFormat::Png => {
            let rows = (num_rounds + 1).div_ceil(cols); // levels 0..=num_rounds, ceil-divided
            let gap = 0.05 * cell_w;
            let cell_px_w = (img_dims.0 as f64) / cols as f64;
            let marker_size = marker_size_for(cell_px_w);
            let mut scene = Scene::new().with_background(Color::WHITE);
            for (i, pts) in points.iter().enumerate() {
                let col = i % cols;
                let row = i / cols;
                // Place cells with row 0 at the TOP (highest math y).
                let dx = col as f64 * (cell_w + gap) - mn_x;
                let dy = (rows - 1 - row) as f64 * (cell_h + gap) - mn_y;
                add_cell(&mut scene, pts, bounds, (dx, dy), marker_size, palette[i]);
            }
            let total_w = cols as f64 * cell_w + (cols.saturating_sub(1)) as f64 * gap;
            let total_h = rows as f64 * cell_h + (rows.saturating_sub(1)) as f64 * gap;
            let total_w_px = img_dims.0;
            let total_h_px = ((total_h / total_w) * total_w_px as f64).round().max(1.0) as u32;
            let vp =
                Viewport::rect_for(total_w_px, total_h_px, ((0.0, 0.0), (total_w, total_h)), 16);
            let png = scene.to_png(&vp).expect("render PNG");
            std::fs::write(filename, png).expect("write PNG");
        }
        OutputFormat::Gif => {
            // Each frame fills the whole image, so the effective
            // cell-pixel-width is the full image width.
            let marker_size = marker_size_for(img_dims.0 as f64);
            let mut frames: Vec<Scene> = Vec::with_capacity(n);
            for (i, pts) in points.iter().enumerate() {
                let mut frame = Scene::new().with_background(Color::WHITE);
                add_cell(
                    &mut frame,
                    pts,
                    bounds,
                    (-mn_x, -mn_y),
                    marker_size,
                    palette[i],
                );
                frames.push(frame);
            }
            let vp = Viewport::square_for(img_dims.0, ((0.0, 0.0), (cell_w, cell_h)), 16);
            let gif_bytes = render_gif(&frames, &vp, delay_ms as u16).expect("render GIF");
            std::fs::write(filename, gif_bytes).expect("write GIF");
        }
    }
}
