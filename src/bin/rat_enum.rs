use std::collections::HashSet;
use std::sync::Mutex;
use std::time::Instant;

use clap::{Parser, ValueEnum};

use tilezz::cyclotomic::linalg::norm_sq;
use tilezz::cyclotomic::*;
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::snake::{Snake, Turtle};
use tilezz::vis::animation::render_gif;
use tilezz::vis::draw::{MarkerStyle, TileStyle};
use tilezz::vis::plotutils::P64;
use tilezz::vis::scene::{Color, Fill, Scene, Stroke, TextStyle, Viewport};

static VERBOSE: Mutex<bool> = Mutex::new(false);

// --------

/// Canonical-rotation prune for the open walk `prefix` virtually
/// extended by one more angle `new`. Lets the caller skip walks
/// that cannot be the lex-min cyclic rotation of any closure they
/// could grow into, *before* paying the cyclotomic intersect cost
/// of `Snake::add`.
///
/// Returns `false` when there's a rotation index `k > 0` such that
/// the rotation `[prefix++new][k..]` is strictly lex-less than
/// `[prefix++new][0..]` within the wrap-free comparable region; that
/// decision is permanent (no future angles can flip it) so the
/// eventual closed polygon's canonical rotation cannot start at
/// position 0 -- pruning is safe.
fn is_canonical_extended(prefix: &[i8], new: i8) -> bool {
    let base = prefix.len();
    let d = base + 1;
    let get = |i: usize| -> i8 {
        if i < base { prefix[i] } else { new }
    };
    for k in 1..d {
        for i in 0..(d - k) {
            match get(k + i).cmp(&get(i)) {
                std::cmp::Ordering::Less => return false,
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => continue,
            }
        }
    }
    true
}

/// Enumerate every simple polygon with boundary length up to
/// `max_steps` over the cyclotomic ring `ZZ`, in canonical-CCW form.
///
/// Depth-first walk: a single `Snake` is mutated in place via
/// [`Snake::add`] / [`Snake::pop`] across the recursion, so the inner
/// loop allocates only when a polygon actually closes (the canonical
/// [`Rat`] needs a fresh vec for hashing). Each polygon's `n` cyclic
/// walks collapse to a single canonical walk via the lex-min
/// rotation prune in [`is_canonical_extended`].
pub fn rat_enum<ZZ: ZZType + Units>(max_steps: usize) -> Vec<Vec<i8>> {
    let mut result: HashSet<Vec<i8>> = HashSet::new();
    let mut snake: Snake<ZZ> = Snake::new();

    println!("-------- enumeration started --------");
    rat_enum_step::<ZZ>(&mut snake, max_steps, &mut result);
    println!(
        "-------- enumeration completed --------\n{} rats found",
        result.len()
    );

    let mut result: Vec<Vec<i8>> = result.into_iter().collect();
    result.sort_by_key(|x| x.len());
    result
}

fn rat_enum_step<ZZ: ZZType + Units>(
    snake: &mut Snake<ZZ>,
    max_steps: usize,
    result: &mut HashSet<Vec<i8>>,
) {
    let depth = snake.angles().len();
    if depth >= max_steps {
        return;
    }
    // Heuristic: the snake's current head must lie within `remaining`
    // unit steps of the origin, else the path can never close in time.
    let remaining = (max_steps - depth) as i64;
    let remaining_steps_sq = norm_sq(&ZZ::from(remaining));

    for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
        // Canonical-rotation prune (see `is_canonical_extended`):
        // pre-check *before* calling `Snake::add` so that
        // canonical-rejected branches don't pay the cyclotomic
        // intersect cost.
        if !is_canonical_extended(snake.angles(), direction) {
            continue;
        }
        if !snake.add(direction) {
            continue;
        }
        if snake.is_closed() {
            // polygon completed -> record canonical ccw description.
            let r = {
                let tmp = Rat::from_unchecked(snake);
                if tmp.chirality() > 0 {
                    tmp
                } else {
                    tmp.reversed()
                }
                .canonical()
            };
            let seq = r.seq().to_vec();
            if result.insert(seq.clone()) {
                println!("RAT {seq:?}");
            }
        } else if norm_sq(&snake.head().pos) <= remaining_steps_sq {
            // Still reachable -- recurse to extend further.
            rat_enum_step::<ZZ>(snake, max_steps, result);
        }
        // Backtrack: restore the snake to the pre-add state.
        snake.pop();
    }
}

fn polygons<ZZ: ZZType + Units>(rats: Vec<Vec<i8>>) -> Vec<Vec<P64>> {
    rats.into_iter()
        .map(|seq| Rat::<ZZ>::from_slice_unchecked(&seq).to_polyline_f64(Turtle::default()))
        .collect()
}

fn run_rat_enum_polylines(ring: u8, max_steps: usize) -> Vec<Vec<P64>> {
    match ring {
        4 => polygons::<ZZ4>(rat_enum::<ZZ4>(max_steps)),
        8 => polygons::<ZZ8>(rat_enum::<ZZ8>(max_steps)),
        12 => polygons::<ZZ12>(rat_enum::<ZZ12>(max_steps)),
        16 => polygons::<ZZ16>(rat_enum::<ZZ16>(max_steps)),
        20 => polygons::<ZZ20>(rat_enum::<ZZ20>(max_steps)),
        24 => polygons::<ZZ24>(rat_enum::<ZZ24>(max_steps)),
        60 => polygons::<ZZ60>(rat_enum::<ZZ60>(max_steps)),
        _ => panic!("invalid ring selected"),
    }
}

fn run_rat_enum_seqs(ring: u8, max_steps: usize) -> Vec<Vec<i8>> {
    let f: fn(usize) -> Vec<Vec<i8>> = match ring {
        4 => rat_enum::<ZZ4>,
        8 => rat_enum::<ZZ8>,
        12 => rat_enum::<ZZ12>,
        16 => rat_enum::<ZZ16>,
        20 => rat_enum::<ZZ20>,
        24 => rat_enum::<ZZ24>,
        60 => rat_enum::<ZZ60>,
        _ => panic!("invalid ring selected"),
    };
    f(max_steps)
}

// --------

#[derive(Copy, Clone, Debug, ValueEnum)]
enum Mode {
    /// Enumerate and render output (GIF)
    Render,
    /// Enumerate only and report elapsed time
    Bench,
}

#[derive(Parser, Debug)]
#[command(version, about = "Compute all simple polygons over a ring with a maximal boundary length", long_about = None)]
struct Cli {
    #[arg(short = 'r', long)]
    ring: u8,

    #[arg(short = 'n', long)]
    max_steps: usize,

    #[arg(long, value_enum, default_value_t = Mode::Render)]
    mode: Mode,

    #[arg(short = 'o', long, help = "Output GIF filename (render mode only)")]
    filename: Option<String>,

    #[arg(short, long)]
    verbose: bool,
}

fn main() {
    let cli = Cli::parse();
    if cli.verbose {
        let mut verbose = VERBOSE.lock().unwrap();
        *verbose = true;
    }

    match cli.mode {
        Mode::Bench => {
            let t0 = Instant::now();
            let rats: Vec<Vec<i8>> = run_rat_enum_seqs(cli.ring, cli.max_steps);
            let dt = t0.elapsed();

            // Use the result so it can't be trivially optimized away.
            let total_boundary_len: usize = rats.iter().map(|s| s.len()).sum();

            println!(
                "benchmark: ring={} max_steps={} -> {} unique rats (total boundary len={}) in {:?}",
                cli.ring,
                cli.max_steps,
                rats.len(),
                total_boundary_len,
                dt
            );
        }
        Mode::Render => {
            let rats: Vec<Vec<P64>> = run_rat_enum_polylines(cli.ring, cli.max_steps);

            let Some(filename) = cli.filename else {
                return;
            };

            // Symmetric square viewport sized to fit every polygon
            // centered at the origin. Each polygon goes through
            // Rat::to_polyline_f64 with the default Turtle (origin,
            // direction 0), so we just need the max |x|/|y| across
            // all polylines.
            let mut max_abs: f64 = 1.0;
            for poly in &rats {
                for (x, y) in poly {
                    max_abs = max_abs.max(x.abs()).max(y.abs());
                }
            }
            let pad = 0.15 * max_abs;
            let r = max_abs + pad;
            let bounds = ((-r, -r), (r, r));

            // Style: filled yellow tile with thin black border,
            // red vertex markers labelled with their index.
            let style = TileStyle::filled(
                Fill::solid(Color::YELLOW.with_alpha(80)),
                Stroke::solid(Color::BLACK, 0.01 * r),
            )
            .with_vertex_marker(MarkerStyle::filled_circle(
                0.06 * r,
                Color::RED,
            ))
            .with_vertex_labels(TextStyle::new(0.04 * r, Color::WHITE).bold());

            let frames: Vec<Scene> = rats
                .iter()
                .enumerate()
                .map(|(ix, tile)| {
                    println!("rendering frame {}/{}", ix + 1, rats.len());
                    let mut scene = Scene::new().with_background(Color::WHITE);
                    scene.draw_tile(tile, &style);
                    scene
                })
                .collect();

            let w = 500u32;
            let vp = Viewport::square_for(w, bounds, 8);
            let gif_bytes = render_gif(&frames, &vp, 500).expect("render GIF");
            std::fs::write(&filename, gif_bytes).expect("write GIF");
            println!("wrote {filename}");
        }
    }
}
