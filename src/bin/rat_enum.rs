use std::sync::Mutex;

use clap::Parser;

use std::collections::BTreeSet;

use plotters::prelude::*;

use tilezz::cyclotomic::linalg::norm_sq;
use tilezz::cyclotomic::*;
use tilezz::rat::Rat;
use tilezz::snake::{Snake, Turtle};
use tilezz::vis::plotutils::P64;

static VERBOSE: Mutex<bool> = Mutex::new(false);

// --------

pub fn rat_enum<ZZ: ZZType>(max_steps: usize) -> Vec<Vec<i8>> {
    let mut result: BTreeSet<Vec<i8>> = BTreeSet::new();
    let mut snakes: Vec<Vec<i8>> = vec![vec![]];

    for step in 0..max_steps {
        println!(
            "{}: {} alive snakes, {} unique rats",
            step,
            snakes.len(),
            result.len()
        );

        let remaining = (max_steps - step) as i64;
        let remaining_steps_sq = norm_sq(&ZZ::from(remaining));

        let mut next: Vec<Vec<i8>> = Vec::new();
        for seq in snakes.iter() {
            for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
                let mut s: Snake<ZZ> = Snake::try_from(seq.as_slice()).unwrap();
                if !s.add(direction) {
                    continue; // invalid step (self-crossing)
                }
                if s.is_closed() {
                    // polygon completed -> get canonical ccw description
                    let r = {
                        let tmp = Rat::from_unchecked(&s);
                        if tmp.chirality() > 0 {
                            tmp
                        } else {
                            tmp.reversed()
                        }
                        .canonical()
                    };
                    result.insert(r.seq().to_vec());
                }

                // heuristic - path must not be too far from origin,
                // i.e. the way back to origin must be doable
                // within the remaining number of steps.
                if norm_sq(&s.head().pos) > remaining_steps_sq {
                    continue; // no way to close path
                }

                // add the prolonged path
                next.push(s.angles().to_vec());
            }
        }
        snakes = next;
    }

    let mut result: Vec<Vec<i8>> = result.into_iter().collect();
    result.sort_by(|x, y| x.len().cmp(&y.len()));
    result
}

fn polygons<ZZ: ZZType>(rats: Vec<Vec<i8>>) -> Vec<Vec<P64>> {
    rats.into_iter()
        .map(|seq| Rat::<ZZ>::from_slice_unchecked(&seq).to_polyline_f64(Turtle::default()))
        .collect()
}

fn run_rat_enum(ring: u8, max_steps: usize) -> Vec<Vec<P64>> {
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

// --------

#[derive(Parser, Debug)]
#[command(version, about = "Compute all simple polygons over a ring with a maximal boundary length", long_about = None)]
struct Cli {
    #[arg(short = 'r', long)]
    ring: u8,

    #[arg(short = 'n', long)]
    max_steps: usize,

    #[arg(short = 'o', long, help = "Output GIF filename")]
    filename: Option<String>,

    #[arg(short, long)]
    verbose: bool,
}

#[cfg(feature = "examples")]
fn main() {
    use tilezz::vis::plotters::{plot_tile, TileStyle};

    let cli = Cli::parse();
    if cli.verbose {
        let mut verbose = VERBOSE.lock().unwrap();
        *verbose = true;
    }

    let rats: Vec<Vec<P64>> = run_rat_enum(cli.ring, cli.max_steps);

    if cli.filename.is_none() {
        return;
    }
    let filename = cli.filename.unwrap();

    let mut root = BitMapBackend::gif(filename, (500, 500), 500)
        .unwrap()
        .into_drawing_area();

    for (ix, tile) in rats.iter().enumerate() {
        println!("plotting tile {}/{}", ix + 1, rats.len());

        let _ = root.fill(&WHITE);
        plot_tile(&mut root, &tile, &TileStyle::default());
        root.present().unwrap();
    }
}
