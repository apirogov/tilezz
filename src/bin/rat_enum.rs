use std::sync::Mutex;

use clap::Parser;

use std::collections::BTreeSet;

use itertools::Itertools;
use plotters::prelude::*;

use tilezz::cyclotomic::linalg::norm_sq;
use tilezz::cyclotomic::*;
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::snake::{Snake, Turtle};
use tilezz::vis::plotutils::P64;

static VERBOSE: Mutex<bool> = Mutex::new(false);

// --------

// this seems to be slower... maybe should optimize by using a trie of snakes
// if one shorter snake does not work, no point using one with same prefix!
pub fn rat_enum_alt<ZZ: ZZType>(max_steps: usize) -> Vec<Vec<i8>> {
    let mut result: BTreeSet<Vec<i8>> = BTreeSet::new();

    // simple snakes of a fixed length (initialized with length 0 (dummy) and 1 (start))
    let mut k_snakes: Vec<Vec<Vec<i8>>> = vec![vec![vec![]]];
    k_snakes.push(((-ZZ::hturn() + 1)..ZZ::hturn()).map(|i| vec![i]).collect());
    let mut total_snakes: usize = (ZZ::turn() - 1) as usize;

    println!("-------- enumeration started --------");
    let status_log = |step, num_snakes: usize, result: &BTreeSet<_>| {
        println!(
            "round {}: {} snakes, {} rats",
            step,
            num_snakes,
            result.len()
        );
    };

    status_log(1, total_snakes, &result);
    for len in 2..=max_steps {
        let remaining = (max_steps - len) as i64;
        let remaining_steps_sq = norm_sq(&ZZ::from(remaining));

        let mut next: Vec<Vec<i8>> = Vec::new();

        // for snakes of length n, we use length x=floor(n/2) and y=n-x
        // that way we have confirmed simple snakes to combine,
        // instead of naively recomputing all possibilities stepwise
        let i1 = len / 2;
        let i2 = len - i1;
        for (s1, s2) in k_snakes[i1].iter().cartesian_product(k_snakes[i2].iter()) {
            // load first half without checks, then trace out second half with checks
            let mut s = Snake::<ZZ>::from_slice_unsafe(s1);
            let r = s.extend_from_slice(s2);
            if r.is_err() {
                continue; // self-crossing
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
                let seq = r.seq().to_vec();
                if result.insert(seq.clone()) {
                    total_snakes += 1;
                    println!("RAT {seq:?}");
                }
            }

            // heuristic - path must not be too far from origin,
            // i.e. the way back to origin must be doable
            // within the remaining number of steps.
            if norm_sq(&s.head().pos) > remaining_steps_sq {
                continue; // no way to close path
            }

            // add the prolonged path

            next.push(s.angles().to_vec());
            println!("SNAKE {:?}", s.angles());
        }

        k_snakes.push(next);
        status_log(len, total_snakes, &result);
    }
    println!("-------- enumeration completed --------");

    let mut result: Vec<Vec<i8>> = result.into_iter().collect();
    result.sort_by(|x, y| x.len().cmp(&y.len()));
    result
}

pub fn rat_enum<ZZ: ZZType>(max_steps: usize) -> Vec<Vec<i8>> {
    let mut result: BTreeSet<Vec<i8>> = BTreeSet::new();

    let mut snakes: Vec<Vec<i8>> = vec![vec![]];

    println!("-------- enumeration started --------");
    let status_log = |step, snakes: &Vec<_>, result: &BTreeSet<_>| {
        println!(
            "round {}: {} alive snakes, {} unique rats",
            step,
            snakes.len(),
            result.len()
        );
    };

    for step in 0..max_steps {
        status_log(step, &snakes, &result);

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
                    let seq = r.seq().to_vec();
                    if result.insert(seq.clone()) {
                        println!("RAT {seq:?}");
                    }
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
    status_log(max_steps, &snakes, &result);
    println!("-------- enumeration completed --------");

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

    let w = 500;
    // let root = BitMapBackend::new(&filename, (w, w)).into_drawing_area();
    let root = BitMapBackend::gif(filename, (w, w), 500)
        .unwrap()
        .into_drawing_area();

    let _ = root.fill(&WHITE);

    // GIF
    let grid = (1, 1);
    // PNG
    // let row = 3;
    // let grid = ((rats.len() / row + rats.len() % row) as usize, row as usize);

    let areas: Vec<_> = root.split_evenly(grid);

    let mut style = TileStyle::default();
    style.node_size = 0;
    style.node_labels = false;
    for (ix, tile) in rats.iter().enumerate() {
        println!("plotting tile {}/{}", ix + 1, rats.len());

        // let area = &areas[ix]; // PNG
        let area = &areas[0]; // GIF

        let _ = area.fill(&WHITE);

        let pad = w * 2 / 100;
        let mut area = area.margin(pad, pad, pad, pad);

        plot_tile(&mut area, &tile, &style);

        area.present().unwrap();
    }
}
