use crossbeam;
use num_cpus;
use std::sync::{Arc, Mutex};

use std::sync::atomic::{AtomicUsize, Ordering};

#[cfg(feature = "examples")]
fn main() {
    let values: Vec<i32> = (0..1000).into_iter().collect();
    let num_values = values.len();

    let num_threads = num_cpus::get();
    println!("using {num_threads} threads");
    let per_thread = num_values / num_threads;

    let protected_result: Arc<Mutex<Vec<i32>>> = Arc::new(Mutex::new(Vec::new()));

    // let protected_set: Arc<Mutex<HashSet<i32>>> = Arc::new(Mutex::new(HashSet::new()));

    let thread_id_counter = Arc::new(AtomicUsize::new(0));

    crossbeam::scope(|s| {
        values.chunks(per_thread).for_each(|chunk| {
            let thread_id = thread_id_counter.fetch_add(1, Ordering::SeqCst);
            let result_arc = Arc::clone(&protected_result);

            s.spawn(move |_| {
                for val in chunk.iter() {
                    println!("[{thread_id}] read {val}");
                    result_arc.lock().unwrap().push(*val);
                    println!("[{thread_id}] pushed {val}");
                }
            });
        });
    })
    .unwrap();

    let mutex = Arc::into_inner(protected_result).unwrap();
    let data: Vec<i32> = mutex.into_inner().unwrap();
    for val in data {
        println!("recv {val}");
    }
}
