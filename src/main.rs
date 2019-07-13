extern crate rand;
extern crate num_cpus;

pub mod system;
use system::System;

use rand::Rng;
use std::fs::File;
use std::io::Write;

use std::collections::VecDeque;

use std::sync::{Arc, Mutex};
use std::sync::mpsc::channel;
use std::thread;

/// Linear dimension of the Lattice.
pub const N: i32 = 64;

/// Lowest temperature to consider.
pub const T0: f64 = 0.1;
/// Highest temperature to consider.
pub const TF: f64 = 5.;
/// Number of temperatures.
pub const STEPS: usize = 400;

fn push_neighbors(sys: &System,
                  i: usize, j: usize,
                  queue: &mut VecDeque<(usize, usize)>) {

    queue.push_back((i, sys.nnplus(j)));
    queue.push_back((i, sys.nnminus(j)));
    queue.push_back((sys.nnplus(i), j));
    queue.push_back((sys.nnminus(i), j));
}

#[allow(non_snake_case)]
/// Perform one step of the Metropolis algorithm
///
/// We pick a random site on the lattice and propose
/// flipping it.
/// This change is accepted with probability `exp(-dE/T)`
/// where `dE` is the change in the energy due to the flip
/// and `T` is the temperature of the system.
pub fn metropolis_step(sys: &mut System,
                   beta: f64) -> (i32, i32) {
    let i = sys.random_site();
    let j = sys.random_site();

    let mut neighbors: VecDeque<(usize, usize)> = VecDeque::new();
    let mut dE = sys.deltaE(i, j);
    let s = sys.flip(i, j);
    let mut dM = -2*s;
    let flip_prob = 1. - (-2.*beta).exp();

    push_neighbors(sys, i, j, &mut neighbors);

    while let Some((i, j)) = neighbors.pop_front() {

        if sys.at(i, j) == s && sys.rng.next_f64() < flip_prob {
            dM -= 2*s;
            dE += sys.deltaE(i, j);

            sys.flip(i, j);

            push_neighbors(sys, i, j, &mut neighbors);
        }

    }

    (dE, dM)
}

/// Evolve the system for `n` steps.
///
/// We propose `n` changes to the system which are
/// accepted or rejected according to their Boltzmann factor.
/// Observables are not tracked.
pub fn evolve(sys: &mut System,
          n: i32,
          beta: f64) {
    for _ in 0..n { metropolis_step(sys, beta); }
}

#[allow(non_snake_case)]
/// Accumulate the energy and magenization of the system over `n` steps.
///
/// We perform `n` [Metropolis steps](metropolis_step), noting the energy and
/// magnetization at each step. We then return the average of these quantities.
pub fn time_average(sys: &mut System,
                n: i32,
                beta: f64) -> (f64, f64) {

    let mut U = sys.evaluate_energy();
    let mut M = sys.magnetization();


    let mut M_tot: i64 = 0;
    let mut U_tot: i64 = 0;

    for _ in 0..n {
        let (dE, dM) = metropolis_step(sys, beta);

        U += dE;
        M += dM;

        U_tot += U as i64;
        M_tot += (M*M) as i64;
    }

    let en = (U_tot as f64) / (n as f64);
    let mag = ((M_tot as f64) / (n as f64)).sqrt();

    (en, mag)
}

/// Perform a Monte carlo average of magnetization and energy after a burn in phase.
///
/// We burn in the Markov chain by first [evolve]()ing `n_evolve` steps, then average the
/// magnetization over `n_average` steps.
/// We employ the [Metropolis-Hastings
/// algorithm](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm) to compute the
/// averages.
pub fn ensemble_average(sys: &mut System,
                    beta: f64,
                    n_evolve: i32,
                    n_average: i32) -> (f64, f64) {
    evolve(sys, n_evolve, beta);
    time_average(sys, n_average, beta)
}

/// Fetch the next temperature to be calculated
fn next_t<T>(ts: &Arc<Mutex<VecDeque<T>>>)
    -> Option<T> {
        ts.lock().unwrap().pop_front()
    }

#[allow(non_snake_case)]
fn main() {

    /// Number of threads to parallelize computation over.
    let NUM_THREADS: usize = num_cpus::get();


    // Note the -1. This ensures that we actually reach TF
    let dt = (TF - T0)/(STEPS as f64 - 1.);

    // This is a Deque so that threads can grab the next
    // temperature as needed.
    //
    // We want to go from low temperature to high, so using
    // a Deque means we don't have to reverse `ts` for
    // proper pop semantics.
    let ts = (0..STEPS)
        .map(|i| T0 + (i as f64)*dt)
        .collect::<VecDeque<_>>();

    // Don't be intimidated
    // An Arc just allows us to safely share a
    // pointer between threads.
    // We then have a mutex inside so we can be sure
    // only one thread is accessing the deque at once.
    let ts = Arc::new(Mutex::new(ts));

    let (data_s, data_r) = channel();

    for _ in 0..NUM_THREADS{

        let ts = ts.clone();
        let data_s = data_s.clone();

        thread::spawn(move || {
            let mut sys = System::new(N as usize);

            while let Some(t) = next_t(&ts) {
                if t % 0.1 < 0.01 {
                    println!("{}", t);
                }
                let (U, M) =
                    ensemble_average(&mut sys, 1./t, 1000, 100);

                data_s.send((t, M, U)).unwrap();
            }

        });
    }

    let mut f = File::create("met.dat")
        .expect("Unable to open data file.");

    writeln!(&mut f, "T,M,U").unwrap();

    for _ in 0..STEPS {
        let (t, M, U) = data_r.recv().unwrap();

        writeln!(f, "{}, {}, {}", t, M, U)
            .ok()
            .expect("Unable to write to file.");
    }

}
