extern crate rand;

mod system;

use system::System;

use rand::Rng;
use std::fs::File;
use std::io::Write;



#[allow(non_snake_case)]
fn metropolis_step(sys: &mut System,
                   beta: f64) -> (i32, i32) {
    let i = sys.random_site();
    let j = sys.random_site();

    let dE = sys.deltaE(i, j);

    if dE < 0 || (sys.rng.next_f64() < (-beta*(dE as f64)).exp() ){
        return (dE, -2*sys.flip(i, j));
    }

    (0, 0)
}

fn evolve(sys: &mut System,
          n: i32,
          beta: f64) {

    for _i in 0..n {
        metropolis_step(sys, beta);
    }
}

#[allow(non_snake_case)]
fn time_average(sys: &mut System,
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
        M_tot += M as i64;
    }

    let en = (U_tot as f64) / (n as f64);
    let mag = (M_tot as f64) / (n as f64);

    (en, mag)
}

fn ensemble_average(sys: &mut System,
                    beta: f64,
                    n_evolve: i32,
                    n_average: i32) -> (f64, f64) {
    evolve(sys, n_evolve, beta);
    time_average(sys, n_average, beta)
}

#[allow(non_snake_case)]
fn main() {
    let N: i32 = 64;
    let mut sys = System::new(N as usize);
    let dt: f64 = (5. - 0.1)/400.;

    let ts = (0..400).map( |i|
                           0.1 + (i as f64)*dt
                           ).collect::<Vec<_>>();


    let data = ts.iter().map( |t| {
        let (U, M) = ensemble_average(&mut sys, 1./t, 1000*N*N, 100*N*N);

        if t % 0.1 < 0.01 {
            println!("T: {}", t);
        }

        (t.clone(), M, U)
    }).collect::<Vec<_>>();

    let mut f = File::create("met.dat")
        .ok()
        .expect("Unable to open data file.");

    writeln!(&mut f, "T,M,U").unwrap();

    for point in &data {
        let (t, M, U) = *point;
        writeln!(f, "{}, {}, {}", t, M, U).unwrap();
    }

}
