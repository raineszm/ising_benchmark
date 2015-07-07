extern crate rand;

use rand::{XorShiftRng, Rng};
use rand::distributions::{Range, IndependentSample};

use std::fs::File;
use std::io::Write;

struct System {
    n: usize,
    sites: Vec<i32>,
    nnplus: Vec<usize>,
    nnminus: Vec<usize>,
    rng: XorShiftRng,
    site_range: Range<usize>,
}

impl System {
    fn new(n: usize) -> System {
        System {
            n: n,
            sites: vec![1; n*n],
            nnplus: (0..n).map(|i| (i + 1) % n).collect(),
            nnminus: (0..n).map(|i| (n + i - 1) % n).collect(),
            rng: rand::weak_rng(),
            site_range: Range::new(0, n),
        }
    }

    fn at(&self, i: usize, j: usize) -> i32 {
        self.sites[j + self.n*i]
    }

    fn flip(&mut self, i: usize, j: usize) -> i32 {
        let s = self.at(i, j);
        self.sites[j + self.n*i] = -s;
        s
    }

    fn sum_neighbors(&self, i: usize, j:usize) -> i32 {
       self.at(i, self.nnplus[j])
            + self.at(i, self.nnminus[j])
            + self.at(self.nnplus[i], j)
            + self.at(self.nnminus[i], j)
    }

    fn magnetization(&self) -> i32 {
        self.sites.iter().fold(0, |acc, &s| acc + s)
    }

    fn random_site(&mut self) -> usize {
        self.site_range.ind_sample(&mut self.rng)
    }

    #[allow(non_snake_case)]
    fn deltaE(&self, i: usize, j:usize) -> i32 {
        2*self.at(i, j)*self.sum_neighbors(i, j)
    }
}


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

fn evaluate_energy(sys: &System) -> i32 {
    let n = sys.n;
    let mut total = 0;

    for i in 0..n {
        for j in 0..n {
            total -= sys.at(i,j)*(
                sys.at(sys.nnplus[i], j) +
                sys.at(i, sys.nnplus[j]))
        }
    }
    total
}

#[allow(non_snake_case)]
fn time_average(sys: &mut System,
                n: i32,
                beta: f64) -> (f64, f64) {

    let mut U = evaluate_energy(sys);
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
