extern crate rand;

use rand::{XorShiftRng, Rng};
use rand::distributions::{Range, IndependentSample};

struct Lattice {
    n: usize,
    sites: Vec<i32>,
    nnplus: Vec<usize>,
    nnminus: Vec<usize>,
    rng: XorShiftRng,
    site_range: Range<usize>,
}

impl Lattice {
    fn new(n: usize) -> Lattice {
        Lattice {
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

}

#[allow(non_snake_case)]
fn deltaE(lat: &Lattice, i: usize, j:usize) -> i32 {
    2*lat.at(i, j)*lat.sum_neighbors(i, j)
}

#[allow(non_snake_case)]
fn metropolis_step(lat: &mut Lattice,
                   beta: f64) -> (i32, i32) {
    let i = lat.random_site();
    let j = lat.random_site();

    let dE = deltaE(lat, i, j);

    if dE < 0 || (lat.rng.next_f64() < (-beta*(dE as f64)).exp() ){
        return (dE, -2*lat.flip(i, j));
    }

    (0, 0)
}

fn evolve(lat: &mut Lattice,
          n: i32,
          beta: f64) {

    for _i in 0..n {
        metropolis_step(lat, beta);
    }
}

fn evaluate_energy(lat: &Lattice) -> i32 {
    let n = lat.n;
    let mut total = 0;

    for i in 0..n {
        for j in 0..n {
            total -= lat.at(i,j)*(
                lat.at(lat.nnplus[i], j) +
                lat.at(i, lat.nnplus[j]))
        }
    }
    total
}

#[allow(non_snake_case)]
fn time_average(lat: &mut Lattice,
                n: i32,
                beta: f64) -> (f64, f64) {

    let mut U = evaluate_energy(lat);
    let mut M = lat.magnetization();


    let mut M_tot: i64 = 0;
    let mut U_tot: i64 = 0;
    let mut dE: i32 = 0;
    let mut dM: i32 = 0;

    for _ in 0..n {
        let (dE, dM) = metropolis_step(lat, beta);

        U += dE;
        M += dM;

        U_tot += U as i64;
        M_tot += M as i64;
    }

    let en = (U_tot as f64) / (n as f64);
    let mag = (M_tot as f64) / (n as f64);

    (en, mag)
}

fn ensemble_average(lat: &mut Lattice,
                    beta: f64,
                    n_evolve: i32,
                    n_average: i32) -> (f64, f64) {
    evolve(lat, n_evolve, beta);
    time_average(lat, n_average, beta)
}

fn main() {
    let N: i32 = 64;
    let mut lat = Lattice::new(N as usize);
    let mut t = 0.1;
    let dt: f64 = (5. - 0.1)/400.;


    for _ in 0..400 {

        let (U, M) = ensemble_average(&mut lat, 1./t, 1000*N*N, 100*N*N);

        if t % 0.1 < 0.01 {
            println!("T: {}", t);
            println!("U: {}", U);
            println!("M: {}", M);
        }

        t += dt;
    }
}
