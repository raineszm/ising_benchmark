extern crate rand;

use self::rand::XorShiftRng;
use self::rand::distributions::{Range, IndependentSample};

pub struct System {
    n: usize,
    sites: Vec<i32>,
    nnplus: Vec<usize>,
    nnminus: Vec<usize>,
    pub rng: XorShiftRng,
    site_range: Range<usize>,
}

impl System {
    pub fn new(n: usize) -> System {
        System {
            n: n,
            sites: vec![1; n*n],
            nnplus: (0..n).map(|i| (i + 1) % n).collect(),
            nnminus: (0..n).map(|i| (n + i - 1) % n).collect(),
            rng: rand::weak_rng(),
            site_range: Range::new(0, n),
        }
    }

    pub fn at(&self, i: usize, j: usize) -> i32 {
        self.sites[j + self.n*i]
    }

    pub fn flip(&mut self, i: usize, j: usize) -> i32 {
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

    pub fn evaluate_energy(&self) -> i32 {
        let mut total = 0;

        for i in 0..self.n {
            for j in 0..self.n {
                total -= self.at(i, j)*(
                    self.at(self.nnplus[i], j) +
                    self.at(i, self.nnplus[j]))
            }
        }
        total
    }

    pub fn magnetization(&self) -> i32 {
        self.sites.iter().fold(0, |acc, &s| acc + s)
    }

    pub fn random_site(&mut self) -> usize {
        self.site_range.ind_sample(&mut self.rng)
    }

    #[allow(non_snake_case)]
    pub fn deltaE(&self, i: usize, j:usize) -> i32 {
        2*self.at(i, j)*self.sum_neighbors(i, j)
    }
}
