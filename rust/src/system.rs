extern crate rand;

use rand::{SeedableRng, RngCore};
use rand::rngs::SmallRng;

/// Description of the state of a 2D Ising System
///
/// Describes a 2D square lattice with classical Ising spins
/// interacting via a nearest neighbor ferromagnetic exchange.
pub struct System {
    /// The lattice is `n` x `n` sites.
    n: usize,
    /// State of the system. Each spin is either up (+1) or down (-1)
    sites: Vec<i32>,
    /// Lookup table of nearest neighbors
    nnplus: Vec<usize>,
    /// Lookup table of nearest neighbors
    nnminus: Vec<usize>,
    /// Random number generator for picking random sites.
    pub rng: SmallRng,
}

impl System {
    /// Create a new `n` by `n` [System]()
    pub fn new(n: usize) -> System {
        System {
            n: n,
            sites: vec![1; n*n],
            nnplus: (0..n).map(|i| (i + 1) % n).collect(),
            nnminus: (0..n).map(|i| (n + i - 1) % n).collect(),
            rng: SmallRng::from_entropy(),
        }
    }

    /// Get the value of the spin located at `(i, j)`
    ///
    /// Spins are either up (+1) or down (-1)
    pub fn at(&self, i: usize, j: usize) -> i32 {
        self.sites[j + self.n*i]
    }

    /// Flip the value of the spin located at `(i, j)`
    ///
    /// Spins are either up (+1) or down (-1).
    /// The old value of the spin is returned.
    pub fn flip(&mut self, i: usize, j: usize) -> i32 {
        let s = self.at(i, j);
        self.sites[j + self.n*i] = -s;
        s
    }

    pub fn nnplus(&self, i: usize) -> usize {
        self.nnplus[i]
    }

    pub fn nnminus(&self, i: usize) -> usize {
        self.nnminus[i]
    }

    fn sum_neighbors(&self, i: usize, j:usize) -> i32 {
       self.at(i, self.nnplus[j])
            + self.at(i, self.nnminus[j])
            + self.at(self.nnplus[i], j)
            + self.at(self.nnminus[i], j)
    }

    /// Evaluate the total energy of the system.
    ///
    /// The energy of the system is the product of neighboring sites
    /// summed over all links.
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

    /// Calculate the total magnetization of the system.
    ///
    /// This is just the sum of all spins.
    pub fn magnetization(&self) -> i32 {
        self.sites.iter().fold(0, |acc, &s| acc + s)
    }

    /// Returns a random index along one linear dimension.
    ///
    /// If one wants to pick a random spin in the lattice it can
    /// be done by picking a random row and column
    ///
    /// ```
    /// let lat = Lattice(64);
    /// let i = lat.random_site();
    /// let j = lat.random_site();
    /// let site = lat.at(i, j);
    /// ```
    pub fn random_site(&mut self) -> usize {
        (self.rng.next_u32()/(u32::MAX/(self.n as u32) + 1)) as usize
    }


    #[allow(non_snake_case)]
    /// Calculate the change in energy due to flipping the spin at `(i, j)`.
    pub fn deltaE(&self, i: usize, j:usize) -> i32 {
        2*self.at(i, j)*self.sum_neighbors(i, j)
    }
}
