from collections import OrderedDict

import numpy as np
from numba import int64, jitclass, njit, typeof

from .fifo import Fifo
from .lattice import Lattice


@njit
def randint(n):
    return np.random.randint(0, n)


@njit
def randfloat():
    return np.random.random()


@njit
def accepted(dE, beta):
    """Determine whether a positive energy change is accepted"""
    return randfloat() < np.exp(-beta * dE)


@jitclass(
    OrderedDict(
        {
            "N": int64,
            "lattice": typeof(Lattice(1)),
            "energy": int64,
            "magnetization": int64,
            "neighbors": typeof(Fifo(1)),
        }
    )
)
class Simulation:
    def __init__(self, N):
        self.N = N
        self.lattice = Lattice(N)
        self.energy = self.lattice.energy()
        self.magnetization = self.lattice.magnetization()
        self.neighbors = Fifo(N ** 2)

    def metropolis_step(self, beta):
        i = randint(self.N)
        j = randint(self.N)

        flip_prob = 1 - np.exp(-2 * beta)

        s = self.lattice.flip(i, j)
        dE = self.lattice.delta_E(i, j)
        dM = -2 * s

        self.lattice.push_neighbors(i, j, self.neighbors)

        while not self.neighbors.empty:
            i, j = self.neighbors.pop()

            if self.lattice.s[i, j] == s and randfloat() < flip_prob:
                dM -= 2 * s
                dE += self.lattice.delta_E(i, j)

                self.lattice.flip(i, j)

                self.lattice.push_neighbors(i, j, self.neighbors)

        self.magnetization += dM
        self.energy += dE

    def evolve(self, n, beta):
        for i in range(n):
            self.metropolis_step(beta)

    def time_average(self, n, beta):
        """Calculates the thermal average of the magnetization and energy."""
        mag = 0
        en = 0

        for i in range(n):
            self.metropolis_step(beta)
            mag += self.magnetization ** 2
            en += self.energy

        return np.sqrt(mag / float(n)), en / float(n)

    def ensemble_av(self, beta, n_evolve, n_average):
        # Lets the system equilibrate for a while
        self.evolve(n_evolve, beta)

        # Now take thermal averages
        return self.time_average(n_average, beta)
