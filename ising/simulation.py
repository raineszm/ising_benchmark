from collections import OrderedDict

import numpy as np
from numba import int64, jitclass, njit, typeof

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
        }
    )
)
class Simulation:
    def __init__(self, N):
        self.N = N
        self.lattice = Lattice(N)
        self.energy = self.lattice.energy()
        self.magnetization = self.lattice.magnetization()

    def metropolis_step(self, beta):
        i = randint(self.N)
        j = randint(self.N)

        dE = self.lattice.delta_E(i, j)
        if dE < 0 or accepted(dE, beta):
            spin = self.lattice.s[i, j]
            self.lattice.s[i, j] = -spin
            self.magnetization -= 2 * spin
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
            mag += self.magnetization
            en += self.energy

        return mag / float(n), en / float(n)

    def ensemble_av(self, beta, n_evolve, n_average):
        # Lets the system equilibrate for a while
        self.evolve(n_evolve, beta)

        # Now take thermal averages
        return self.time_average(n_average, beta)
