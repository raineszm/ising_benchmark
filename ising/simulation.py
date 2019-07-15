from collections import deque

import random
import math
from .lattice import Lattice


def randint(n):
    return random.randrange(0, n)


def randfloat():
    return random.random()


def accepted(dE, beta):
    """Determine whether a positive energy change is accepted"""
    return randfloat() < math.exp(-beta * dE)


class Simulation:
    def __init__(self, N):
        self.N = N
        self.lattice = Lattice(N)
        self.energy = self.lattice.energy()
        self.magnetization = self.lattice.magnetization()

    def metropolis_step(self, beta):
        i = randint(self.N)
        j = randint(self.N)

        neighbors = deque()

        flip_prob = 1 - math.exp(-2 * beta)

        s = self.lattice.flip(i, j)
        dE = self.lattice.delta_E(i, j)
        dM = -2 * s

        self.lattice.push_neighbors(i, j, neighbors)

        while neighbors:
            i, j = neighbors.popleft()

            if self.lattice[i, j] == s and randfloat() < flip_prob:
                dM -= 2 * s
                dE += self.lattice.delta_E(i, j)

                self.lattice.flip(i, j)

                self.lattice.push_neighbors(i, j, neighbors)

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

        return math.sqrt(mag / float(n)), en / float(n)

    def ensemble_av(self, beta, n_evolve, n_average):
        # Lets the system equilibrate for a while
        self.evolve(n_evolve, beta)

        # Now take thermal averages
        return self.time_average(n_average, beta)
