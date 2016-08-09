import numpy as np
from numba import jitclass, int64, jit

from .lattice import Lattice, make_lattice

simulations_spec = [
    ('N', int64), ('energy', int64), ('magnetization', int64),
    ('lattice', Lattice.class_type.instance_type)
]


@jitclass(simulations_spec)
class Simulation(object):
    def __init__(self, N, lat):
        self.N = N
        self.lattice = lat
        self.energy = self.lattice.energy
        self.magnetization = self.lattice.magnetization

    def metropolis_step(self, beta, i, j, r):
        # i = randint(self.N)
        # j = randint(self.N)

        dE = self.lattice.delta_E(i, j)
        if dE < 0 or (r < np.exp(-beta * dE)):
            spin = self.lattice.s[i, j]
            self.lattice.s[i, j] = -spin
            self.magnetization -= 2 * spin
            self.energy += dE

    def evolve(self, n, beta, random_sites, rand_accept):
        for i in range(n):
            self.metropolis_step(beta,
                                 random_sites[i, 0],
                                 random_sites[i, 1],
                                 rand_accept[i])

    def time_average(self, n, beta, random_sites, rand_accept):
        '''Calculates the thermal average of the magnetization and energy.'''
        mag = 0
        en = 0

        for i in range(n):
            self.metropolis_step(beta,
                                 random_sites[i, 0],
                                 random_sites[i, 1],
                                 rand_accept[i])
            mag += self.magnetization
            en += self.energy

        return mag / float(n), en / float(n)

    def ensemble_av(self, beta, n_evolve, n_average, r1, r2):
        # Lets the system equilibrate for a while
        self.evolve(n_evolve, beta, r1[:n_evolve], r2[:n_evolve])

        # Now take thermal averages
        return self.time_average(n_average, beta, r1[n_evolve:], r2[n_evolve:])


@jit
def make_simulation(N):
    l = make_lattice(N)
    return Simulation(N, l)
