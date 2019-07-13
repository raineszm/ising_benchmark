from collections import OrderedDict

import numpy as np
from numba import int64, jitclass, njit, typeof

from .lattice import Lattice


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


@njit
def metropolis_step(sim, beta):
    i = randint(sim.N)
    j = randint(sim.N)

    dE = sim.lattice.delta_E(i, j)
    if dE < 0 or accepted(dE, beta):
        spin = sim.lattice.s[i, j]
        sim.lattice.s[i, j] = -spin
        sim.magnetization -= 2 * spin
        sim.energy += dE


@njit
def evolve(sim, n, beta):
    for i in range(n):
        metropolis_step(sim, beta)


@njit
def time_average(sim, n, beta):
    """Calculates the thermal average of the magnetization and energy."""
    mag = 0
    en = 0

    for i in range(n):
        metropolis_step(sim, beta)
        mag += sim.magnetization
        en += sim.energy

    return mag / float(n), en / float(n)


@njit
def ensemble_av(sim, beta, n_evolve, n_average):
    # Lets the system equilibrate for a while
    evolve(sim, n_evolve, beta)

    # Now take thermal averages
    return time_average(sim, n_average, beta)


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
