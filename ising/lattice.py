"""
.. module:: lattice

.. moduleauthor:: Zachary Raines

"""
from collections import namedtuple

import numpy as np
from numba import int64, jit, njit

Lattice = namedtuple("Lattice", ["N", "s", "nn1", "nn2"])


@njit
def make_lattice(N):
    """
    Creates a new lattice.

    :param N: The size of the lattice.
    :type N: int
    :param nn1: Array of nearest neighbor indices
    :param nn2: Array of nearest neighbor indices
    :returns: :class:`Lattice`
    """
    r = np.arange(N)
    s = np.ones((N, N), np.int_)
    nn1 = np.roll(r, 1)
    nn2 = np.roll(r, -1)
    return Lattice(N, s, nn1, nn2)


@njit
def delta_E(lat, i, j):
    """Calculate the energy change due to a spin flip"""
    return (
        2
        * lat.s[i, j]
        * (
            lat.s[lat.nn1[i], j]
            + lat.s[i, lat.nn1[j]]
            + lat.s[lat.nn2[i], j]
            + lat.s[i, lat.nn2[j]]
        )
    )


@njit
def magnetization(lat):
    return np.sum(lat.s)


@njit
def energy(lat):
    total = 0

    for i in range(lat.N):
        for j in range(lat.N):
            total -= lat.s[i, j] * (lat.s[lat.nn1[i], j] + lat.s[i, lat.nn2[j]])
    return total
