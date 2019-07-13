"""
.. module:: lattice

.. moduleauthor:: Zachary Raines

"""
from collections import OrderedDict

import numpy as np
from numba import int64, jit, jitclass, njit


@jitclass({"N": int64, "s": int64[:, :], "nn1": int64[:], "nn2": int64[:]})
class Lattice:
    def __init__(self, N):
        """
        Creates a new lattice.

        :param N: The size of the lattice.
        :type N: int
        :returns: :class:`Lattice`
        """
        r = np.arange(N)
        self.N = N
        self.s = np.ones((N, N), np.int_)
        self.nn1 = np.roll(r, 1)
        self.nn2 = np.roll(r, -1)

    def delta_E(self, i, j):
        """Calcuselfe the energy change due to a spin flip"""
        return (
            2
            * self.s[i, j]
            * (
                self.s[self.nn1[i], j]
                + self.s[i, self.nn1[j]]
                + self.s[self.nn2[i], j]
                + self.s[i, self.nn2[j]]
            )
        )

    def magnetization(self):
        return np.sum(self.s)

    def energy(self):
        total = 0

        for i in range(self.N):
            for j in range(self.N):
                total -= self.s[i, j] * (
                    self.s[self.nn1[i], j] + self.s[i, self.nn2[j]]
                )
        return total

    def flip(self, i, j):
        spin = self.s[i, j]
        self.s[i, j] = -spin
        return spin

    def push_neighbors(self, i, j, queue):
        queue.append((i, self.nn1[j]))
        queue.append((i, self.nn2[j]))
        queue.append((self.nn1[i], i))
        queue.append((self.nn2[i], i))
