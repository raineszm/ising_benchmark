from numba import jit, jitclass, int64
import numpy as np

lattice_fields = [
    ('N', int64), ('s', int64[:, :]), ('nn1', int64[:]), ('nn2', int64[:])
]


@jitclass(lattice_fields)
class Lattice(object):
    def __init__(self, N, nn1, nn2):
        self.N = N
        self.s = np.ones((N, N), np.int64)
        self.nn1 = nn1
        self.nn2 = nn2

    def delta_E(self, i, j):
        '''Calculate the energy change due to a spin flip'''
        return 2 * self.s[i, j] * (
            self.s[self.nn1[i], j] + self.s[i, self.nn1[j]] +
            self.s[self.nn2[i], j] + self.s[i, self.nn2[j]])

    @property
    def magnetization(self):
        return np.sum(self.s)

    @property
    def energy(self):
        total = 0

        for i in range(self.N):
            for j in range(self.N):
                total -= self.s[i, j] * (
                    self.s[self.nn1[i], j] + self.s[i, self.nn2[j]])
        return total


@jit(Lattice.class_type.instance_type(int64))
def make_lattice(N):
    r = np.arange(N)
    nn1 = np.roll(r, 1)
    nn2 = np.roll(r, -1)
    return Lattice(N, nn1, nn2)
