"""
.. module:: lattice

.. moduleauthor:: Zachary Raines

"""
import array
class Lattice:
    def __init__(self, N):
        """
        Creates a new lattice.

        :param N: The size of the lattice.
        :type N: int
        :returns: :class:`Lattice`
        """
        r = list(range(N))
        self.N = N
        self.s = array.array('i', [1]*(N**2))
        self.nn1 = r[1:] + [r[0]]
        self.nn2 = [r[-1]] + r[:-1]

    def __getitem__(self, key):
        i, j = key
        return self.s[j + self.N*i]


    def __setitem__(self, key, value):
        i, j = key
        self.s[j + self.N*i] = value

    def delta_E(self, i, j):
        """Calculate the energy change due to a spin flip"""
        return (
            2
            * self[i, j]
            * (
                self[self.nn1[i], j]
                + self[i, self.nn1[j]]
                + self[self.nn2[i], j]
                + self[i, self.nn2[j]]
            )
        )

    def magnetization(self):
        return sum(self.s)

    def energy(self):
        total = 0

        for i in range(self.N):
            for j in range(self.N):
                total -= self[i, j] * (
                    self[self.nn1[i], j] + self[i, self.nn2[j]]
                )
        return total

    def flip(self, i, j):
        spin = self[i, j]
        self[i, j] = -spin
        return spin

    def push_neighbors(self, i, j, queue):
        queue.append((i, self.nn1[j]))
        queue.append((i, self.nn2[j]))
        queue.append((self.nn1[i], i))
        queue.append((self.nn2[i], i))
