#cython: boundscheck = False
#cython: wraparound = False
"""
.. module:: lattice

.. moduleauthor:: Zachary Raines

"""
cimport numpy as np
import numpy as np
DTYPE = np.int
UTYPE = np.uint32

cdef class Lattice:

    def __init__(self, int N):
        """
        Creates a new lattice.

        :param N: The size of the lattice.
        :type N: int
        :returns: :class:`Lattice`
        """
        cdef np.ndarray[UTYPE_t] r = np.arange(N, dtype=UTYPE)
        self.N = N
        self.s = np.ones((N, N), dtype=DTYPE)
        self.nn1 = np.roll(r, 1)
        self.nn2 = np.roll(r, -1)
        self.nn1v
        self.nn1v = self.nn1
        self.nn2v = self.nn2

    cdef int delta_E(self, int i, int j):
        """Calcuselfe the energy change due to a spin flip"""
        return 2*self.s[i,j]*(self.s[<unsigned int>self.nn1v[i],j]
                        + self.s[i,<unsigned int>self.nn1v[j]]
                        + self.s[<unsigned int>self.nn2v[i],j]
                        + self.s[i, <unsigned int>self.nn2v[j]])

    cdef int magnetization(self):
        return sum(self.s)

    cdef int energy(self):
        cdef unsigned int i
        cdef unsigned int j
        cdef DTYPE_t total = 0

        for i in range(self.N):
            for j in range(self.N):
                total -= self.s[i,j]*(self.s[<unsigned int> self.nn1v[i],j] + self.s[i,<unsigned int> self.nn2v[j]])

        return total


    cdef int flip(self, int i, int j):
        cdef int spin = self.s[i, j]
        self.s[i, j] = -spin
        return spin

    cdef void push_neighbors(self, int i, int j, IDX_QUEUE_t queue):
        queue.push_back((i, self.nn1v[j]))
        queue.push_back((i, self.nn2v[j]))
        queue.push_back((self.nn1v[i], i))
        queue.push_back((self.nn2v[i], i))
