#cython: boundscheck = False
#cython: wraparound = False
cimport cython
cimport numpy as np
from .lattice cimport Lattice

#C Library Imports
from libc.math cimport exp, sqrt
from libc.stdlib cimport rand, RAND_MAX, srand

from libcpp.deque cimport deque

#Python libarary imports
from time import time

#Type definitions
ctypedef np.int_t DTYPE_t
ctypedef np.uint32_t UTYPE_t
ctypedef deque[(int, int)] IDX_QUEUE_t


#---------------- RANDOM NUMBER GENERATION ----------------------------
cdef double fRAND_MAX = <double> RAND_MAX

#Generate a random integer less than n
cdef inline unsigned int randint(int n):
    return rand() % n

#Generate a random number between 0 and 1
cdef inline double randdouble():
    return <double>rand()/fRAND_MAX

#------------------------------- END ----------------------------------


cdef class Simulation:
    cdef Lattice lattice
    #Thermodynamic variables
    cdef DTYPE_t energy
    cdef DTYPE_t magnetization
    cdef IDX_QUEUE_t neighbors

    def __init__(self, int N):
        srand(<unsigned int>time())
        self.N = N
        self.lattice = Lattice(N)
        self.energy = self.lattice.energy()
        self.magnetization = self.lattice.magnetization()

    def metropolis_step(self, double beta):
        cdef UTYPE_t i = randint(self.N)
        cdef UTYPE_t j = randint(self.N)

        cdef double flip_prob = 1.0 - exp(-2 * beta)

        cdef int s = self.lattice.flip(i, j)
        cdef int dE = self.lattice.delta_E(i, j)
        cdef int dM = -2 * s

        self.lattice.push_neighbors(i, j, self.neighbors)

        while not self.neighbors.empty():
            i, j = self.neighbors.pop_front()

            if self.lattice.s[i, j] == s and randdouble() < flip_prob:
                dM -= 2 * s
                dE += self.lattice.delta_E(i, j)

                self.lattice.flip(i, j)

                self.lattice.push_neighbors(i, j, self.neighbors)

        self.magnetization += dM
        self.energy += dE

    def evolve(self, int n, double beta):
        for i in range(n):
            self.metropolis_step(beta)

    def time_average(self, int n, double beta):
        """Calculates the thermal average of the magnetization and energy."""
        cdef int mag = 0
        cdef int en = 0

        for i in range(n):
            self.metropolis_step(beta)
            mag += self.magnetization ** 2
            en += self.energy

        return sqrt(mag / <double>(n)), en / <double>(n)

    def ensemble_av(self, double beta, int n_evolve, int n_average):
        # Lets the system equilibrate for a while
        self.evolve(n_evolve, beta)

        # Now take thermal averages
        return self.time_average(n_average, beta)
