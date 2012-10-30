#cython: boundscheck = False
#cython: wraparound = False
cimport cython
cimport numpy as np
import numpy as np
from time import time

from libc.math cimport exp, fmod
from libc.stdio cimport printf
from libc.stdlib cimport rand, RAND_MAX, srand


cdef float fRAND_MAX = <float> RAND_MAX

DTYPE = np.int
ctypedef np.int_t DTYPE_t

utype = np.uint32
ctypedef np.uint32_t UTYPE_t

cdef int N = 50

cdef DTYPE_t energy = 0
cdef DTYPE_t magnetization = 0

cdef np.ndarray nn1, nn2
cdef UTYPE_t [::1] nn1t
cdef UTYPE_t [::1] nn2t

cdef float mean_energy = 0
cdef float mean_magnetization = 0

#Spins are 1 or -1
cdef np.ndarray spins = np.ones((N,N), dtype=DTYPE)

cdef inline unsigned int randint(int n):
    return rand() % n

cdef inline float randfloat():
    return <float>rand()/fRAND_MAX


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int delta_E(DTYPE_t [:,::1] s, unsigned int i, unsigned int j):    
    return 2*s[i,j]*(s[<unsigned int>nn1t[i],j] 
                       + s[i,<unsigned int>nn1t[j]]
                       + s[<unsigned int>nn2t[i],j]
                       + s[i, <unsigned int>nn2t[j]])
    
cdef inline bint accepted(int dE, float beta):
    global rnd_j
    cdef float P = exp(-beta*dE)
    return randfloat() < P

@cython.boundscheck(False)
cdef metropolis_step(DTYPE_t [:,::1] s, float beta):
    cdef unsigned int i, j
    cdef DTYPE_t spin
    global magnetization
    global energy

    i = randint(N)
    j = randint(N)

    cdef int dE = delta_E(s,i,j)
    if dE < 0 or accepted(dE, beta):
        spin = s[i,j]
        s[i,j] = - spin
        magnetization -= 2*spin
        energy += dE


cdef evolve(DTYPE_t [:, ::1] s, int n, float beta):
    cdef int i
    for i in range(n):
        metropolis_step(s, beta)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef time_average(DTYPE_t [:,::1] s, int n, float beta):
    cdef float mag = 0
    cdef float en = 0
    cdef int i
    for i in range(n):
        metropolis_step(s, beta)
        mag += magnetization
        en += energy

    return mag/n, en/n


cdef int evaluate_energy(DTYPE_t [:, ::1] s):
    cdef unsigned int i
    cdef unsigned int j
    cdef DTYPE_t total = 0

    for i in range(N):
        for j in range(N):
            total -= s[i,j]*(s[<unsigned int> nn1t[i],j] + s[i,<unsigned int> nn2t[j]])

    return total


def _ensemble_av(float beta, int n_evolve=1000, int n_average=100):
    cdef float T = 1/beta
    if fmod(T, 0.5) < 0.01:
        printf("%f\n", T)

    evolve(spins, n_evolve, beta)
    return time_average(spins, n_average, beta)

def initialize(int n):
    global spins
    global magnetization
    global energy
    global N
    N = n
    srand(<unsigned int>time())

    spins = np.ones((N,N), dtype=DTYPE)
    populate_nn()

    magnetization = np.sum(spins)
    energy = evaluate_energy(spins)

cdef populate_nn():
    global nn1, nn2, nn1t, nn2t
    cdef np.ndarray[UTYPE_t] r = np.arange(N, dtype=utype)
    nn1 = np.roll(r, 1)
    nn2 = np.roll(r, -1)
    nn1t = nn1
    nn2t = nn2

ensemble_av = np.frompyfunc(_ensemble_av, 3, 2)







