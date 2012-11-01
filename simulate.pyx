#cython: boundscheck = False
#cython: wraparound = False
cimport cython
cimport numpy as np

#C Library Imports
from libc.math cimport exp, fmod
from libc.stdio cimport printf
from libc.stdlib cimport rand, RAND_MAX, srand

#Python libarary imports
import numpy as np
from time import time

#Type definitions
DTYPE = np.int
ctypedef np.int_t DTYPE_t
utype = np.uint32
ctypedef np.uint32_t UTYPE_t

#Lattice size: Lattice is NxN
cdef int N = 50

#The lattice
#Spins are 1 or -1
cdef np.ndarray spins = np.ones((N,N), dtype=DTYPE)

#Thermodynamic variables
cdef DTYPE_t energy = 0
cdef DTYPE_t magnetization = 0

#Nearest neighbor maps
#These nn1 contains positive nearest neighbor indices
#and nn2 negative
cdef np.ndarray nn1, nn2
cdef UTYPE_t [::1] nn1v
cdef UTYPE_t [::1] nn2v


#---------------- RANDOM NUMBER GENERATION ----------------------------
cdef float fRAND_MAX = <float> RAND_MAX

#Generate a random integer less than n
cdef inline unsigned int randint(int n):
    return rand() % n

#Generate a random number between 0 and 1
cdef inline float randfloat():
    return <float>rand()/fRAND_MAX

#------------------------------- END ----------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline DTYPE_t delta_E(DTYPE_t [:,::1] s, unsigned int i, unsigned int j):    
    '''Calculate the energy change due to a spin flip'''
    return 2*s[i,j]*(s[<unsigned int>nn1v[i],j] 
                       + s[i,<unsigned int>nn1v[j]]
                       + s[<unsigned int>nn2v[i],j]
                       + s[i, <unsigned int>nn2v[j]])
    
cdef inline bint accepted(DTYPE_t dE, float beta):
    '''Determine whether a positive energy change is accepted'''
    return randfloat() < exp(-beta*dE) 

@cython.boundscheck(False)
@cython.wraparound(False)
cdef metropolis_step(DTYPE_t [:,::1] s, float beta):
    cdef unsigned int i, j
    cdef DTYPE_t spin
    global magnetization
    global energy

    i = randint(N)
    j = randint(N)

    cdef DTYPE_t dE = delta_E(s,i,j)
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
cdef time_average(DTYPE_t [:,::1] s, int n, beta):
    '''Calculates the thermal average of the magnetization and energy.'''
    cdef double mag = 0
    cdef double en = 0
    cdef double fn = <double> n
    cdef int i
    cdef int j

    #implement variance checking
    for i in range(n):
        metropolis_step(s, beta)
        mag += magnetization
        en += energy

    return mag/fn, en/fn


cdef int evaluate_energy(DTYPE_t [:, ::1] s):
    cdef unsigned int i
    cdef unsigned int j
    cdef DTYPE_t total = 0

    for i in range(N):
        for j in range(N):
            total -= s[i,j]*(s[<unsigned int> nn1v[i],j] + s[i,<unsigned int> nn2v[j]])

    return total


def _ensemble_av(float beta, int n_evolve=1000, int n_average=100):
    #Updates us periodically on the simulation progress
    cdef float T = 1/beta
    if fmod(T, 0.1) < 0.01:
        printf("%f\n", T)

    #Lets the system equilibrate for a while
    evolve(spins, n_evolve, beta)

    #Now take thermal averages
    return time_average(spins, n_average, beta)

def initialize(int n, double var=VAR):
    global spins
    global magnetization
    global energy
    global N

    N = n
    srand(<unsigned int>time())

    spins = np.ones((N,N), dtype=DTYPE)

    populate_nn()

    #We calculate the total energy and magnetization once,
    #then update it at each iteration step
    magnetization = np.sum(spins)
    energy = evaluate_energy(spins)

cdef populate_nn():
    '''Prepopulates an array of nearest neighbor indices, to speed utype
    calculation of the energy change.'''

    global nn1, nn2, nn1v, nn2v
    cdef np.ndarray[UTYPE_t] r = np.arange(N, dtype=utype)
    nn1 = np.roll(r, 1)
    nn2 = np.roll(r, -1)
    nn1v = nn1
    nn2v = nn2

#Convert our ensemble average to a function on arrays
ensemble_av = np.frompyfunc(_ensemble_av, 3, 2)
