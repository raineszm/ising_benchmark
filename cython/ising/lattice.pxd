cimport numpy as np
from libcpp.deque cimport deque


#Type definitions
ctypedef np.int_t DTYPE_t
ctypedef np.uint32_t UTYPE_t
ctypedef deque[(int, int)] IDX_QUEUE_t

cdef class Lattice:
    cdef int N

    cdef np.ndarray s
    #Nearest neighbor maps
    #These nn1 contains positive nearest neighbor indices
    #and nn2 negative
    cdef np.ndarray nn1, nn2
    cdef UTYPE_t [::1] nn1v
    cdef UTYPE_t [::1] nn2v

    cdef int delta_E(self, int i, int j)

    cdef int magnetization(self)

    cdef int energy(self)

    cdef int flip(self, int i, int j)

    cdef void push_neighbors(self, int i, int j, IDX_QUEUE_t queue)
