import cffi

ffi = cffi.FFI()
ffi.set_source('ising._random', "#include <stdlib.h>")
ffi.cdef("""#define RAND_MAX ...
               int rand(void);
               void srand(int);""")
ffi.compile()
