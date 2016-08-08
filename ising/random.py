from numba import cffi_support

from . import _random
rand = _random.lib.rand
RAND_MAX = _random.lib.RAND_MAX
cffi_support.register_module(_random)
