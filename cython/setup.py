from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    name = "Metropolis",
    include_dirs = [numpy.get_include()],
    ext_modules = cythonize('simulate.pyx'), # accepts a glob pattern
)