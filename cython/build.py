import os
import sys

import numpy as np

from Cython.Build import cythonize

# use cythonize to build the extensions
modules = ["ising/simulate.pyx"]

extensions = cythonize(modules)


def build(setup_kwargs):
    """Needed for the poetry building interface."""

    setup_kwargs.update({"ext_modules": extensions, "include_dirs": [np.get_include()]})
