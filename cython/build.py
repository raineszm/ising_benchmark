import os
import sys

import numpy as np

from Cython.Build import cythonize

extensions = cythonize("ising/*.pyx")


def build(setup_kwargs):
    """Needed for the poetry building interface."""

    setup_kwargs.update({"ext_modules": extensions, "include_dirs": [np.get_include()]})
