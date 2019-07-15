#!/usr/bin/env python
import sys

import numpy as np

import ising.run

N = 64
T = np.linspace(0.1, 5, 400)
fname = "data.csv"

if len(sys.argv) > 1:
    fname = sys.argv[1]

ising.run.psimulate(fname, T, N, 1000, 100)
