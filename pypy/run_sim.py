#!/usr/bin/env python
import sys

import ising.run

N = 64
TMIN = 0.1
TMAX = 5
STEPS = 400
T = [TMIN + i*(TMAX - TMIN)/STEPS for i in range(STEPS)]


fname = "data.csv"

if len(sys.argv) > 1:
    fname = sys.argv[1]

ising.run.psimulate(fname, T, N, 1000, 100)
