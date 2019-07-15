#!/usr/bin/env python
import ising.run

N = 64
TMIN = 0.1
TMAX = 5
STEPS = 400
T = [TMIN + i*(TMAX - TMIN)/STEPS for i in range(STEPS)]
ising.run.psimulate("data.csv", T, N, 1000, 100)
