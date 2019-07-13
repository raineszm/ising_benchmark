#!/usr/bin/env python
import numpy as np

import ising.run

N = 64
T = np.linspace(0.1, 5, 400)
ising.run.psimulate("data.csv", T, N, 1000 * N ** 2, 100 * N ** 2)
