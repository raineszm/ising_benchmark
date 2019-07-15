#!/usr/bin/env python3
import sys

import numpy as np

import simulate

sys.path.append("ising")


N = 64

simulate.initialize(N)
T = np.linspace(0.1, 5, 400)
M, U = simulate.ensemble_av(1 / T, 1000 * N ** 2, 100 * N ** 2)
np.savez("data-vari.npz", T=T, M=M, U=U, N=N)
