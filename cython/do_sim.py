#!/usr/bin/env python3
from numpy import *

import simulate

N = 64

simulate.initialize(N)
T = linspace(0.1, 5, 400)
M, U = simulate.ensemble_av(1 / T, 1000 * N ** 2, 100 * N ** 2)
savez("data-vari.npz", T=T, M=M, U=U, N=N)
