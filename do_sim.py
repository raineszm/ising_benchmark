#!/usr/bin/env python3
import simulate
from numpy import *

N=64

simulate.initialize(N)
T = linspace(0.1, 5, 400)
M, U = simulate.ensemble_av(1/T, 1000*N**2, 100*N**2)
# T1 = linspace(0.1, 2.2, 130, endpoint=False)
# T2 = linspace(2.2, 2.5, 300, endpoint=False)
# T3 = linspace(2.5, 5, 350)
# M1, U1 = simulate.ensemble_av(1/T1, 10*N**2, 10*N**2)
# M2, U2 = simulate.ensemble_av(1/T2, 300*N**2, 50*N**2)
# M3, U3 = simulate.ensemble_av(1/T3, 10*N**2, 10*N**2)
# T = concatenate((T1, T2, T3))
# M = concatenate((M1, M2, M3))
# U = concatenate((U1, U2, U3))
savez("data-vari.npz", T=T, M=M, U=U, N=N)
# savez("data-fine.npz", T1=T1, T2=T2, T3=T3, U1=U1, U2=U2, U3=U3, M1=M1, M2=M2, M3=M3)
