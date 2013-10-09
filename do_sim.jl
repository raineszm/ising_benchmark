#!/usr/bin/env julia
import Metropolis

N=32
Metropolis.initialize(N)
T = linspace(0.1, 5, 400)
M = similar(T)
U = similar(T)

out = open("met.dat", "w")
@printf(out, "#T\tM\tU")

for (i,t) in enumerate(T)
    (M[i], U[i]) = Metropolis.ensemble_av(1/t, 100N^2, 100N^2)
    @printf(out, "%f\t%f\t%f\n", t, M[i], U[i])
end

close(out)

