#!/usr/bin/env julia
import Metropolis
require("setup.jl")


tic()
data = pmap(T) do t
    (M, U) = Metropolis.ensemble_av(1/t, 1000N^2, 100N^2)
    (t, M, U)
end
toc()

sort!(data)

out = open("met.dat", "w")
@printf(out, "#T\tM\tU")
for (t, M, U) in data
    @printf(out, "%f\t%f\t%f\n", t, abs(M), U)
end
close(out)

