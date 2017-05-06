#!/usr/bin/env julia
push!(LOAD_PATH, "src")
import Metropolis.Runner

tic()
Metropolis.Runner.main()
toc()
