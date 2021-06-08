#!/usr/bin/env julia
import Distributed: @everywhere
@everywhere push!(LOAD_PATH, "src")
@everywhere import Ising
using BenchmarkTools: @btime


data_file = "data.csv"
if length(ARGS) > 0
    data_file = ARGS[1]
end

@btime Ising.Metropolis.Runner.main(data_file)
