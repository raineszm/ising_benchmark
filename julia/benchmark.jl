#!/usr/bin/env julia
push!(LOAD_PATH, "src")
import Ising
using BenchmarkTools: @btime


data_file = "data.csv"
if length(ARGS) > 0
    data_file = ARGS[1]
end

@btime Ising.Runner.main(data_file)
