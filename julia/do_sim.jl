#!/usr/bin/env julia
import Distributed: @everywhere
@everywhere push!(LOAD_PATH, "src")
@everywhere import Metropolis

data_file = "data.csv"
println(ARGS)
if length(ARGS) > 0
    data_file = ARGS[1]
end

Metropolis.Runner.main(data_file)
