#!/usr/bin/env julia
push!(LOAD_PATH, "src")
import Ising

data_file = "data.csv"
if length(ARGS) > 0
    data_file = ARGS[1]
end

Ising.Runner.main(data_file)
