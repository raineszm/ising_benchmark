#!/usr/bin/env julia
import Distributed: @everywhere
@everywhere push!(LOAD_PATH, "src")
@everywhere import Metropolis

Metropolis.Runner.main()
