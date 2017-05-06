#!/usr/bin/env julia

module Metropolis
    type Lattice
        N :: Int
        nn_plus :: Vector{Int}
        nn_minus :: Vector{Int}
        spins :: Array{Int, 2}
    end

    function Lattice(N :: Int)
        r = 1:N
        nn_plus = [r[end], r[1:(end-1)]...]
        nn_minus = [r[2:end]..., r[1]]
        spins = ones(Int, N, N)
        Lattice(N, nn_plus, nn_minus, spins)
    end

    crand(n :: Int) = ccall(:rand, Cint, ()) % n + 1

    #Calculate the energy change due to a spin flip
    function delta_E(lattice :: Lattice,
                     i :: Int, j :: Int)
        2*lattice.spins[i,j]*(lattice.spins[lattice.nn_plus[i],j]
            +lattice.spins[lattice.nn_minus[i],j]
            +lattice.spins[i, lattice.nn_plus[j]]
            +lattice.spins[i,lattice.nn_minus[j]])
    end

    function metropolis_step!(lattice :: Lattice, beta :: Float64)
        i = crand(lattice.N)
        j = crand(lattice.N)

        dE = delta_E(lattice, i, j)

        if dE < 0 || rand() < ccall(:exp, Cdouble, (Cdouble,), -beta*dE)
            @inbounds spin = lattice.spins[i,j]
            @inbounds lattice.spins[i,j] = -spin
            return dE, -2spin
        end
        0,0
    end

    function evolve!(lattice :: Lattice, n :: Int, beta :: Float64)
        for i=1:n
            metropolis_step!(lattice, beta)
        end
    end

    function time_average!(lattice :: Lattice, n :: Int, beta :: Float64)
        U = evaluate_energy(lattice)
        M = sum(lattice.spins)
        mag = 0
        en = 0

        for i=1:n
            (dE, dM) = metropolis_step!(lattice, beta)
            U += dE
            M += dM
            mag += M
            en += U
        end

        return mag/n, en/n
    end

    function evaluate_energy(lattice :: Lattice)
        total = 0
        N = lattice.N

        for i=1:N
            @simd for j=1:N
                @inbounds total -=
                    lattice.spins[i,j]*(
                                    lattice.spins[lattice.nn_plus[i],j]
                                    + lattice.spins[i,lattice.nn_minus[j]])
            end
        end
        total
    end

    function ensemble_av(lattice, beta, n_evolve, n_average)
        srand(floor(Int, time()))

        T = 1/beta
        #Update us on the simulation progress
        if mod(T, 0.1) < 0.01
            println(T)
        end

        evolve!(lattice, n_evolve, beta)
        time_average!(lattice, n_average, beta)
    end

    include("Runner.jl")


end
