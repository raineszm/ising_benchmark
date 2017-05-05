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
        nn_plus = [r[end], r[1:(end-1)]]
        nn_minus = [r[2:end], r[1]]
        spins = ones(Int, N, N)
        Lattice(N, nn_plus, nn_minus, spins)
    end

    const lattice = Lattice(64)

    crand(n :: Int) = ccall(:rand, Cint, ()) % n + 1


    #Calculate the energy change due to a spin flip
    function delta_E(s :: Array{Int, 2}, i :: Int, j :: Int)
        2*s[i,j]*(s[lattice.nn_plus[i],j]
            +s[lattice.nn_minus[i],j]
            +s[i, lattice.nn_plus[j]]
            +s[i,lattice.nn_minus[j]])
    end

    function metropolis_step!(s :: Array{Int, 2}, beta :: Float64)
        i = crand(lattice.N)
        j = crand(lattice.N)

        dE = delta_E(s, i, j)

        if dE < 0 || rand() < ccall(:exp, Cdouble, (Cdouble,), -beta*dE)
            @inbounds spin = s[i,j]
            @inbounds s[i,j] = -spin
            return dE, -2spin
        end
        0,0
    end

    function evolve!(s :: Array{Int,2}, n :: Int, beta :: Float64)
        for i=1:n
            metropolis_step!(s, beta)
        end
    end

    function time_average!(s :: Array{Int, 2}, n :: Int, beta :: Float64)
        U = evaluate_energy(s)
        M = sum(s)
        mag = 0
        en = 0

        for i=1:n
            (dE, dM) = metropolis_step!(s, beta)
            U += dE
            M += dM
            mag += M
            en += U
        end

        return mag/n, en/n
    end

    function evaluate_energy(s :: Array{Int, 2})
        total = 0
        N = lattice.N

        for i=1:N
            @simd for j=1:N
                @inbounds total -= s[i,j]*(s[lattice.nn_plus[i],j]
                + s[i,lattice.nn_minus[j]])
            end
        end
        total
    end

    function ensemble_av(beta, n_evolve, n_average)
        srand(ifloor(time()))
        spins = lattice.spins

        T = 1/beta
        #Update us on the simulation progress
        if mod(T, 0.1) < 0.01
            println(T)
        end

        evolve!(spins, n_evolve, beta)
        time_average!(spins, n_average, beta)
    end


    function init_lattice!(n :: Int)
        l = Lattice(n)
        global lattice
        lattice.N = l.N
        lattice.nn_plus = l.nn_plus
        lattice.nn_minus = l.nn_minus
        lattice.spins = l.spins
        true
    end

end
