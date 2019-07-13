#!/usr/bin/env julia

module Metropolis
using DataStructures: Deque

mutable struct Lattice
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

function flip!(lattice :: Lattice, i:: Int, j:: Int)
    @inbounds spin = lattice.spins[i,j]
    @inbounds lattice.spins[i,j] = -spin
    spin
end

function push_neighbors!(lattice :: Lattice, i :: Int, j :: Int, queue)
    push!(queue, (i, lattice.nn_plus[i]))
    push!(queue, (i, lattice.nn_minus[i]))
    push!(queue, (lattice.nn_plus[j], j))
    push!(queue, (lattice.nn_minus[j], j))
end


function metropolis_step!(lattice :: Lattice, beta :: Float64)
    i = crand(lattice.N)
    j = crand(lattice.N)

    neighbors = Deque{Tuple{Int, Int}}()

    # Flip the picked spin
    spin = flip!(lattice, i, j)

    dE = delta_E(lattice, i, j)
    dM = -2*spin

    flip_prob = 1 - exp(-2*beta)

    push_neighbors!(lattice, i, j, neighbors)

    while !isempty(neighbors)
        i, j = popfirst!(neighbors)
        if lattice.spins[i, j] == spin && rand() < flip_prob
            dM -= 2*spin
            dE += delta_E(lattice, i, j)

            flip!(lattice, i, j)

            push_neighbors!(lattice, i, j, neighbors)

        end
    end

    dE, dM
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
        mag += M * M
        en += U
    end

    return sqrt(mag/n), en/n
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
