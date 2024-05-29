#!/usr/bin/env julia

module Metropolis


mutable struct Lattice
    N::Int
    nn_plus::Vector{Int}
    nn_minus::Vector{Int}
    spins::Array{Int,2}
    queue::Vector{Tuple{Int,Int}}
end

function Lattice(N::Int)
    r = 1:N
    nn_plus = [r[end], r[1:(end-1)]...]
    nn_minus = [r[2:end]..., r[1]]
    spins = ones(Int, N, N)
    queue = Vector{Tuple{Int,Int}}()
    sizehint!(queue, N * N)
    Lattice(N, nn_plus, nn_minus, spins, queue)
end

#Calculate the energy change due to a spin flip
function delta_E(lattice::Lattice,
    i::Int, j::Int)
    2 * lattice.spins[i, j] * (lattice.spins[lattice.nn_plus[i], j]
                               + lattice.spins[lattice.nn_minus[i], j]
                               + lattice.spins[i, lattice.nn_plus[j]]
                               + lattice.spins[i, lattice.nn_minus[j]])
end

function flip!(lattice::Lattice, i::Int, j::Int)
    @inbounds spin = lattice.spins[i, j]
    @inbounds lattice.spins[i, j] = -spin
    spin
end

function push_neighbors!(lattice::Lattice, i::Int, j::Int, queue)
    push!(queue, (i, lattice.nn_plus[j]))
    push!(queue, (i, lattice.nn_minus[j]))
    push!(queue, (lattice.nn_plus[i], j))
    push!(queue, (lattice.nn_minus[i], j))
end


function metropolis_step!(lattice::Lattice, beta::Float64)
    i = rand(1:lattice.N)
    j = rand(1:lattice.N)


    empty!(lattice.queue)

    # Flip the picked spin
    spin = flip!(lattice, i, j)

    dE = delta_E(lattice, i, j)
    dM = -2 * spin

    flip_prob = 1 - exp(-2 * beta)

    push_neighbors!(lattice, i, j, lattice.queue)

    while !isempty(lattice.queue)
        i, j = popfirst!(lattice.queue)
        if lattice.spins[i, j] == spin && rand() < flip_prob
            dM -= 2 * spin
            dE += delta_E(lattice, i, j)

            flip!(lattice, i, j)

            push_neighbors!(lattice, i, j, lattice.queue)

        end
    end

    dE, dM
end

function evolve!(lattice::Lattice, n::Int, beta::Float64)
    for i = 1:n
        metropolis_step!(lattice, beta)
    end
end

function time_average!(lattice::Lattice, n::Int, beta::Float64)
    U = evaluate_energy(lattice)
    M = sum(lattice.spins)
    mag = 0
    en = 0

    for i = 1:n
        (dE, dM) = metropolis_step!(lattice, beta)
        U += dE
        M += dM
        mag += M * M
        en += U
    end

    return sqrt(mag / n), en / n
end

function evaluate_energy(lattice::Lattice)
    total = 0
    N = lattice.N

    for i = 1:N
        @simd for j = 1:N
            @inbounds total -=
                lattice.spins[i, j] * (
                    lattice.spins[lattice.nn_plus[i], j]
                    +
                    lattice.spins[i, lattice.nn_minus[j]])
        end
    end
    total
end

function ensemble_av(lattice, beta, n_evolve, n_average)
    T = 1 / beta

    evolve!(lattice, n_evolve, beta)
    time_average!(lattice, n_average, beta)
end

include("Runner.jl")

end
