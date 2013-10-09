#!/usr/bin/env julia

module Metropolis

N = 64 # Lattice Size
#Holds the nearest neighbor indices for each index
global nn1 
global nn2


#Calculate the energy change due to a spin flip
function delta_E(s :: Array{Int8, 2}, i, j)
    2*s[i,j]*(s[nn1[i],j]
        +s[nn2[i],j]
        +s[i,nn1[j]]
        +s[i,nn2[j]])
end

#Determine whether to accept a proposed new configuration
accepted(dE, beta) = rand() < exp(-beta*dE)

function metropolis_step!(s :: Array{Int8, 2}, beta :: Real)
    (i, j) = rand(1:(N::Int), 2)

    dE = delta_E(s, i, j)

    if dE < 0 || accepted(dE, beta)
        spin = s[i,j]
        s[i,j] = - s[i,j]
        return dE, -2spin
    end
    0,0
end

function evolve(s :: Array{Int8,2}, n, beta :: Real)
    for i=1:n
        metropolis_step!(s, beta)
    end
end

function time_average(s :: Array{Int8, 2}, n, beta :: Real)
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

function evaluate_energy(s :: Array{Int8, 2})
    total = 0

    for i=1:(N::Int)
        for j=1:(N::Int)
            total -= s[i,j]*(s[nn1[i],j]+ s[i,nn2[j]])
        end
    end
    total
end

function ensemble_av(beta, n_evolve, n_average)
    srand(ifloor(time()))
    spins = ones(Int8, N, N)

    T = 1/beta
    #Update us on the simulation progress
    if mod(T, 0.1) < 0.01
        println(T)
    end

    evolve(spins, n_evolve, beta)
    time_average(spins, n_average, beta)
end

function set_lattice_size(n :: Integer)
    global N
    global nn1, nn2
    
    N = n

    r = 1:N
    nn1 = [r[end], r[1:(end-1)]]
    nn2 = [r[2:end], r[1]]

    N
end

end
