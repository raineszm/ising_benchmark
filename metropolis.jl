#!/usr/bin/env julia

module Metropolis

N = 64 # Lattice Size
#Holds the nearest neighbor indices for each index
global nn1 
global nn2
global magnetization = 0
global energy = 0

#The lattices
#Spins are 1 or -1
spins = ones(Int8, N, N)


#---------------- RANDOM NUMBER GENERATION ----------------------------
#Generate a random integer less than n
rand_int(n :: Int) = rand(Uint) % n + 1
#------------------------------- END ----------------------------------


#Calculate the energy change due to a spin flip
function delta_E(s :: Array{Int8, 2}, i, j)
    2*s[i,j]*(s[nn1[i],j]+s[nn2[i],j]
        +s[i,nn1[j]]+s[i,nn2[j]])
end

#Determine whether to accept a proposed new configuration
accepted(dE :: Real, beta :: Real) = rand() < exp(-beta*dE)

function metropolis_step(s :: Array{Int8, 2}, beta :: Real)
    global energy, magnetization
    i = rand_int(N)
    j = rand_int(N)

    dE = delta_E(s, i, j)

    if dE < 0 || accepted(dE, beta)
        spin = s[i,j]
        s[i,j] = - s[i,j]
        energy += dE
        magnetization -= 2spin
    end
    s
end

function evolve(s :: Array{Int8,2}, n, beta :: Real)
    for i=1:n
        s=metropolis_step(s, beta)
    end
end

function time_average(s :: Array{Int8, 2}, n, beta :: Real)
    mag = 0
    en = 0

    for i=1:N
        metropolis_step(s, beta)
        mag += magnetization
        en += energy
    end

    return mag/n, en/n
end

function evaluate_energy(s :: Array{Int8, 2})
    total = 0

    for i=1:N
        for j=1:N
            total -= s[i,j]*(s[nn1[i],j]+ s[i,nn2[j]])
        end
    end
    total
end

function ensemble_av(beta, n_evolve :: Integer, n_average :: Integer)
    T = 1/beta
    #Update us on the simulation progress
    if mod(T, 0.1) < 0.02
        println(T)
    end

    evolve(spins, n_evolve, beta)
    time_average(spins, n_average, beta)
end

function initialize(n :: Integer)
    global N
    global spins
    global magnetization
    global energy

    N = n

    srand(ifloor(time()))

    spins = ones(Int8, N, N)

    populate_nn()

    magnetization = sum(spins)
    energy = evaluate_energy(spins)
    true
end



function populate_nn()
    r = 1:N
    global nn1, nn2
    nn1 = [r[end], r[1:(end-1)]]
    nn2 = [r[2:end], r[1]]
    true
end
end
