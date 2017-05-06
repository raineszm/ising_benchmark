module Runner
import ..Metropolis

const N = 64
const STEPS = 400
const T0 = 0.1
const TF = 5

function run_sim(N, c_in, c_out, LOCK)
    lat = Metropolis.Lattice(N)

    while isready(c_in)
        put!(LOCK, true)
        if !isready(LOCK)
            take!(LOCK)
            continue
        end

        t = take!(c_in)
        take!(LOCK)

        (M, U) = Metropolis.ensemble_av(lat, 1/t, 1000N^2, 100N^2)
        put!(c_out, (t, M, U))
    end
end

function main()
    T = linspace(T0, TF, STEPS)
    c_in = RemoteChannel(() -> Channel{Float64}(STEPS))
    c_out = RemoteChannel(() -> Channel{Tuple{Float64, Float64, Float64}}(STEPS))
    c_lock = RemoteChannel(() -> Channel{Bool}(1))

    for t in T
        put!(c_in, t)
    end
    close(c_in)

    futures = Array{Future}(nworkers())

    for i in 1:nworkers()
        futures[i] =  @spawn run_sim(N, c_in, c_out, c_lock)
    end

    out = open("met.dat", "w")
    @printf(out, "#T\tM\tU\n")

    for i in 1:STEPS
        (t, M, U) = take!(c_out)
        @printf(out, "%f\t%f\t%f\n", t, abs(M), U)
    end
    close(out)

    for f in futures
        wait(f)
    end
end


end
