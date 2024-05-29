module Runner
using Distributed: RemoteChannel, RemoteException, workers, remote_do
using Printf: @printf
import ..Metropolis

const N = 64
const STEPS = 400
const T0 = 0.1
const TF = 5

function run_sim(N, c_in, c_out)
    lat = Metropolis.Lattice(N)

    println("Started worker")

    try
        while true
            t = take!(c_in)
            (M, U) = Metropolis.ensemble_av(lat, 1 / t, 1000, 100)
            put!(c_out, (t, M, U))
        end
    catch err
        if isa(err, RemoteException)
            captured = err.captured.ex
            if isa(captured, InvalidStateException)
                return
            end
        end
        throw(err)
    end
end

function feed_jobs(c_in, T)
    for t in T
        put!(c_in, t)
    end
    close(c_in)
end

function main(data_file)
    T = LinRange(T0, TF, STEPS)
    c_in = RemoteChannel(() -> Channel{Float64}(STEPS))
    c_out = RemoteChannel(() -> Channel{Tuple{Float64,Float64,Float64}}(STEPS))

    @async feed_jobs(c_in, T)

    for p in workers()
        remote_do(run_sim, p, N, c_in, c_out)
    end

    println("Starting")


    open(data_file, "w") do out
        @printf(out, "#T\tM\tU\n")
        flush(out)

        for i in 1:STEPS
            (t, M, U) = take!(c_out)

            ##Update us on the simulation progress
            if mod(t, 0.1) < 0.01
                println(t)
            end

            @printf(out, "%f\t%f\t%f\n", t, abs(M), U)
            flush(out)
        end
    end
    close(c_out)
end


end
