module Runner
using Base.Threads
using Printf: @printf
import ..Metropolis

const N = 128
const STEPS = 400
const T0 = 0.1
const TF = 5

function run_sim(N, c_in, c_out)
    lat = Metropolis.Lattice(N)

    while isready(c_in)
        t = take!(c_in)
        (M, U) = Metropolis.ensemble_av(lat, 1 / t, 1000, 100)
        put!(c_out, (t, M, U))
    end
end


function main(data_file)
    T = LinRange(T0, TF, STEPS)
    c_in = Channel{Float64}(STEPS) do ch
        for t in T
            put!(ch, t)
        end
    end
    c_out = Channel{Tuple{Float64,Float64,Float64}}(STEPS)

    for _ in 1:nthreads(:default)
        @spawn run_sim(N, c_in, c_out)
    end

    open(data_file, "w") do out
        @printf(out, "#T\tM\tU\n")
        flush(out)

        for i in 1:STEPS
            (t, M, U) = take!(c_out)

            ##Update us on the simulation progress
            if mod(t, 0.1) < 0.01
                println(t)
                flush(stdout)
            end

            @printf(out, "%f\t%f\t%f\n", t, abs(M), U)
            flush(out)
        end
    end
    close(c_out)
end


end
