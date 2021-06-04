module Runner
using Printf: @printf
import ..Metropolis

const N = 64
const STEPS = 400
const T0 = 0.1
const TF = 5

function main(data_file)
    T = LinRange(T0, TF, STEPS)
    lat = Metropolis.Lattice(N)

    open(data_file, "w") do out
        @printf(out, "#T\tM\tU\n")
        flush(out)

        for t in T
            (M, U) = Metropolis.ensemble_av(lat, 1/t, 1000, 100)

            ##Update us on the simulation progress
            if mod(t, 0.1) < 0.01
                println(t)
            end

            @printf(out, "%f\t%f\t%f\n", t, abs(M), U)
            flush(out)
        end
    end
end


end
