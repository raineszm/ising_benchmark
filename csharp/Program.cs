using System;
using System.Collections.Generic;

namespace Ising
{

    class Program
    {
        public const double T0 = 0.1;
        public const double TF = 5.0;
        public const int N = 64;
        public const int STEPS = 400;
        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
        }


        // const int NUM_THREADS = static_cast<int>(std::thread::hardware_concurrency());

        // inline double
        // rand_double(std::minstd_rand& rng)
        // {
        //   return std::generate_canonical<double, 8>(rng);
        // }

        public static void PushNeighbors(Lattice lattice, int i, int j, Queue<(int, int)> queue)
        {
            queue.Enqueue((i, lattice.nearestNeighborsPlus[j]));
            queue.Enqueue((i, lattice.nearestNeighborsMinus[j]));
            queue.Enqueue((lattice.nearestNeighborsMinus[i], j));
            queue.Enqueue((lattice.nearestNeighborsPlus[i], j));
        }

        public static (int, int) MetropolisStep(Lattice lattice, double beta)
        {
            var i = lattice.RandomSite();
            var j = lattice.RandomSite();
            var neighbors = new Queue<(int, int)>();

            var energyChange = lattice.EnergyChange(i, j);

            var s = lattice.Flip(i, j);
            var magnetizationChange = 0;

            var flipProbability = 1 - Math.Exp(-2 * beta);

            PushNeighbors(lattice, i, j, neighbors);

            while (neighbors.Count > 0)
            {
                (i, j) = neighbors.Dequeue();

                if (lattice[i, j] == s && lattice.random.NextDouble() < flipProbability)
                {
                    magnetizationChange -= 2 * s;
                    energyChange += lattice.EnergyChange(i, j);
                    lattice.Flip(i, j);

                    PushNeighbors(lattice, i, j, neighbors);
                }

            }

            return (energyChange, magnetizationChange);

        }

        public static void Evolve(Lattice lattice, int numberOfSteps, double beta)
        {
            for (int i = 0; i < numberOfSteps; i++)
            {
                MetropolisStep(lattice, beta);
            }
        }

        // inline double
        // average(long agg, int n)
        // {
        //   return static_cast<double>(agg) / n;
        // }

        public static (double, double) TimeAverage(Lattice lattice, int numberOfSteps, double beta)
        {
            var energy = lattice.Energy;
            var magnetization = lattice.Magnetization;

            var totalEnergy = 0;
            var totalMagnetization = 0;
            var totalMagnetizationSquared = 0;


            for (var i = 0; i < numberOfSteps; i++)
            {
                var (energyChange, magnetizationChange) = MetropolisStep(lattice, beta);

                energy += energyChange;
                magnetization += magnetizationChange;

                totalEnergy += energy;
                totalMagnetization += magnetization;
                totalMagnetizationSquared += magnetization * magnetization;

            }

            return ((double)totalEnergy / numberOfSteps, Math.Sqrt((double)totalMagnetizationSquared / numberOfSteps));

        }

        public static (double, double) EnsembleAverage(Lattice lattice,
                                                       double beta,
                                                       int numberOfEvolveSteps,
                                                       int numberOfAverageSteps)
        {
            Evolve(lattice, numberOfEvolveSteps, beta);
            return TimeAverage(lattice, numberOfAverageSteps, beta);
        }


        // void
        // worker(std::queue<double>& ts, std::mutex& mtx, Channel<data_type>& chan)
        // {
        //   Lattice<N> lat;

        //   while (true) {
        //     double t;

        //     {
        //       std::unique_lock<std::mutex> lck(mtx);

        //       if (ts.empty()) {

        //         break;
        //       }

        //       t = ts.front();
        //       ts.pop();

        //       if (fmod(t, 0.1) < 0.01) {
        //         std::cout << "T: " << t << std::endl;
        //       }
        //     }

        //     auto observables = ensemble_average(lat, 1. / t, 1000, 100);

        //     chan.put(std::tuple_cat(std::make_tuple(t), observables));
        //   }
        // }

        // int
        // main(int argc, char** argv)
        // {
        //   Lattice<N> lat;

        //   std::queue<double> ts;
        //   Channel<data_type> chan;

        //   // Build vector of temperatures
        //   double dt = (TF - T0) / (STEPS - 1);
        //   double t = T0;
        //   for (auto i = 0; i < STEPS; i++) {
        //     ts.push(t);
        //     t += dt;
        //   }

        //   std::vector<std::thread> threads;
        //   std::mutex mtx;

        //   for (auto i = 0; i < NUM_THREADS; i++) {
        //     threads.push_back(
        //       std::thread([&ts, &mtx, &chan]() { worker(ts, mtx, chan); }));
        //   }

        //   std::string file_path = "data.csv";
        //   if (argc > 1) {
        //     file_path = argv[1];

        //   }
        //   std::ofstream data_file(file_path);

        //   if (!data_file.is_open()) {
        //     std::cerr << "Could not open data file" << std::endl;
        //     return -1;
        //   } else {
        //     data_file << "T,U,M_rms" << std::endl;
        //     for (int i = 0; i < STEPS; i++) {
        //       data_type row = chan.take();
        //       data_file << std::get<0>(row) << ", " << std::get<1>(row) << ", "
        //                 << std::get<2>(row) << std::endl;
        //     }
        //   }

        //   for (std::thread& t : threads)
        //     t.join();

        //   return 0;
        // }

    }
}
