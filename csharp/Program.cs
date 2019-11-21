using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace Ising
{

    class Program
    {
        public const float T0 = 0.1f;
        public const float TF = 5.0f;
        public const int N = 64;
        public const int STEPS = 400;
        public static void PushNeighbors(Lattice lattice, int i, int j, Queue<(int, int)> queue)
        {
            queue.Enqueue((i, lattice.nearestNeighborsPlus[j]));
            queue.Enqueue((i, lattice.nearestNeighborsMinus[j]));
            queue.Enqueue((lattice.nearestNeighborsMinus[i], j));
            queue.Enqueue((lattice.nearestNeighborsPlus[i], j));
        }

        public static (int, int) MetropolisStep(Lattice lattice, float beta)
        {
            var i = lattice.RandomSite();
            var j = lattice.RandomSite();
            var neighbors = new Queue<(int, int)>();

            var energyChange = lattice.EnergyChange(i, j);

            var s = lattice.Flip(i, j);
            var magnetizationChange = 0;

            var flipProbability = 1 - MathF.Exp(-2 * beta);

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

        public static void Evolve(Lattice lattice, int numberOfSteps, float beta)
        {
            for (int i = 0; i < numberOfSteps; i++)
            {
                MetropolisStep(lattice, beta);
            }
        }

        // inline float
        // average(long agg, int n)
        // {
        //   return static_cast<float>(agg) / n;
        // }

        public static (float, float) TimeAverage(Lattice lattice, int numberOfSteps, float beta)
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

            return ((float)totalEnergy / numberOfSteps, MathF.Sqrt((float)totalMagnetizationSquared / numberOfSteps));

        }

        public static (float, float) EnsembleAverage(Lattice lattice,
                                                       float beta,
                                                       int numberOfEvolveSteps,
                                                       int numberOfAverageSteps)
        {
            Evolve(lattice, numberOfEvolveSteps, beta);
            return TimeAverage(lattice, numberOfAverageSteps, beta);
        }

        public static void Worker(BlockingCollection<float> ts, BlockingCollection<(float, float, float)> output)
        {
            Lattice lattice = new Lattice(N);

            while (ts.TryTake(out var t))
            {
                var (E, M) = EnsembleAverage(lattice, 1 / t, 1000, 100);

                if (t % 0.1 < 0.01)
                    Console.WriteLine($"T: {t}");


                output.Add((t, E, M));
            }
        }

        static void Main(string[] args)
        {
            var ts = new BlockingCollection<float>(new ConcurrentQueue<float>());
            var output = new BlockingCollection<(float, float, float)>(new ConcurrentQueue<(float, float, float)>());


            float dt = (TF - T0) / (STEPS - 1);
            foreach (var t in Enumerable.Range(0, STEPS).Select(i => T0 + i * dt))
            {
                ts.Add(t);
            }
            ts.CompleteAdding();

            Task[] tasks =
                Enumerable.Range(0, Environment.ProcessorCount - 1).Select(
                _ => Task.Factory.StartNew(() => Worker(ts, output),
                        creationOptions: TaskCreationOptions.LongRunning)
                ).ToArray();

            string filePath = args.Length > 1 ? args[1] : @"data.csv";

            using (var dataFile = new StreamWriter(filePath))
            {
                dataFile.WriteLine("T,U,M_rms");

                for (var i = 0; i < STEPS; i++)
                {
                    var (t, E, M) = output.Take();
                    dataFile.WriteLine($"{t},{E},{M}");

                }
            }

            foreach (var task in tasks)
                task.Wait();

        }



    }
}
