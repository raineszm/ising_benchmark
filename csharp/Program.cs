using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Threading.Tasks.Dataflow;

namespace Ising
{

    class IsingProgram
    {
        public const float T0 = 0.1f;
        public const float TF = 5.0f;
        public const int N = 64;
        public const int STEPS = 400;

        public Lattice lattice;
        public Queue<(int, int)> neighbors;

        public IsingProgram(int size)
        {
            lattice = new Lattice(size);
            neighbors = new Queue<(int, int)>(size * size);
        }

        public void PushNeighbors(int i, int j)
        {
            neighbors.Enqueue((i, lattice.nearestNeighborsPlus[j]));
            neighbors.Enqueue((i, lattice.nearestNeighborsMinus[j]));
            neighbors.Enqueue((lattice.nearestNeighborsMinus[i], j));
            neighbors.Enqueue((lattice.nearestNeighborsPlus[i], j));
        }

        public (int, int) MetropolisStep(float beta)
        {
            var i = lattice.RandomSite();
            var j = lattice.RandomSite();

            var energyChange = lattice.EnergyChange(i, j);

            var s = lattice.Flip(i, j);
            var magnetizationChange = -2 * s;

            var flipProbability = 1 - MathF.Exp(-2 * beta);

            PushNeighbors(i, j);

            while (neighbors.Count > 0)
            {
                (i, j) = neighbors.Dequeue();

                if (lattice[i, j] == s && lattice.random.NextDouble() < flipProbability)
                {
                    magnetizationChange -= 2 * s;
                    energyChange += lattice.EnergyChange(i, j);
                    lattice.Flip(i, j);

                    PushNeighbors(i, j);
                }

            }

            return (energyChange, magnetizationChange);

        }

        public void Evolve(int numberOfSteps, float beta)
        {
            for (int i = 0; i < numberOfSteps; i++)
            {
                MetropolisStep(beta);
            }
        }


        public (float, float) TimeAverage(int numberOfSteps, float beta)
        {
            var energy = lattice.Energy;
            var magnetization = lattice.Magnetization;

            var totalEnergy = 0;
            var totalMagnetization = 0;
            var totalMagnetizationSquared = 0;


            for (var i = 0; i < numberOfSteps; i++)
            {
                var (energyChange, magnetizationChange) = MetropolisStep(beta);

                energy += energyChange;
                magnetization += magnetizationChange;

                totalEnergy += energy;
                totalMagnetization += magnetization;
                totalMagnetizationSquared += magnetization * magnetization;

            }

            return ((float)totalEnergy / numberOfSteps, MathF.Sqrt((float)totalMagnetizationSquared / numberOfSteps));

        }

        public (float, float) EnsembleAverage(
            float beta,
            int numberOfEvolveSteps,
            int numberOfAverageSteps)
        {
            Evolve(numberOfEvolveSteps, beta);
            return TimeAverage(numberOfAverageSteps, beta);
        }

        public (float, float, float) Work(float t)
        {
            if (t % 0.1 < 0.01)
                Console.WriteLine($"From {Thread.CurrentThread.ManagedThreadId} T: {t}");
            var (E, M) = EnsembleAverage(1 / t, 1000, 100);
            return (t, E, M);
        }

        public static TransformBlock<float, (float, float, float)> WorkerBlock()
        {
            var program = new ThreadLocal<IsingProgram>(() => new IsingProgram(N));
            var options = new ExecutionDataflowBlockOptions
            {
                MaxDegreeOfParallelism = Environment.ProcessorCount
            };
            return new TransformBlock<float, (float, float, float)>(
                t => program.Value.Work(t), options);
        }

        static void Main(string[] args)
        {
            var input = new BufferBlock<float>();

            float dt = (TF - T0) / (STEPS - 1);
            foreach (var t in Enumerable.Range(0, STEPS).Select(i => T0 + i * dt))
            {
                input.Post(t);
            }

            string filePath = args.Length > 1 ? args[1] : @"data.csv";

            using var dataFile = new StreamWriter(filePath);

            var linkOptions = new DataflowLinkOptions
            {
                PropagateCompletion = true
            };

            dataFile.WriteLine("T,U,M_rms");

            var writer = new ActionBlock<(float, float, float)>(
                data =>
                {
                    var (t, E, M) = data;
                    dataFile.WriteLine($"{t},{E},{M}");
                });

            var worker = WorkerBlock();
            input.LinkTo(worker, linkOptions);
            worker.LinkTo(writer, linkOptions);

            input.Complete();
            writer.Completion.Wait();
        }

    }
}
