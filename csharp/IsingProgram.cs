using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks.Dataflow;

namespace Ising
{
    public class IsingProgram
    {
        private const float T0 = 0.1f;
        private const float TF = 5.0f;
        private const int N = 64;
        private const int STEPS = 400;

        private readonly Lattice _lattice;
        private readonly Queue<(int, int)> _neighbors;

        private IsingProgram(int size)
        {
            _lattice = new Lattice(size);
            _neighbors = new Queue<(int, int)>(size * size);
        }

        private void PushNeighbors(int i, int j)
        {
            _neighbors.Enqueue((i, _lattice.NearestNeighborsPlus[j]));
            _neighbors.Enqueue((i, _lattice.NearestNeighborsMinus[j]));
            _neighbors.Enqueue((_lattice.NearestNeighborsMinus[i], j));
            _neighbors.Enqueue((_lattice.NearestNeighborsPlus[i], j));
        }

        private (int, int) MetropolisStep(float beta)
        {
            var i = _lattice.RandomSite();
            var j = _lattice.RandomSite();

            var energyChange = _lattice.EnergyChange(i, j);

            var s = _lattice.Flip(i, j);
            var magnetizationChange = -2 * s;

            var flipProbability = 1 - MathF.Exp(-2 * beta);

            PushNeighbors(i, j);

            while (_neighbors.Count > 0)
            {
                (i, j) = _neighbors.Dequeue();

                if (_lattice[i, j] == s && _lattice.Random.NextDouble() < flipProbability)
                {
                    magnetizationChange -= 2 * s;
                    energyChange += _lattice.EnergyChange(i, j);
                    _lattice.Flip(i, j);

                    PushNeighbors(i, j);
                }
            }

            return (energyChange, magnetizationChange);
        }

        private void Evolve(int numberOfSteps, float beta)
        {
            for (var i = 0; i < numberOfSteps; i++) MetropolisStep(beta);
        }


        private (float, float) TimeAverage(int numberOfSteps, float beta)
        {
            var energy = _lattice.Energy;
            var magnetization = _lattice.Magnetization;

            var totalEnergy = 0;
            var totalMagnetizationSquared = 0;


            for (var i = 0; i < numberOfSteps; i++)
            {
                var (energyChange, magnetizationChange) = MetropolisStep(beta);

                energy += energyChange;
                magnetization += magnetizationChange;

                totalEnergy += energy;
                totalMagnetizationSquared += magnetization * magnetization;
            }

            return ((float) totalEnergy / numberOfSteps, MathF.Sqrt((float) totalMagnetizationSquared / numberOfSteps));
        }

        private (float, float) EnsembleAverage(
            float beta,
            int numberOfEvolveSteps,
            int numberOfAverageSteps)
        {
            Evolve(numberOfEvolveSteps, beta);
            return TimeAverage(numberOfAverageSteps, beta);
        }

        private (float, float, float) Work(float t)
        {
            if (t % 0.1 < 0.01)
                Console.WriteLine($"From {Thread.CurrentThread.ManagedThreadId} T: {t}");
            var (E, M) = EnsembleAverage(1 / t, 1000, 100);
            return (t, E, M);
        }

        private static TransformBlock<float, (float, float, float)> WorkerBlock()
        {
            var program = new ThreadLocal<IsingProgram>(() => new IsingProgram(N));
            var options = new ExecutionDataflowBlockOptions
            {
                MaxDegreeOfParallelism = Environment.ProcessorCount
            };
            return new TransformBlock<float, (float, float, float)>(
                t => program.Value.Work(t), options);
        }

        private static void Main(string[] args)
        {
            var input = new BufferBlock<float>();

            const float dt = (TF - T0) / (STEPS - 1);
            foreach (var t in Enumerable.Range(0, STEPS).Select(i => T0 + i * dt)) input.Post(t);

            var filePath = args.Length > 1 ? args[1] : @"data.csv";

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
                    dataFile?.WriteLine($"{t},{E},{M}");
                });

            var worker = WorkerBlock();
            input.LinkTo(worker, linkOptions);
            worker.LinkTo(writer, linkOptions);

            input.Complete();
            writer.Completion.Wait();
        }
    }
}