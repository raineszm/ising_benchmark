using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks.Dataflow;

namespace Ising
{
    public class Program
    {
        private const float T0 = 0.1f;
        private const float Tf = 5.0f;
        private const int N = 64;
        private const int Steps = 400;

        private readonly Lattice _lattice;
        private readonly Queue<int> _neighbors;

        private Program(int size)
        {
            _lattice = new Lattice(size);
            _neighbors = new Queue<int>(size * size);
        }

        private void PushNeighbors(int i)
        {
            _neighbors.Enqueue(_lattice.NearestNeighborsMinusX[i]);
            _neighbors.Enqueue(_lattice.NearestNeighborsMinusY[i]);
            _neighbors.Enqueue(_lattice.NearestNeighborsPlusX[i]);
            _neighbors.Enqueue(_lattice.NearestNeighborsPlusY[i]);
        }

        private (int, int) MetropolisStep(float beta)
        {
            var site = _lattice.RandomSite();

            var energyChange = _lattice.EnergyChange(site);

            var s = _lattice.Flip(site);
            var magnetizationChange = -2 * s;

            var flipProbability = 1 - MathF.Exp(-2 * beta);

            PushNeighbors(site);

            while (_neighbors.Count > 0)
            {
                site = _neighbors.Dequeue();

                if (_lattice[site] != s || !(_lattice.Random.NextDouble() < flipProbability)) continue;

                magnetizationChange -= 2 * s;
                energyChange += _lattice.EnergyChange(site);
                _lattice.Flip(site);

                PushNeighbors(site);
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
            var (energy, rmsMagnetization) = EnsembleAverage(1 / t, 1000, 100);
            return (t, energy, rmsMagnetization);
        }

        private static TransformBlock<float, (float, float, float)> WorkerBlock()
        {
            var program = new ThreadLocal<Program>(() => new Program(N));
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

            const float dt = (Tf - T0) / (Steps - 1);
            foreach (var t in Enumerable.Range(0, Steps).Select(i => T0 + i * dt)) input.Post(t);

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
                    var (t, energy, rmsMagnetization) = data;
                    // ReSharper disable once AccessToDisposedClosure
                    dataFile?.WriteLine($"{t},{energy},{rmsMagnetization}");
                });

            var worker = WorkerBlock();
            input.LinkTo(worker, linkOptions);
            worker.LinkTo(writer, linkOptions);

            input.Complete();
            writer.Completion.Wait();
        }
    }
}