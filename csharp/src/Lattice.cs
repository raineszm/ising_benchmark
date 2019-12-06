using System;
using System.Collections.Generic;
using System.Linq;

namespace Ising
{
    public class Lattice
    {
        private readonly int[] _sites;
        private readonly int _size;
        public readonly int[] NearestNeighborsMinusX;
        public readonly int[] NearestNeighborsPlusX;
        public readonly int[] NearestNeighborsMinusY;
        public readonly int[] NearestNeighborsPlusY;

        public readonly Random Random = new Random();

        public Lattice(int size)
        {
            _size = size;
            _sites = new int[size * size];
            Array.Fill(_sites, 1);
            var indices = Enumerable.Range(0, size).ToArray();
            var nearestNeighborsMinus = indices.Select(i => (i + size - 1) % size).ToArray();
            var nearestNeighborsPlus = indices.Select(i => (i + 1) % size).ToArray();
            NearestNeighborsMinusX = LinearProduct(nearestNeighborsMinus, indices);
            NearestNeighborsPlusX = LinearProduct(nearestNeighborsPlus, indices);
            NearestNeighborsMinusY = LinearProduct(indices, nearestNeighborsMinus);
            NearestNeighborsPlusY = LinearProduct(indices, nearestNeighborsPlus);
        }

        private int[] LinearProduct(IEnumerable<int> xs, IEnumerable<int> ys)
        {
            return (from j in ys
                    from i in xs
                    select LinearIndex(i, j)
                ).ToArray();
        }

        private int LinearIndex(int i, int j) => _size * j + i;

        public int this[int i]
        {
            get => _sites[i];
            private set => _sites[i] = value;
        }

        public int Magnetization => _sites.Sum();

        public int Energy
        {
            get
            {
                var total = 0;
                for (var i = 0; i < _size * _size; i++)
                    total -= this[i] * (this[NearestNeighborsPlusX[i]] + this[NearestNeighborsMinusY[i]]);
                return total;
            }
        }

        private int SumNeighbors(int i)
        {
            return this[NearestNeighborsMinusX[i]]
                   + this[NearestNeighborsPlusX[i]]
                   + this[NearestNeighborsPlusY[i]]
                   + this[NearestNeighborsMinusY[i]];
        }

        public int Flip(int i) => -(this[i] *= -1);

        public int RandomSite() => Random.Next() / (int.MaxValue / (_size * _size) + 1);

        public int EnergyChange(int i) => 2 * this[i] * SumNeighbors(i);
    }
}