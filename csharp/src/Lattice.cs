using System;
using System.Linq;

namespace Ising
{
    public class Lattice
    {
        private readonly int[] _sites;
        private readonly int _size;
        public readonly int[] NearestNeighborsMinus;
        public readonly int[] NearestNeighborsPlus;

        public readonly Random Random = new Random();

        public Lattice(int size)
        {
            _size = size;
            _sites = new int[size * size];
            Array.Fill(_sites, 1);
            NearestNeighborsMinus = Enumerable.Range(0, size).Select(i => (i + size - 1) % size).ToArray();
            NearestNeighborsPlus = Enumerable.Range(0, size).Select(i => (i + 1) % size).ToArray();
        }


        public int this[int i, int j]
        {
            get => _sites[_size * j + i];
            private set => _sites[_size * j + i] = value;
        }

        public int Magnetization => _sites.Sum();

        public int Energy
        {
            get
            {
                var total = 0;
                for (var i = 0; i < _size; i++)
                for (var j = 0; j < _size; j++)
                    total -= this[i, j] * (this[NearestNeighborsPlus[i], j] + this[i, NearestNeighborsMinus[j]]);
                return total;
            }
        }

        private int SumNeighbors(int i, int j)
        {
            return this[i, NearestNeighborsMinus[j]]
                   + this[i, NearestNeighborsPlus[j]]
                   + this[NearestNeighborsPlus[i], j]
                   + this[NearestNeighborsMinus[i], j];
        }

        public int Flip(int i, int j) => -(this[i, j] *= -1);

        public int RandomSite() => Random.Next() / (int.MaxValue / _size + 1);

        public int EnergyChange(int i, int j) => 2 * this[i, j] * SumNeighbors(i, j);
    }
}