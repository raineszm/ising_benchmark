using System;
using System.Linq;

namespace Ising
{
    public class Lattice
    {
        protected int[] sites;
        public int[] nearestNeighborsPlus;
        public int[] nearestNeighborsMinus;

        public readonly Random random = new Random();
        public readonly int size;

        public Lattice(int size)
        {
            this.size = size;
            sites = new int[size * size];
            Array.Fill(sites, 1);
            nearestNeighborsMinus = Enumerable.Range(0, size).Select(i => (i + size - 1) % size).ToArray();
            nearestNeighborsPlus = Enumerable.Range(0, size).Select(i => (i + 1) % size).ToArray();
        }


        public int this[int i, int j]
        {
            get => sites[size * j + i];
            set => sites[size * j + i] = value;
        }

        public int SumNeighbors(int i, int j)
        {
            return this[i, nearestNeighborsMinus[j]]
                + this[i, nearestNeighborsPlus[j]]
                + this[nearestNeighborsPlus[i], j]
                + this[nearestNeighborsMinus[i], j];

        }

        public int Flip(int i, int j) => -(this[i, j] *= -1);

        public int RandomSite() => random.Next() / (int.MaxValue / size + 1);

        public int Magnetization => sites.Sum();

        public int Energy
        {
            get
            {
                int total = 0;
                for (var i = 0; i < size; i++)
                {
                    for (var j = 0; j < size; j++)
                    {
                        total -= this[i, j] * (this[nearestNeighborsPlus[i], j] + this[i, nearestNeighborsMinus[j]]);
                    }
                }
                return total;

            }

        }


        public int EnergyChange(int i, int j) => 2 * this[i, j] * SumNeighbors(i, j);


    }
}