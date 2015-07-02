#include "Lattice.h"

#include <algorithm>

Lattice::Lattice(int N_) :
    N(N_),
    random_site(0, N - 1),
    sites(N*N)
{
    fill(1);
}

int Lattice::at(int i, int j) const {
    return sites[j + N*i];
}

void Lattice::fill(int s) {
    std::fill(sites.begin(), sites.end(), s);
}

const std::vector<int>& Lattice::get_sites() const {
    return sites;
}

const std::vector<int> Lattice::neighbors(int i, int j) const {
    std::vector<int> neigh;
    neigh.push_back(at(i, (j+1)%N));
    neigh.push_back(at(i, (N + j - 1)%N));
    neigh.push_back(at((i + 1)%N, j));
    neigh.push_back(at((N + i -1) % N, j));
    return neigh;
}

int Lattice::flip(int i, int j) {
    int s = at(i, j);
    sites[j + N*i] *= -s;
    return s;
}
