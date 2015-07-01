#include "Lattice.h"

#include <algorithm>

Lattice::Lattice(int N_) : N(N_), sites(N*N) {
    fill(1);
}

int Lattice::at(int i, int j) const {
    return sites[i + N*j];
}

void Lattice::fill(int s) {
    std::fill(sites.begin(), sites.end(), s);
}

const std::vector<int>& Lattice::get_sites() const {
    return sites;
}
