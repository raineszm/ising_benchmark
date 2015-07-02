#include <stddef.h>
#include <iostream>
#include <array>

#include "Lattice.h"

template <int N>
inline int deltaE(const Lattice<N>& lat, int i, int j) {
    return 2*lat.at(i, j)*lat.sum_neighbors(i, j);
}

template <class URNG>
inline double randpp(URNG& g) {
    return std::generate_canonical<double, std::numeric_limits<double>::digits>(g);
}

template <class URNG, int N>
bool metropolis_step(URNG& g,
        Lattice<N>& lat,
        double beta,
        int& dE,
        int& dM) {
    int i = lat.random_site(g);
    int j = lat.random_site(g);

    dE = deltaE(lat, i, j);

    if (dE < 0 || randpp(g) < exp(-beta*dE)) {
        dM = -2*lat.flip(i, j);
        return true;
    } else {
        dE = 0;
        dM = 0;
        return false;
    }
}

template <class URNG, int N>
void evolve(
        URNG& g,
        Lattice<N>& lat,
        int n,
        double beta) {
    int dE, dM;

    for (auto i = 0; i < n; i++) {
        metropolis_step(g, lat, beta, dE, dM);
    }
}

template <int N>
int evaluate_energy(const Lattice<N>& lat) {
    int total = 0;

    for (auto i = 0; i < N; i++) {
        for (auto j = 0; j < N; j++) {
            total -= lat.at(i,j)*(
                    lat.at((i + 1) % N, j) +
                    lat.at(i, (j + 1) % N));
        }
    }

    return total;
}

template <class URNG, int N>
void time_average(
        URNG& g,
        Lattice<N>& lat,
        int n,
        double beta,
        double& en,
        double& mag) {
    int U = evaluate_energy(lat);
    const std::array<int, N*N>& sites = lat.get_sites();
    int M = std::accumulate(sites.begin(), sites.end(), 0);
    int U_tot = 0;
    int M_tot = 0;
    int dE, dM;

    for (auto i = 0; i < n; i++) {
        metropolis_step(g, lat, beta, dE, dM);
        U += dE;
        M += dM;
        U_tot += U;
        M_tot += M;
    }

    en = (float) U/n;
    mag = (float) M/n;
}

template <class URNG, int N>
void ensemble_average(
        URNG& g,
        Lattice<N>& lat,
        double beta,
        int n_evolve,
        int n_average,
        double& en,
        double& mag) {
    evolve(g, lat, n_evolve, beta);
    time_average(g, lat, n_average, beta, en, mag);
}



int main(int, char**) {

    const int N = 64;
    Lattice<N> lat;
    std::random_device rd;
    std::minstd_rand rng(rd());

    double en, mag;

    double dt = (5 - 0.1)/400;
    double t = 0.1;
    for (int i = 0; i < 400; i ++) {
        ensemble_average(rng, lat, 1/t, 1000*N*N, 100*N*N, en, mag);

        if (fmod(t, 0.1) < 0.01) {
            std::cout << "T: " << t << std::endl;
            std::cout << "U: " << en << std::endl;
            std::cout << "M: " << mag<< std::endl;
        }
        t += dt;
    }
}
