#include <stddef.h>
#include <iostream>
#include <random>

#include "Lattice.h"

int deltaE(const Lattice& lat, int i, int j) {
    int dE = 0;
    for(auto neighbor : lat.neighbors(i, j)) {
        dE += neighbor;
    }
    return 2*dE*lat.at(i, j);
}

template <class URNG>
inline double randpp(URNG& g) {
    return std::generate_canonical<double, std::numeric_limits<double>::digits>(g);
}

template <class URNG>
bool metropolis_step(URNG& g,
        Lattice& lat,
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

template <class URNG>
void evolve(
        URNG& g,
        Lattice& lat,
        int n,
        double beta) {
    int dE, dM;

    for (auto i = 0; i < n; i++) {
        metropolis_step(g, lat, beta, dE, dM);
    }
}

int evaluate_energy(const Lattice& lat) {
    int total = 0;
    int N = lat.N;

    for (auto i = 0; i < N; i++) {
        for (auto j = 0; j < lat.N; j++) {
            total -= lat.at(i,j)*(
                    lat.at((i + 1) % N, j) +
                    lat.at(i, (j + 1) % N));
        }
    }

    return total;
}

template <class URNG>
void time_average(
        URNG& g,
        Lattice& lat,
        int n,
        double beta,
        double& en,
        double& mag) {
    int U = evaluate_energy(lat);
    const std::vector<int>& sites = lat.get_sites();
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

template <class URNG>
void ensemble_average(
        URNG& g,
        Lattice& lat,
        double beta,
        int n_evolve,
        int n_average,
        double& en,
        double& mag) {
    evolve(g, lat, n_evolve, beta);
    time_average(g, lat, n_average, beta, en, mag);
}



int main(int argc, char** argv) {

    int N = 64;
    Lattice lat(N);
    std::random_device rd;
    std::default_random_engine rng(rd());

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
