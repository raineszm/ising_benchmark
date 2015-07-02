#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <array>
#include <tuple>

#include <random>

#include <fstream>
#include <iostream>

#include "Lattice.h"

template <int N>
inline int deltaE(const Lattice<N>& lat, int i, int j) {
    return 2*lat.at(i, j)*lat.sum_neighbors(i, j);
}

inline double rand_double() {
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

template <int N>
void metropolis_step(
        Lattice<N>& lat,
        double beta,
        int& dE,
        int& dM) {
    int i = lat.random_site();
    int j = lat.random_site();

    dE = deltaE(lat, i, j);

    if (dE < 0 || (rand_double() < std::exp(-beta*dE))) {
        dM = -2*lat.flip(i, j);
    } else {
        dE = 0;
        dM = 0;
    }
}

template <int N>
void evolve(
        Lattice<N>& lat,
        int n,
        double beta) {
    int dE, dM;

    for (auto i = 0; i < n; i++) {
        metropolis_step(lat, beta, dE, dM);
    }
}

template <int N>
int evaluate_energy(const Lattice<N>& lat) {
    int total = 0;

    for (auto i = 0; i < N; i++) {
        for (auto j = 0; j < N; j++) {
            total -= lat.at(i,j)*(
                    lat.at(lat.nnplus(i), j) +
                    lat.at(i, lat.nnminus(j)));
        }
    }

    return total;
}

template <int N>
void time_average(
        Lattice<N>& lat,
        int n,
        double beta,
        double& en,
        double& mag) {
    int U = evaluate_energy(lat);
    int M = lat.magnetization();

    long U_tot = 0;
    long M_tot = 0;
    int dE, dM;

    for (auto i = 0; i < n; i++) {
        metropolis_step(lat, beta, dE, dM);
        U += dE;
        M += dM;
        U_tot += U;
        M_tot += M;
    }

    en = static_cast<double>(U_tot)/static_cast<double>(n);
    mag = static_cast<double>(M_tot)/static_cast<double>(n);
}

template <int N>
void ensemble_average(
        Lattice<N>& lat,
        double beta,
        int n_evolve,
        int n_average,
        double& en,
        double& mag) {
    evolve(lat, n_evolve, beta);
    time_average(lat, n_average, beta, en, mag);
}


int main(int, char**) {

    const int N = 64;
    Lattice<N> lat;
    std::random_device rd;
    srand(rd());

    double en, mag;
    std::vector<double> ts;
    std::vector<std::tuple<double, double, double>> data;

    double dt = (5 - 0.1)/400;
    double t = 0.1;
    for (int i = 0; i < 400; i ++) {
        ts.push_back(t);
        t += dt;
    }

    for (auto t: ts) {
        ensemble_average(lat, 1./t, 1000*N*N, 100*N*N, en, mag);

        if (fmod(t, 0.1) < 0.01) {
            std::cout << "T: " << t << std::endl;
        }
        data.push_back(std::make_tuple(t, en, mag));
    }

    std::ofstream data_file("met.dat");

    if (!data_file.is_open()) {
        std::cerr << "Could not open data file" << std::endl;
        return -1;
    } else {
        data_file << "T,M,U" << std::endl;
        for(auto row : data) {
            data_file << std::get<0>(row) << ", " <<
                std::get<1>(row) << ", " << std::get<2>(row) << std::endl;
        }
    }

    return 0;
}
