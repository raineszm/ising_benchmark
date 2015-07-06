#include <cstddef>
#include <cmath>
#include <array>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <queue>
#include <mutex>
#include <future>

#include <fstream>
#include <iostream>

#include <random>
#include <pthread.h>

#include "ising.h"
#include "Lattice.h"

const int N = 64;
const int NUM_THREADS = 4;

template <int N>
inline int deltaE(const Lattice<N>& lat, int i, int j) {
    return 2*lat.at(i, j)*lat.sum_neighbors(i, j);
}

inline double rand_double() {
    std::minstd_rand *rng =
        static_cast<std::minstd_rand*>(pthread_getspecific(rng_key));
    return std::generate_canonical<double, 8>(*rng);
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

using data_vector = std::vector<std::tuple<double, double, double>>;

data_vector metropolis_subset(std::queue<double>& ts, std::mutex& mtx) {
        Lattice<N> lat;
        std::random_device rd;
        std::minstd_rand rng(rd());
        pthread_setspecific(rng_key, &rng);

        double en, mag, t;
        data_vector local_data;

        while(true) {
            {
                std::unique_lock<std::mutex> lck(mtx);

                if (ts.empty()) {

                    break;
                }

                t = ts.front();
                ts.pop();

                if (fmod(t, 0.1) < 0.01) {
                    std::cout << "T: " << t << std::endl;
                }

            }

            ensemble_average(lat, 1./t, 1000*N*N, 100*N*N, en, mag);



            local_data.push_back(std::make_tuple(t, en, mag));
        }

        return local_data;
}

int main(int, char**) {

    Lattice<N> lat;
    std::random_device rd;
    pthread_key_create(&rng_key, NULL);

    std::queue<double> ts;
    data_vector data;

    //Build vector of temperatures
    double dt = (5 - 0.1)/400;
    double t = 0.1;
    for (int i = 0; i < 400; i ++) {
        ts.push(t);
        t += dt;
    }

    std::vector<std::future<data_vector>> futures;
    std::mutex mtx;

    for (auto i = 0; i < NUM_THREADS; i++) {
        futures.push_back(std::async(
            std::launch::async,
            [&ts, &mtx]() {
                return metropolis_subset(ts, mtx);
            }));
    }

    for (auto i = 0; i < NUM_THREADS; i++) {
        data_vector local_data = futures[i].get();
        std::move(local_data.begin(), local_data.end(), std::back_inserter(data));
    }

    pthread_key_delete(rng_key);

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
