#include <cstddef>
#include <cmath>
#include <tuple>
#include <queue>
#include <mutex>
#include <thread>

#include <fstream>
#include <iostream>

#include <random>

#include "Lattice.h"
#include "Channel.h"

const double T0 = 0.1;
const double TF = 5.;
const int N = 64;
const int STEPS = 400;
const int NUM_THREADS = static_cast<int>(std::thread::hardware_concurrency());

inline double rand_double(std::minstd_rand& rng) {
    return std::generate_canonical<double, 8>(rng);
}

template <int N>
auto metropolis_step(
        Lattice<N>& lat,
        double beta) {
    int i = lat.random_site();
    int j = lat.random_site();

    int dE = lat.deltaE(i, j);

    if (dE < 0 || (rand_double(lat.rng) < std::exp(-beta*dE))) {
        return std::make_tuple(dE,
                -2*lat.flip(i, j));
    } else {
        return std::make_tuple(0, 0);
    }

}

template <int N>
void evolve(
        Lattice<N>& lat,
        int n,
        double beta) {

    for (auto i = 0; i < n; i++) {
        metropolis_step(lat, beta);
    }
}

inline
double average(long agg, int n) {
    return static_cast<double>(agg) / n;
}

template <int N>
auto time_average(
        Lattice<N>& lat,
        int n,
        double beta) {
    int U = lat.energy();
    int M = lat.magnetization();

    long U_tot = 0;
    long M_tot = 0;
    int dE, dM;

    for (auto i = 0; i < n; i++) {
        std::tie(dE, dM) = metropolis_step(lat, beta);
        U += dE;
        M += dM;
        U_tot += U;
        M_tot += M;
    }

    return std::make_tuple(average(M_tot, n), average(U_tot, n));
}

template <int N>
auto ensemble_average(
        Lattice<N>& lat,
        double beta,
        int n_evolve,
        int n_average) {
    evolve(lat, n_evolve, beta);
    return time_average(lat, n_average, beta);
}

using data_type = std::tuple<double, double, double>;

void metropolis_subset(std::queue<double>& ts, std::mutex& mtx,
        Channel<data_type>& chan) {
        Lattice<N> lat;

        while(true) {
            double t;

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

            auto observables = ensemble_average(lat, 1./t, 1000*N*N, 100*N*N);

            chan.put(std::tuple_cat(std::make_tuple(t), observables));
        }
}

int main(int, char**) {

    Lattice<N> lat;

    std::queue<double> ts;
    Channel<data_type> chan;

    //Build vector of temperatures
    double dt = (TF - T0)/(STEPS - 1);
    double t = T0;
    for (auto i = 0; i < STEPS; i++) {
        ts.push(t);
        t += dt;
    }

    std::vector<std::thread> threads;
    std::mutex mtx;

    for (auto i = 0; i < NUM_THREADS; i++) {
        threads.push_back(std::thread(
            [&ts, &mtx, &chan]() {
                return metropolis_subset(ts, mtx, chan);
            }));
    }

    std::ofstream data_file("met.dat");

    if (!data_file.is_open()) {
        std::cerr << "Could not open data file" << std::endl;
        return -1;
    } else {
        data_file << "T,M,U" << std::endl;
        for (int i = 0; i < STEPS; i++) {
            data_type row = chan.take();
            data_file << std::get<0>(row) << ", " <<
                std::get<1>(row) << ", " << std::get<2>(row) << std::endl;
        }
    }

    for (std::thread& t: threads)
        t.join();

    return 0;
}
