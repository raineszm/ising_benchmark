#pragma once
#include <cstdlib>
#include <array>
#include <algorithm>
#include <random>

#include "ising.h"

template <int N>
class Lattice {

    public:

        Lattice() {
            sites.fill(1);
            nnplus_[N-1] = 0;
            nnminus_[0] = N - 1;

            for(auto i = 0; i < N - 1; i++) {
                nnplus_[i] = i + 1;
                nnminus_[N - i] = N - i - 1;
            }
        }

        int at(int i, int j) const {
            return sites[j + N*i];
        }


        void fill(int s) {
            sites.fill(s);
        }


        const std::array<int, N*N>& get_sites() const {
            return sites;
        }


        int sum_neighbors(int i, int j) const {
            return (at(i, nnplus_[j])
                    + at(i, nnminus_[j])
                    + at(nnplus_[i], j)
                    + at(nnminus_[i], j));
        }

        int flip(int i, int j) {
            int s = at(i, j);
            sites[j + N*i] *= -1;
            return s;
        }

        int get_N() const { return N; }

        inline int random_site() const {
            return rng()/(rng.max()/N + 1);
        }

        int magnetization() const {
            return std::accumulate(sites.begin(), sites.end(), 0);
        }

        int nnplus(int k) const {
            return nnplus_[k];
        }

        int nnminus(int k) const {
            return nnminus_[k];
        }


    protected:
        std::array<int, N*N> sites;

        std::array<int, N> nnplus_;
        std::array<int, N> nnminus_;

};
