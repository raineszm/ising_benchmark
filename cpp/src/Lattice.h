#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <random>

template<int N>
class Lattice
{

public:
  Lattice()
    : rng(std::chrono::system_clock::now().time_since_epoch().count())
  {
    sites.fill(1);
    for (auto i = 0; i < N; i++) {
      nnplus_[i] = (i + 1) % N;
      nnminus_[i] = (i + N - 1) % N;
    }
  }

  int at(int i, int j) const { return sites[j + N * i]; }

  void fill(int s) { sites.fill(s); }

  const std::array<int, N * N>& get_sites() const { return sites; }

  int sum_neighbors(int i, int j) const
  {
    return (at(i, nnplus_[j]) + at(i, nnminus_[j]) + at(nnplus_[i], j) +
            at(nnminus_[i], j));
  }

  int flip(int i, int j)
  {
    int s = at(i, j);
    sites[j + N * i] *= -1;
    return s;
  }

  int get_N() const { return N; }

  inline int random_site() { return rng() / (rng.max() / N + 1); }

  int magnetization() const
  {
    return std::accumulate(sites.begin(), sites.end(), 0);
  }

  int energy() const
  {
    int total = 0;

    for (auto i = 0; i < N; i++) {
      for (auto j = 0; j < N; j++) {
        total -= at(i, j) * (at(nnplus(i), j) + at(i, nnminus(j)));
      }
    }

    return total;
  }

  int nnplus(int k) const { return nnplus_[k]; }

  int nnminus(int k) const { return nnminus_[k]; }

  inline int deltaE(int i, int j) const
  {
    return 2 * at(i, j) * sum_neighbors(i, j);
  }

  std::minstd_rand rng;

protected:
  std::array<int, N * N> sites;

  std::array<int, N> nnplus_;
  std::array<int, N> nnminus_;
};
