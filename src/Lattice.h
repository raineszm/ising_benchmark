#include <cstdlib>
#include <array>
#include <algorithm>

template <int N>
class Lattice {

    public:

        Lattice() {
            sites.fill(1);
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
            return (at(i, (j+1)%N)
                    + at(i, (N + j - 1)%N)
                    + at((i + 1)%N, j)
                    + at((N + i -1) % N, j));
        }

        int flip(int i, int j) {
            int s = at(i, j);
            sites[j + N*i] *= -s;
            return s;
        }

        int get_N() const { return N; }

        int random_site() {
            return rand() % N;
        }

        int magnetization() {
            return std::accumulate(sites.begin(), sites.end(), 0);
        }

    protected:
        std::array<int, N*N> sites;

};
