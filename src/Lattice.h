#include <cstdlib>
#include <array>
#include <algorithm>

template <int N>
class Lattice {

    public:

        Lattice() {
            sites.fill(1);
            nnplus[N-1] = 0;
            nnminus[0] = N - 1;

            for(auto i = 0; i < N - 1; i++) {
                nnplus[i] = i + 1;
                nnminus[N - i] = N - i - 1;
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
            return (at(i, nnplus[j])
                    + at(i, nnminus[j])
                    + at(nnplus[i], j)
                    + at(nnminus[i], j));
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

        std::array<int, N> nnplus;
        std::array<int, N> nnminus;

};
