#include <vector>
#include <random>

class Lattice {

    public:

        const int N;

        Lattice(int N_);

        int at(int i, int j) const;
        void fill(int s);

        const std::vector<int>& get_sites() const;

        const std::vector<int> neighbors(int i, int j) const;

        int flip(int i, int j);

        std::uniform_int_distribution<int> random_site;

    protected:
        std::vector<int> sites;

};
