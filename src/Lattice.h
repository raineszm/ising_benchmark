#include <vector>

class Lattice {

    public:

        const int N;

        Lattice(int N_);

        int at(int i, int j) const;
        void fill(int s);

        const std::vector<int>& get_sites() const;


    protected:
        std::vector<int> sites;

};
