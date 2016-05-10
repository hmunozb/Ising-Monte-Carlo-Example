#include "topologies.h"
using namespace std;
using namespace boost;
ising_graph<int> make_2D_ising(int l, int J){
	ising_graph<int> g(l*l);
	function < size_t(int, int)> nd = [l](int i, int j) {
		return j * l + i;
	};

	for (int i = 0; i < l; ++i) {
		for (int j = 0; j < l; ++j) {
			if (i > 0)
				add_edge(nd(i - 1, j), nd(i, j), J, g);
			if (i < l - 1)
				add_edge(nd(i, j), nd(i + 1, j), J, g);
			if (j > 0)
				add_edge(nd(i, j - 1), nd(i, j), J, g);
			if (j < l - 1)
				add_edge(nd(i, j), nd(i, j + 1), J, g);
		}
	}

	return g;
}

ising_graph<int> make_2D_ising_periodic(int l, int J) {
    ising_graph<int> g(l * l);
    function < size_t(int, int)> nd = [l](int i, int j) {
        return j * l + i;
    };
    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < l; ++j) {
            if (i > 0)
                add_edge(nd(i - 1, j), nd(i, j), J, g);
            add_edge(nd(i, j), nd((i + 1) % l, j), J, g);
            if (j > 0)
                add_edge(nd(i, j - 1), nd(i, j), J, g);
            add_edge(nd(i, j), nd(i, (j + 1) % l), J, g);
        }
    }

    return g;
}

ising_graph<int> make_2D_ising_spin_glass(int l, int m, unsigned int sd){
    randutils::random_generator<std::mt19937_64> rng;
    std::normal_distribution<double> normd(m, sd);
    auto rand_norm = [&](){ return normd(rng.engine());};
    
    ising_graph<int> g(l * l);
    function < size_t(int, int)> nd = [l](int i, int j) {
        return j * l + i;
    };
    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < l; ++j) {
            if (i > 0)
                add_edge(nd(i - 1, j), nd(i, j), (int)std::round(rand_norm()), g);
            add_edge(nd(i, j), nd((i + 1) % l, j), (int)std::round(rand_norm()), g);
            if (j > 0)
                add_edge(nd(i, j - 1), nd(i, j), (int)std::round(rand_norm()), g);
            add_edge(nd(i, j), nd(i, (j + 1) % l), (int)std::round(rand_norm()), g);
        }
    }

    return g;
    
}
