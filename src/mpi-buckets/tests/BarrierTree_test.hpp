#include <mpi/BarrierTree.hpp>
#include <algorithm>

using std::vector;

struct Server {
     
};

template<class Gen>
vector<int> getDestRanks(int n, double p, Gen&& g) {
    auto dist = std::uniform_real_distribution<>(); 
    vector<int> ranks(n); 
    std::iota(ranks.begin(), ranks.end(), 0); 
    auto new_end = std::remove_if(ranks.begin(), ranks.end(), [&](int) { return dist(g) < p; }); 
    ranks.erase(new_end, ranks.end()); 
    return ranks; 
}

void runTest(MPI_Comm comm, int tag) {
    BarrierTree tree(comm, tag); 
    
    
}


int main() {
    
    

}
