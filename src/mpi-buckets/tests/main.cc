#include "CollisionManager.hpp"
#include "kokkos.hpp"

#include <mpi.h>
#include <random>
#include <vector>

template <class RNG>
std::vector<double> fillRandomPoints(double centerX,
                                     double centerY,
                                     double centerZ,
                                     double scale,
                                     size_t points,
                                     RNG&&  rng)
{
    std::normal_distribution<> dist(0, scale);
    std::vector<double>        pts;
    for (size_t i = 0; i < points; i++)
    {
        pts.push_back(dist(rng) + centerX);
        pts.push_back(dist(rng) + centerY);
        pts.push_back(dist(rng) + centerZ);
    }
    return pts;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);


    KokkosMock coord_d_;
    int        rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    coord_d_.values
        = fillRandomPoints(rank, 0.0, 0.0, 0.3, 100, std::mt19937_64{});



    double           background_grid_cell_size = 1.0;
    std::vector<int> exchange_members          = getExchangeMembers(
        coord_d_, background_grid_cell_size, MPI_COMM_WORLD, 0, 1, 2);

    for (int i = 0; i < exchange_members.size(); i++)
    {
        std::cout << i << std::endl;
    }

    MPI_Finalize();
}
