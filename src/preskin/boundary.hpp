#include <Ioss_FaceGenerator.h>
#include <Ioss_Region.h>

#include <Kokkos_Core.hpp>
#include <string>
#include <unordered_map>
#include <vector>

namespace ps {
std::size_t
compute_boundary_nodes(
    const Ioss::Region&                                             region,
    const std::unordered_map<std::string, std::vector<Ioss::Face>>& boundary_map,
    Kokkos::View<double*>&                                          coords,
    Kokkos::View<int*>&                                             ids);
}