#include "boundary.hpp"

#include <Ioss_NodeBlock.h>

#include <Kokkos_Core_fwd.hpp>

namespace ps {
std::size_t
compute_boundary_nodes(
    const Ioss::Region&                                             region,
    const std::unordered_map<std::string, std::vector<Ioss::Face>>& boundary_map,
    Kokkos::View<double*>&                                          coords,
    Kokkos::View<int*>&                                             ids)
{
  auto&             db        = *region.get_database();
  const std::size_t num_nodes = static_cast<std::size_t>(region.get_property("node_count").get_int());

  // Initialized to zero...
  Kokkos::View<std::size_t*> node_referenced("node_referenced", num_nodes);
  auto                       node_referenced_h = Kokkos::create_mirror_view(node_referenced);

  for (auto&& [_, boundary] : boundary_map) {
    for (const auto& face : boundary) {
      for (auto global_node_id : face.connectivity_) {
        if (global_node_id != 0) {
          const auto local_id         = db.node_global_to_local(global_node_id, true) - 1;
          node_referenced_h(local_id) = 1;
        }
      }
    }
  }

  Kokkos::deep_copy(node_referenced, node_referenced_h);
  Kokkos::View<std::size_t*> node_ref_idx("node_ref_idx", num_nodes);

  // Stream compression
  std::size_t           count;
  Kokkos::RangePolicy<> rp(0, num_nodes);
  Kokkos::parallel_scan(
      "compression_reindex",
      rp,
      KOKKOS_LAMBDA(std::int64_t i, std::size_t & index, bool final) {
        auto val = node_referenced(i);

        if (final) node_ref_idx(i) = index;

        index += val;
      },
      count);

  // Get the total number
  count += node_referenced_h(num_nodes - 1);

  std::cout << "got " << count << " total nodes referenced\n";

  auto* node_block = region.get_node_blocks()[0];
  if (!node_block) throw std::runtime_error("mesh must have a node block");

  Kokkos::View<double*> coords_in("coords_in", 3 * num_nodes);
  auto                  coords_in_h = Kokkos::create_mirror_view(coords_in);

  auto num_coords_read =
      node_block->get_field_data("mesh_model_coordinates", coords_in_h.data(), coords_in_h.extent(0) * sizeof(double));
  std::cout << "got " << num_coords_read << " out of " << num_nodes << " coords\n";
  if (num_nodes != num_coords_read) throw std::runtime_error("invalid number of nodes");
  Kokkos::deep_copy(coords_in, coords_in_h);

  Kokkos::View<int*> ids_in("ids_in", num_nodes);
  auto               ids_in_h = Kokkos::create_mirror_view(ids_in);

  auto num_ids_read =
      node_block->get_field_data("ids", ids_in_h.data(), ids_in_h.extent(0) * sizeof(int));
  std::cout << "got " << num_ids_read << " out of " << num_nodes << " coords\n";
  if (num_nodes != num_ids_read) throw std::runtime_error("invalid number of ids");
  Kokkos::deep_copy(ids_in, ids_in_h);

  Kokkos::resize(coords, 3 * count);
  Kokkos::resize(ids, count);
  Kokkos::parallel_for(
      rp, KOKKOS_LAMBDA(std::int64_t i) {
        if (node_referenced(i) > 0) {
          const auto j      = node_ref_idx(i);
          coords(3 * j + 0) = coords_in(3 * i + 0);
          coords(3 * j + 1) = coords_in(3 * i + 1);
          coords(3 * j + 2) = coords_in(3 * i + 2);
          ids(j)            = ids_in(i);
        }
      });

  Kokkos::fence();

  std::cout << "will write " << count << " nodes\n";
  return count;
}
}  // namespace ps
