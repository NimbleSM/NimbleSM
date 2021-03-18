#ifndef NIMBLE_BVH_CONTACT_MANAGER_H
#define NIMBLE_BVH_CONTACT_MANAGER_H

#include "parallel_contact_manager.h"
#include <bvh/collision_object.hpp>
#include <bvh/collision_world.hpp>
#include <memory>

namespace nimble {
struct NarrowphaseResult {
  std::size_t local_index;
  double gap;
  double contact_force[3];
  double bary[3];
  bool node;
};

class BvhContactManager : public ParallelContactManager {
public:
  BvhContactManager(std::shared_ptr<ContactInterface> interface,
                    nimble::DataManager &data_manager,
                    std::size_t _overdecomposition);
  BvhContactManager(const BvhContactManager &) = delete;
  BvhContactManager(BvhContactManager &&);

  BvhContactManager &operator=(const BvhContactManager &) = delete;
  BvhContactManager &operator=(BvhContactManager &&);

  ~BvhContactManager();

  void ComputeParallelContactForce(int step, bool debug_output,
                                   nimble::Viewify<2> contact_force) override;

private:
  bvh::collision_world m_world;
  bvh::collision_object *m_nodes;
  bvh::collision_object *m_faces;

  std::vector<NarrowphaseResult> m_last_results;
};
} // namespace nimble

#endif // NIMBLESM_BVH_CONTACT_MANAGER_H
