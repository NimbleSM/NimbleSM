#ifndef NIMBLE_BVH_CONTACT_MANAGER_H
#define NIMBLE_BVH_CONTACT_MANAGER_H

#include "parallel_contact_manager.h"
#include <bvh/collision_world.hpp>
#include <bvh/collision_object.hpp>
#include <memory>


namespace nimble {
  struct NarrowphaseResult
  {
    std::size_t first_global_id;
    std::size_t second_global_id;
    double gap;
  };

  class BvhContactManager : public ParallelContactManager
  {
  public:

    BvhContactManager(std::shared_ptr<ContactInterface> interface, std::size_t _overdecomposition);
    BvhContactManager(const BvhContactManager &) = delete;
    BvhContactManager(BvhContactManager &&) noexcept;

    BvhContactManager &operator=( const BvhContactManager &) = delete;
    BvhContactManager &operator=( BvhContactManager &&) noexcept;

    ~BvhContactManager();


    void ComputeParallelContactForce(int step, bool debug_output) override;

  private:

    bvh::collision_world m_world;
    bvh::collision_object *m_nodes;
    bvh::collision_object *m_faces;

    std::vector< NarrowphaseResult > m_last_results;
  };
}

#endif // NIMBLESM_BVH_CONTACT_MANAGER_H
