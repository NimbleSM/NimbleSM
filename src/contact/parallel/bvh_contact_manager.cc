#include "bvh_contact_manager.h"
#include <vt/transport.h>

namespace nimble {
  struct NarrowphaseResult
  {
    std::size_t idx;
  };

  struct NarrowphaseFunc {
    bvh::narrowphase_result operator()(const bvh::broadphase_collision< ContactEntity > &_a,
        const bvh::broadphase_collision< ContactEntity > &_b )
    {
      auto res = bvh::typed_narrowphase_result< NarrowphaseResult >();
      auto tree = build_snapshot_tree_top_down( _b.elements );

      for ( auto &&el : _a.elements )
        query_tree( tree, el, [&_a, &_b, &res]( std::size_t _ai, std::size_t _bi ) {
          res.emplace_back().idx = _b.elements[_bi].contact_entity_global_id();
        } );

      return res;
    }
  };

  BvhContactManager::BvhContactManager(std::shared_ptr<ContactInterface> interface, std::size_t _overdecomposition)
    : ParallelContactManager(interface),
      m_world{_overdecomposition},
      m_nodes{&m_world.create_collision_object()},
      m_faces{&m_world.create_collision_object()}
  {
    m_world.set_narrowphase_functor< ContactEntity >( NarrowphaseFunc{} );
  }

  BvhContactManager::BvhContactManager(BvhContactManager &&) noexcept = default;

  BvhContactManager &
  BvhContactManager::operator=( BvhContactManager &&) noexcept = default;

  BvhContactManager::~BvhContactManager() = default;

  void BvhContactManager::ComputeParallelContactForce(int step, bool debug_output) {
    m_world.start_iteration();

    // Update collision objects, this will build the trees
    m_nodes->set_entity_data(contact_nodes_);
    m_faces->set_entity_data(contact_faces_);

    m_nodes->broadphase(*m_faces);

    m_nodes->for_each_result< NarrowphaseResult >( []( const NarrowphaseResult &_res ){
      std::cout << "collision with " << _res.idx << '\n';
    } );

    m_world.finish_iteration();
  }

}