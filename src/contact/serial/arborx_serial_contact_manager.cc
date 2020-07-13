#include "arborx_serial_contact_manager.h"
#include <fstream>

namespace nimble {

/*  struct NarrowphaseFunc {
    bvh::narrowphase_result operator()(const bvh::broadphase_collision< ContactEntity > &_a,
        const bvh::broadphase_collision< ContactEntity > &_b )
    {
      auto res = bvh::typed_narrowphase_result< NarrowphaseResult >();
      auto tree = build_snapshot_tree_top_down( _b.elements );

      for ( auto &&el : _a.elements )
        query_tree( tree, el, [&_a, &_b, &res]( std::size_t _ai, std::size_t _bi ) {
          auto &entry = res.emplace_back();
          entry.first_global_id = _ai;
          entry.second_global_id = _bi;
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
    for ( auto &&node : contact_nodes_ )
      node.RecomputeKdop();
    for ( auto &&face : contact_faces_ )
      face.RecomputeKdop();

    m_nodes->set_entity_data(contact_nodes_);
    m_faces->set_entity_data(contact_faces_);

    m_nodes->broadphase(*m_faces);

    m_last_results.clear();
    m_nodes->for_each_result< NarrowphaseResult >( [this]( const NarrowphaseResult &_res ){
      m_last_results.emplace_back( _res );
    } );

    m_world.finish_iteration();
  }
*/
}