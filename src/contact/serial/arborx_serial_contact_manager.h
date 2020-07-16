#ifndef NIMBLE_ARBORX_SERIAL_CONTACT_MANAGER_H
#define NIMBLE_ARBORX_SERIAL_CONTACT_MANAGER_H

#include "serial_contact_manager.h"
#include <memory>

#include <ArborX.hpp>

namespace nimble {
/*  struct NarrowphaseResult
  {
    std::size_t first_global_id;
    std::size_t second_global_id;
  };
*/

  class ArborXSerialContactManager: public SerialContactManager
  {
  public:
    ArborXSerialContactManager(std::shared_ptr<ContactInterface> interface);
    ArborXSerialContactManager(const ArborXSerialContactManager &) = delete;
    ArborXSerialContactManager(ArborXSerialContactManager &&) noexcept;

    ArborXSerialContactManager &operator=( const ArborXSerialContactManager &) = delete;
    ArborXSerialContactManager &operator=( ArborXSerialContactManager &&) noexcept;

//    void InitializeContactVisualization(std::string const & contact_visualization_exodus_file_name) override {}
//    void ContactVisualizationWriteStep( double time ) override;

    ~ArborXSerialContactManager();

//    void ComputeContactForce(int step, bool debug_output) override;

  private:
      using bvh_type = ArborX::BVH< nimble_kokkos::kokkos_device_memory_space >;

/*    std::vector< NarrowphaseResult > m_last_results;
*/
  };
}

#endif // NIMBLE_ARBORX_SERIAL_CONTACT_MANAGER_H
