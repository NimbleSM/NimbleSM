#include "arborx_serial_contact_manager.h"
#include <fstream>

#include <ArborX.hpp>

namespace nimble {

  ArborXSerialContactManager::ArborXSerialContactManager(std::shared_ptr<ContactInterface> interface):
  SerialContactManager(interface)
  {
      // Here will be initialisation stuffs
      ArborX::Point test;
      //Kokkos::View<>


  }

  //ArborXSerialContactManager::ArborXSerialContactManager(ArborXSerialContactManager &&) noexcept = default;
  //ArborXSerialContactManager &
  //ArborXSerialContactManager::operator=( ArborXSerialContactManager &&) noexcept = default;
  ArborXSerialContactManager::~ArborXSerialContactManager() = default;

//  void ArborXSerialContactManager::ComputeContactForce(int step, bool debug_output) {
//  }
}