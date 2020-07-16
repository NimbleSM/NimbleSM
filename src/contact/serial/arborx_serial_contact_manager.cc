#include "arborx_serial_contact_manager.h"
#include <fstream>

namespace nimble {

  ArborXSerialContactManager::ArborXSerialContactManager(std::shared_ptr<ContactInterface> interface):
  SerialContactManager(interface)
  {
      // Here will be initialisation stuffs
      ArborX::Point point;
  }

  void ArborXSerialContactManager::ComputeSerialContactForce(int step, bool debug_output) {
  }

//  ArborXSerialContactManager::ArborXSerialContactManager(ArborXSerialContactManager &&) noexcept = default;
//  ArborXSerialContactManager & ArborXSerialContactManager::operator=( ArborXSerialContactManager &&) noexcept = default;
//  ArborXSerialContactManager::~ArborXSerialContactManager() = default;

}