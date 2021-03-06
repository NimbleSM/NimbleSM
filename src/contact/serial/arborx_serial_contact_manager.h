/*
//@HEADER
// ************************************************************************
//
//                                NimbleSM
//                             Copyright 2018
//   National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
// retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
// NO EVENT SHALL NTESS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact David Littlewood (djlittl@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef NIMBLE_ARBORX_SERIAL_CONTACT_MANAGER_H
#define NIMBLE_ARBORX_SERIAL_CONTACT_MANAGER_H

#ifdef NIMBLE_HAVE_ARBORX

#include <memory>

#include "serial_contact_manager.h"

namespace nimble_kokkos {

class ModelData;

}

namespace nimble {

class DataManager;

class ArborXSerialContactManager : public SerialContactManager
{
 public:
  ArborXSerialContactManager(std::shared_ptr<ContactInterface> interface, nimble::DataManager& data_manager);
  ArborXSerialContactManager(const ArborXSerialContactManager&) = delete;

  ~ArborXSerialContactManager() override = default;

  void
  ComputeSerialContactForce(int step, bool debug_output, nimble::Viewify<2> contact_force) override;

 private:
  template <typename ContactManagerType>
  friend struct ArborXCallback;

  void
  updateCollisionData(
      Kokkos::View<int*, nimble_kokkos::kokkos_device>& indices,
      Kokkos::View<int*, nimble_kokkos::kokkos_device>& offset);

  nimble_kokkos::ModelData* model_data = nullptr;
};
}  // namespace nimble

#endif  // NIMBLE_HAVE_ARBORX

#endif  // NIMBLE_ARBORX_SERIAL_CONTACT_MANAGER_H
