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

#ifndef NIMBLE_KOKKOS_BLOCK_H
#define NIMBLE_KOKKOS_BLOCK_H

#include <map>
#include <vector>

#include "nimble_block_base.h"

#include "nimble_element.h"
#include "nimble_kokkos_defs.h"

#ifdef NIMBLE_HAVE_DARMA
  #include "darma.h"
#else
  #include <string>
  #include <vector>
  #include <memory>
#endif

namespace nimble { struct NGPLAMEData; }

namespace nimble_kokkos { class MaterialFactory; }

namespace nimble_kokkos {

class Block : public nimble::BlockBase {

  public:

    Block() : BlockBase(),
      elem_conn_d("element_connectivity_d", 0),
      element_device_(0), material_device_(0)
    {}

    ~Block() override {
      /* Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) { */
      /*     element_device_->~Element(); */
      /*   }); */
      if (element_device_ != nullptr) {
        Kokkos::kokkos_free(element_device_);
      }
      if (material_device_ != nullptr) {
        Kokkos::kokkos_free(material_device_);
      }
    }

    void Initialize(std::string const & macro_material_parameters,
                    int num_elements,
                    MaterialFactory& factory);

    void InstantiateMaterialModel(int num_material_points,
                                  MaterialFactory& factory);

    void InstantiateElement() override;

    std::shared_ptr<nimble::Element> GetHostElement() { return element_; }

    nimble::Element* GetDeviceElement() { return element_device_; }

    std::shared_ptr<nimble::Material> GetHostMaterialModel() { return material_; }

    nimble::Material* GetDeviceMaterialModel() { return material_device_; }

    DeviceElementConnectivityView& GetDeviceElementConnectivityView() { return elem_conn_d; }

    double ComputeCriticalTimeStep(const double * const node_reference_coordinates,
                                   const double * const node_displacements,
                                   int num_elem,
                                   const int * const elem_conn) const;

    template <typename MatT>
    void ComputeTangentStiffnessMatrix(int num_global_unknowns,
                                       const double * const reference_coordinates,
                                       const double * const displacement,
                                       int num_elem,
                                       const int * const elem_conn,
                                       const int * const global_node_ids,
                                       MatT & tangent_stiffness) const ;

    std::shared_ptr<nimble::NGPLAMEData> GetNGPLAMEData() { return ngp_lame_data_; }

  private:

    // element connectivity
    DeviceElementConnectivityView elem_conn_d;

    nimble::Element* element_device_;

    nimble::Material* material_device_;

    std::shared_ptr<nimble::NGPLAMEData> ngp_lame_data_;

  };

} // namespace nimble

#endif // NIMBLE_BLOCK_H
