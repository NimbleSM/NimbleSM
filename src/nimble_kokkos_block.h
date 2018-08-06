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

#include "nimble_element.h"
#include "nimble_material.h"
#include "nimble_kokkos_data_manager.h"

#ifdef NIMBLE_HAVE_EXTRAS
  #include "nimble_ngp_lame_material.h"
#endif

#ifdef NIMBLE_HAVE_DARMA
  #include "darma.h"
#else
  #include <string>
  #include <vector>
  #include <memory>
#endif

namespace nimble_kokkos {

  class Block {

  public:

  Block() : macro_material_parameters_("none"), vol_ave_volume_offset_(-1), rve_boundary_condition_strategy_("none"),
      elem_conn_d("element_connectivity_d"), element_device_(0), material_device_(0) {}

    virtual ~Block() {
      /* Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) { */
      /*     element_device_->~Element(); */
      /*   }); */
      if (element_device_ != 0) {
        Kokkos::kokkos_free(element_device_);
      }
      if (material_device_ != 0) {
        Kokkos::kokkos_free(material_device_);
      }
    }

    void Initialize(std::string const & macro_material_parameters);

    /* void Initialize(std::string const & macro_material_parameters, */
    /*                 std::map<int, std::string> const & rve_material_parameters, */
    /*                 GenesisMesh const & rve_mesh, */
    /*                 std::string rve_boundary_condition_strategy); */

    void InstantiateMaterialModel();

    void InstantiateElement();

    void GetIntegrationPointDataLabelsAndLengths(std::vector< std::pair<std::string, nimble::Length> >& data_labels_and_lengths);

    double GetDensity() const {
      return material_host_->GetDensity();
    }

    double GetBulkModulus() const {
      return material_host_->GetBulkModulus();
    }

    std::shared_ptr<nimble::Element> GetHostElement() { return element_host_; }

    nimble::Element* GetDeviceElement() { return element_device_; }

    std::shared_ptr<nimble::Material> GetHostMaterialModel() { return material_host_; }

    nimble::Material* GetDeviceMaterialModel() { return material_device_; }

    DeviceElementConnectivityView& GetDeviceElementConnectivityView() { return elem_conn_d; }

    double ComputeCriticalTimeStep(const double * const node_coordinates,
                                   int num_elem,
                                   const int * const elem_conn) const;

    void InitializeElementData(int num_elem_in_block,
                               std::vector<int> const & elem_global_ids_in_block,
                               std::vector<int> const & rve_output_global_elem_ids,
                               std::vector<std::string> const & elem_data_labels,
                               std::vector<std::string> const & derived_elem_data_labels,
                               std::vector<double>& elem_data_n,
                               std::vector<double>& elem_data_np1,
                               nimble_kokkos::DataManager& data_manager);

    template <typename MatT>
    void ComputeTangentStiffnessMatrix(int num_global_unknowns,
                                       const double * const reference_coordinates,
                                       const double * const displacement,
                                       int num_elem,
                                       const int * const elem_conn,
                                       const int * const global_node_ids,
                                       MatT & tangent_stiffness) const ;

    void ComputeDerivedElementData(const double * const reference_coordinates,
                                   const double * const displacement,
                                   int num_elem,
                                   const int * const elem_conn,
                                   int num_elem_data,
                                   std::vector<double> const & elem_data_np1,
                                   int num_derived_elem_data,
                                   std::vector< std::vector<double> >& derived_elem_data);

  private:

    void DetermineDataOffsets(std::vector<std::string> const & elem_data_labels,
                              std::vector<std::string> const & derived_elem_data_labels);

    std::string macro_material_parameters_;
    std::map<int, std::string> rve_material_parameters_;
    std::string rve_boundary_condition_strategy_;
    std::vector<int> rve_output_global_elem_ids_;
    // todo: can we avoid carrying the rve_mesh around?
    //GenesisMesh rve_mesh_;
    std::vector<int> def_grad_offset_;
    std::vector<int> stress_offset_;
    std::vector<int> state_data_offset_;
    int vol_ave_volume_offset_;
    std::vector<int> vol_ave_offsets_;
    std::map<int, int> vol_ave_index_to_derived_data_index_;

    // element connectivity
    DeviceElementConnectivityView elem_conn_d;

    std::shared_ptr<nimble::Element> element_host_ = nullptr;
    nimble::Element* element_device_;


    std::shared_ptr<nimble::Material> material_host_ = nullptr;
    nimble::Material* material_device_;
  };

} // namespace nimble

#endif // NIMBLE_BLOCK_H
