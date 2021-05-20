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

#include "nimble_kokkos_block_material_interface.h"

#include "nimble_kokkos_model_data.h"
#include "nimble_material.h"

namespace nimble_kokkos {

namespace {

inline void
compute_block_stress(
    const nimble::BlockData&                        block_data,
    const nimble_kokkos::DeviceFullTensorIntPtView& deformation_gradient_step_n_d,
    const nimble_kokkos::DeviceFullTensorIntPtView& deformation_gradient_step_np1_d,
    const nimble_kokkos::DeviceSymTensorIntPtView&  stress_step_n_d,
    nimble_kokkos::DeviceSymTensorIntPtView         stress_step_np1_d,
    const double                                    time_n,
    const double                                    time_np1)
{
  auto      material_device = block_data.material_device;
  auto      mdpolicy_2d     = make_elem_point_range_policy(block_data.num_elems, block_data.num_points_per_elem);
  const int def_grad_len    = 9;
  const int stress_len      = 6;
  Kokkos::parallel_for(
      "Stress", mdpolicy_2d, KOKKOS_LAMBDA(const int i_elem, const int i_ipt) {
        auto element_deformation_gradient_step_n_d =
            Kokkos::subview(deformation_gradient_step_n_d, i_elem, i_ipt, Kokkos::ALL());
        auto element_deformation_gradient_step_np1_d =
            Kokkos::subview(deformation_gradient_step_np1_d, i_elem, i_ipt, Kokkos::ALL());
        auto element_stress_step_n_d   = Kokkos::subview(stress_step_n_d, i_elem, i_ipt, Kokkos::ALL());
        auto element_stress_step_np1_d = Kokkos::subview(stress_step_np1_d, i_elem, i_ipt, Kokkos::ALL());
        material_device->GetStress(
            time_n,
            time_np1,
            element_deformation_gradient_step_n_d,
            element_deformation_gradient_step_np1_d,
            element_stress_step_n_d,
            element_stress_step_np1_d);
      });
}

}  // namespace

BlockMaterialInterface::BlockMaterialInterface(
    double                                time_n_,
    double                                time_np1_,
    const nimble::FieldIds&               field_ids_,
    const std::vector<nimble::BlockData>& blocks_,
    nimble::ModelDataBase*                model_data_)
    : nimble::BlockMaterialInterfaceBase(),
      time_n(time_n_),
      time_np1(time_np1_),
      field_ids(field_ids_),
      model_data(dynamic_cast<nimble_kokkos::ModelData*>(model_data_)),
      blocks(blocks_)
{
}

void
BlockMaterialInterface::ComputeStress() const
{
  for (auto&& block_data : blocks) {
    auto deformation_gradient_step_n_d = model_data->GetDeviceFullTensorIntegrationPointData(
        block_data.id, field_ids.deformation_gradient, nimble::STEP_N);
    auto deformation_gradient_step_np1_d = model_data->GetDeviceFullTensorIntegrationPointData(
        block_data.id, field_ids.deformation_gradient, nimble::STEP_NP1);
    auto stress_step_n_d =
        model_data->GetDeviceSymTensorIntegrationPointData(block_data.id, field_ids.stress, nimble::STEP_N);
    auto stress_step_np1_d =
        model_data->GetDeviceSymTensorIntegrationPointData(block_data.id, field_ids.stress, nimble::STEP_NP1);
    compute_block_stress(
        block_data,
        deformation_gradient_step_n_d,
        deformation_gradient_step_np1_d,
        stress_step_n_d,
        stress_step_np1_d,
        time_n,
        time_np1);
  }
}

}  // namespace nimble_kokkos
