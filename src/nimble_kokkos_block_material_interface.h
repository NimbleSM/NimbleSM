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

#ifndef SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_H_
#define SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_H_

#include <nimble_data_utils.h>
#include <nimble_kokkos_block.h>
#include <nimble_kokkos_data_manager.h>
#include <nimble_kokkos_defs.h>
#include <nimble_utils.h>

#ifdef NIMBLE_HAVE_EXTRAS
#include <nimble_ngp_lame_material.h>
#endif

namespace nimble_kokkos
{

struct BlockData {
  BlockData(nimble_kokkos::Block &block_, nimble::Material &material_d_, const int block_id_,
            const int num_block_elems_, const int num_points_per_block_elem_)
      :
      block(block_),
      material_device(material_d_),
      id(block_id_),
      num_elems(num_block_elems_),
      num_points_per_elem(num_points_per_block_elem_) {
  }

  nimble_kokkos::Block &block;
  nimble::Material &material_device;
  int id;
  int num_elems;
  int num_points_per_elem;
};

class BlockMaterialInterface {
 public:
  BlockMaterialInterface(const double time_n_, const double time_np1_, const FieldIds &field_ids_,
                         const std::vector<BlockData>& blocks_,
                         nimble_kokkos::ModelData &model_data_)
      :
      time_n(time_n_),
      time_np1(time_np1_),
      field_ids(field_ids_),
      model_data(model_data_),
      blocks(blocks_) {
  }

  virtual ~BlockMaterialInterface() = default;

  using ElemPointRangePolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2> >;
  static inline ElemPointRangePolicy make_elem_point_range_policy(const int num_block_elems, const int num_points_per_elem)
  {
    return ElemPointRangePolicy( { 0, 0 }, { num_block_elems, num_points_per_elem });
  }

  virtual void ComputeStress() const {
    for (auto &&block_data : blocks) {
      auto deformation_gradient_step_n_d = model_data.GetDeviceFullTensorIntegrationPointData(
          block_data.id, field_ids.deformation_gradient, nimble::STEP_N);
      auto deformation_gradient_step_np1_d = model_data.GetDeviceFullTensorIntegrationPointData(
          block_data.id, field_ids.deformation_gradient, nimble::STEP_NP1);
      auto stress_step_n_d = model_data.GetDeviceSymTensorIntegrationPointData(block_data.id, field_ids.stress,
                                                                               nimble::STEP_N);
      auto stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_data.id, field_ids.stress,
                                                                                 nimble::STEP_NP1);
      auto mdpolicy_2d = make_elem_point_range_policy(block_data.num_elems, block_data.num_points_per_elem);
      Kokkos::parallel_for(
          "Stress",
          mdpolicy_2d,
          KOKKOS_LAMBDA (const int i_elem, const int i_ipt) {
            auto element_deformation_gradient_step_n_d = Kokkos::subview(deformation_gradient_step_n_d, i_elem, i_ipt,
                                                                         Kokkos::ALL());
            auto element_deformation_gradient_step_np1_d = Kokkos::subview(deformation_gradient_step_np1_d, i_elem,
                                                                           i_ipt, Kokkos::ALL());
            auto element_stress_step_n_d = Kokkos::subview(stress_step_n_d, i_elem, i_ipt, Kokkos::ALL());
            auto element_stress_step_np1_d = Kokkos::subview(stress_step_np1_d, i_elem, i_ipt, Kokkos::ALL());
            block_data.material_device.GetStress(time_n, time_np1, element_deformation_gradient_step_n_d,
                                                 element_deformation_gradient_step_np1_d, element_stress_step_n_d,
                                                 element_stress_step_np1_d);
          });
    }
  }

 protected:
  const double time_n;
  const double time_np1;
  const FieldIds &field_ids;
  nimble_kokkos::ModelData &model_data;
  std::vector<BlockData> blocks;
};

class ExtrasBlockMaterialInterface : public BlockMaterialInterface {
 public:
  ExtrasBlockMaterialInterface(const double time_n_, const double time_np1_, const FieldIds &field_ids_,
                               const std::vector<BlockData> &blocks_, nimble_kokkos::ModelData &model_data_)
      :
      BlockMaterialInterface(time_n_, time_np1_, field_ids_, blocks_, model_data_) {
#ifndef NIMBLE_HAVE_EXTRAS
    throw std::logic_error("NGP LAME interface not available!")
#else
    const double delta_time = time_np1 - time_n;
    for (auto &&block_data : blocks) {
      auto deformation_gradient_step_n_d = model_data.GetDeviceFullTensorIntegrationPointData(
          block_data.id, field_ids.deformation_gradient, nimble::STEP_N);
      auto deformation_gradient_step_np1_d = model_data.GetDeviceFullTensorIntegrationPointData(
          block_data.id, field_ids.deformation_gradient, nimble::STEP_NP1);
      auto unrotated_stress_step_n_d = model_data.GetDeviceSymTensorIntegrationPointData(block_data.id,
                                                                                         field_ids.unrotated_stress,
                                                                                         nimble::STEP_N);
      auto unrotated_stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_data.id,
                                                                                           field_ids.unrotated_stress,
                                                                                           nimble::STEP_NP1);
      auto ngp_lame_data = *(block_data.block.GetNGPLAMEData());
      auto mdpolicy_2d = make_elem_point_range_policy(block_data.num_elems, block_data.num_points_per_elem);
      Kokkos::parallel_for(
          "NGP LAME Kinematics",
          mdpolicy_2d,
          KOKKOS_LAMBDA (const int i_elem, const int i_ipt) {
            auto element_deformation_gradient_step_n_d = Kokkos::subview(deformation_gradient_step_n_d, i_elem, i_ipt,
                                                                         Kokkos::ALL());
            auto element_deformation_gradient_step_np1_d = Kokkos::subview(deformation_gradient_step_np1_d, i_elem,
                                                                           i_ipt, Kokkos::ALL());
            auto element_unrotated_stress_step_n_d = Kokkos::subview(unrotated_stress_step_n_d, i_elem, i_ipt,
                                                                     Kokkos::ALL());
            double disp_grad_step_n[9], disp_grad_step_np1[9];
            for (int i = 0; i < 9; i++) {
              disp_grad_step_n[i] = element_deformation_gradient_step_n_d(i);
              disp_grad_step_np1[i] = element_deformation_gradient_step_np1_d(i);
              if (i < 3) {
                disp_grad_step_n[i] -= 1.0;
                disp_grad_step_np1[i] -= 1.0;
              }
            }
            int i_mat_pt = i_elem * block_data.num_points_per_elem + i_ipt;
            for (int i = 0; i < 9; i++) {
              ngp_lame_data.disp_grad_(i_mat_pt, i) = disp_grad_step_np1[i];
              ngp_lame_data.velo_grad_(i_mat_pt, i) = (disp_grad_step_np1[i] - disp_grad_step_n[i]) / delta_time;
            }
            for (int i = 0; i < 6; i++) {
              ngp_lame_data.stress_old_(i_mat_pt, i) = element_unrotated_stress_step_n_d(i);
            }
            // TODO handle state variables
          });
    }
#endif
  }

#ifdef NIMBLE_HAVE_EXTRAS
  void ComputeStress() const override {
    // Create a view of NGP LAME material models for use in team loop
    Kokkos::View<mtk_ngp::DevicePtr<lame::ngp::Material>*> ngp_lame_materials_d("NGP LAME Material Models",
                                                                                blocks.size());
    auto ngp_lame_materials_h = Kokkos::create_mirror_view(ngp_lame_materials_d);
    int block_index = 0;
    for (auto &&block_data : blocks) {
      auto nimble_ngp_lame_material_model = std::dynamic_pointer_cast<nimble::NGPLAMEMaterial>(
          block_data.block.GetHostMaterialModel());
      ngp_lame_materials_h(block_index++) = nimble_ngp_lame_material_model->GetNGPLAMEMaterialModel();
    }
    Kokkos::deep_copy(ngp_lame_materials_d, ngp_lame_materials_h);

    // Compute stress
    const Kokkos::TeamPolicy<> team_loop(blocks.size(), Kokkos::AUTO);
    Kokkos::parallel_for("NGP LAME Stress", team_loop, KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type &team) {
      int team_id = team.league_rank();
      ngp_lame_materials_d[team_id]->compute_stress(team);
    });
  }

  ~ExtrasBlockMaterialInterface() {
    // Copy stress back to the Nimble data structures
    for (auto &&block_data : blocks) {
      auto ngp_lame_data = *(block_data.block.GetNGPLAMEData());

      nimble_kokkos::DeviceSymTensorIntPtView unrotated_stress_step_np1_d = model_data
          .GetDeviceSymTensorIntegrationPointData(block_data.id, field_ids.unrotated_stress, nimble::STEP_NP1);
      nimble_kokkos::DeviceSymTensorIntPtView stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(
          block_data.id, field_ids.stress, nimble::STEP_NP1);

      auto mdpolicy_2d = make_elem_point_range_policy(block_data.num_elems, block_data.num_points_per_elem);
      Kokkos::parallel_for(
          "NGP LAME Stress Copy",
          mdpolicy_2d,
          KOKKOS_LAMBDA (const int i_elem, const int i_ipt) {
            auto element_unrotated_stress_step_np1_d = Kokkos::subview(unrotated_stress_step_np1_d, i_elem, i_ipt,
                                                                       Kokkos::ALL());
            auto element_stress_step_np1_d = Kokkos::subview(stress_step_np1_d, i_elem, i_ipt, Kokkos::ALL());
            double unrotated_stress[6], def_grad[9], left_stretch[6], rotation[9], rotation_transpose[9],
                temp_full_tensor_1[9], temp_full_tensor_2[9];
            int i_mat_pt = i_elem * block_data.num_points_per_elem + i_ipt;
            for (int i = 0; i < 6; i++) {
              unrotated_stress[i] = ngp_lame_data.stress_new_(i_mat_pt, i);
            }
            // Convert LAME unrotated stress to stress
            // TODO Figure out a way to reuse rotation calculation in the NGP LAME material model, if available
            for (int i = 0; i < 9; i++) {
              def_grad[i] = ngp_lame_data.disp_grad_(i_mat_pt, i);
              if (i < 3) {
                def_grad[i] += 1.0;
              }
            }
            Polar_Decomp(def_grad, left_stretch, rotation);
            rotation_transpose[K_F_XX] = rotation[K_F_XX];
            rotation_transpose[K_F_XY] = rotation[K_F_YX];
            rotation_transpose[K_F_XZ] = rotation[K_F_ZX];
            rotation_transpose[K_F_YX] = rotation[K_F_XY];
            rotation_transpose[K_F_YY] = rotation[K_F_YY];
            rotation_transpose[K_F_YZ] = rotation[K_F_ZY];
            rotation_transpose[K_F_ZX] = rotation[K_F_XZ];
            rotation_transpose[K_F_ZY] = rotation[K_F_YZ];
            rotation_transpose[K_F_ZZ] = rotation[K_F_ZZ];
            Mult_Sym33_Full33(unrotated_stress, rotation_transpose, temp_full_tensor_1);
            Mult_Full33_Full33(rotation, temp_full_tensor_1, temp_full_tensor_2);
            for (int i = 0; i < 6; i++) {
              element_unrotated_stress_step_np1_d(i) = unrotated_stress[i];
              element_stress_step_np1_d(i) = temp_full_tensor_2[i];
            }
            // TODO handle state variables
          });
    }
  }
#endif

};

}

#endif /* SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_H_ */
