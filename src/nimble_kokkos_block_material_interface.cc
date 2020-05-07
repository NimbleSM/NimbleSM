#include <nimble_kokkos_block_material_interface.h>
#include <nimble_material.h>

namespace nimble_kokkos {

namespace {

inline void compute_block_stress(const BlockData& block_data,
                                 nimble_kokkos::DeviceFullTensorIntPtView deformation_gradient_step_n_d,
                                 nimble_kokkos::DeviceFullTensorIntPtView deformation_gradient_step_np1_d,
                                 nimble_kokkos::DeviceSymTensorIntPtView stress_step_n_d,
                                 nimble_kokkos::DeviceSymTensorIntPtView stress_step_np1_d,
                                 const double time_n,
                                 const double time_np1) {
  auto material_device = block_data.material_device;
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
        material_device->GetStress(time_n, time_np1, element_deformation_gradient_step_n_d,
                                   element_deformation_gradient_step_np1_d, element_stress_step_n_d,
                                   element_stress_step_np1_d);
      });
}

}

void BlockMaterialInterface::ComputeStress() const {
  for (auto &&block_data : blocks) {
    auto deformation_gradient_step_n_d = model_data.GetDeviceFullTensorIntegrationPointData(
        block_data.id, field_ids.deformation_gradient, nimble::STEP_N);
    auto deformation_gradient_step_np1_d = model_data.GetDeviceFullTensorIntegrationPointData(
        block_data.id, field_ids.deformation_gradient, nimble::STEP_NP1);
    auto stress_step_n_d = model_data.GetDeviceSymTensorIntegrationPointData(block_data.id, field_ids.stress,
                                                                             nimble::STEP_N);
    auto stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_data.id, field_ids.stress,
                                                                               nimble::STEP_NP1);
    compute_block_stress(block_data, deformation_gradient_step_n_d, deformation_gradient_step_np1_d,
                         stress_step_n_d, stress_step_np1_d, time_n, time_np1);
  }
}

}
