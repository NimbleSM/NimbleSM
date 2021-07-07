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

#ifdef NIMBLE_HAVE_UQ

#include "nimble_uq_block.h"

#include <memory>
#include <string>
#include <vector>

#include "nimble_data_manager.h"
#include "nimble_element.h"

namespace nimble_uq {

void
Block::ComputeInternalForce(
    const double*                   reference_coordinates,
    const double*                   displacement,
    const double*                   velocity,
    double*                         internal_force,
    double                          time_previous,
    double                          time_current,
    int                             num_elem,
    const int*                      elem_conn,
    const int*                      elem_global_ids,
    std::vector<std::string> const& elem_data_labels,
    std::vector<double> const&      elem_data_n,
    std::vector<double>&            elem_data_np1,
    nimble::DataManager&            data_manager,
    bool                            is_output_step,
    const bool&                     is_off_nominal,
    std::vector<double> const&      uq_params_this_sample,
    bool                            compute_stress_only) const
{
  if (!is_off_nominal) {
    nimble::Block::ComputeInternalForce(
        reference_coordinates,
        displacement,
        velocity,
        internal_force,
        time_previous,
        time_current,
        num_elem,
        elem_conn,
        elem_global_ids,
        elem_data_labels,
        elem_data_n,
        elem_data_np1,
        data_manager,
        is_output_step,
        compute_stress_only);
    return;
  }
  //
  //--- Treat case of off-nominal solutions
  //
  int dim                 = element_->Dim();
  int num_node_per_elem   = element_->NumNodesPerElement();
  int num_int_pt_per_elem = element_->NumIntegrationPointsPerElement();
  int num_element_data    = static_cast<int>(elem_data_labels.size());

  int vector_size      = LengthToInt(nimble::VECTOR, dim);
  int full_tensor_size = LengthToInt(nimble::FULL_TENSOR, dim);
  int sym_tensor_size  = LengthToInt(nimble::SYMMETRIC_TENSOR, dim);

  double ref_coord[vector_size * num_node_per_elem];
  double cur_coord[vector_size * num_node_per_elem];
  double def_grad_n[full_tensor_size * num_int_pt_per_elem];
  double def_grad_np1[full_tensor_size * num_int_pt_per_elem];
  double cauchy_stress_n[sym_tensor_size * num_int_pt_per_elem];
  double cauchy_stress_np1[sym_tensor_size * num_int_pt_per_elem];
  double force[vector_size * num_node_per_elem];

  double*             state_data_n   = nullptr;
  double*             state_data_np1 = nullptr;
  std::vector<double> state_data_n_vec;
  std::vector<double> state_data_np1_vec;
  int                 num_state_data = material_->NumStateVariables();
  if (num_state_data > 0) {
    state_data_n_vec.resize(num_state_data * num_int_pt_per_elem);
    state_data_np1_vec.resize(num_state_data * num_int_pt_per_elem);
    state_data_n   = state_data_n_vec.data();
    state_data_np1 = state_data_np1_vec.data();
  }

  std::vector<double> identity(full_tensor_size, 0.0);
  for (int i = 0; i < dim; i++) { identity[i] = 1.0; }

  for (int elem = 0; elem < num_elem; elem++) {
    for (int node = 0; node < num_node_per_elem; node++) {
      int node_id = elem_conn[elem * num_node_per_elem + node];
      for (int i = 0; i < vector_size; i++) {
        ref_coord[node * vector_size + i] = reference_coordinates[vector_size * node_id + i];
        cur_coord[node * vector_size + i] =
            reference_coordinates[vector_size * node_id + i] + displacement[vector_size * node_id + i];
      }
    }

    element_->ComputeDeformationGradients(ref_coord, cur_coord, def_grad_np1);

    material_->GetOffNominalStress(
        uq_params_this_sample[bulk_modulus_uq_index_],
        uq_params_this_sample[shear_modulus_uq_index_],
        num_int_pt_per_elem,
        def_grad_np1,
        cauchy_stress_np1);

    if (!compute_stress_only) {
      // Compute element contribution to nodal force
      element_->ComputeNodalForces(cur_coord, cauchy_stress_np1, force);

      // Copy internal force to global containers
      for (int node = 0; node < num_node_per_elem; node++) {
        int node_id = elem_conn[elem * num_node_per_elem + node];
        for (int i = 0; i < vector_size; i++) {
          internal_force[vector_size * node_id + i] += force[node * vector_size + i];
        }
      }
    }

  }  // for (int elem = 0; elem < num_elem; elem++)
}

}  // namespace nimble_uq

#endif
