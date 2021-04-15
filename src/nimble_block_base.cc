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

#include "nimble_block_base.h"
#include "nimble_element.h"
#include "nimble_linear_solver.h"

namespace nimble {

double
BlockBase::ComputeCriticalTimeStep(
    const nimble::Viewify<2>& node_reference_coordinates,
    const nimble::Viewify<2>& node_displacements,
    int                       num_elem,
    const int* const          elem_conn) const
{
  int    dim                = element_->Dim();
  int    num_node_per_elem  = element_->NumNodesPerElement();
  double sound_speed        = std::sqrt(GetBulkModulus() / GetDensity());
  double critical_time_step = std::numeric_limits<double>::max();

  int vector_size = 0;
  if (dim == 2) {
    vector_size = 2;
  } else if (dim == 3) {
    vector_size = 3;
  }
  double node_coord[vector_size * num_node_per_elem];

  for (int elem = 0; elem < num_elem; elem++) {
    for (int node = 0; node < num_node_per_elem; node++) {
      int node_id = elem_conn[elem * num_node_per_elem + node];
      for (int i = 0; i < vector_size; i++) {
        node_coord[node * vector_size + i] = node_reference_coordinates(node_id, i) + node_displacements(node_id, i);
      }
    }

    double elem_critical_time_step = element_->ComputeCharacteristicLength(node_coord) / sound_speed;
    if (elem_critical_time_step < critical_time_step) { critical_time_step = elem_critical_time_step; }
  }

  return critical_time_step;
}

template <typename MatT>
void
BlockBase::ComputeTangentStiffnessMatrix(
    int                 num_global_unknowns,
    const double* const reference_coordinates,
    const double* const displacement,
    int                 num_elem,
    const int* const    elem_conn,
    const int* const    global_node_ids,
    MatT&               tangent_stiffness) const
{
  int dim                 = element_->Dim();
  int num_node_per_elem   = element_->NumNodesPerElement();
  int num_int_pt_per_elem = element_->NumIntegrationPointsPerElement();

  int vector_size      = LengthToInt(nimble::VECTOR, dim);
  int full_tensor_size = LengthToInt(nimble::FULL_TENSOR, dim);
  int sym_tensor_size  = LengthToInt(nimble::SYMMETRIC_TENSOR, dim);

  double cur_coord[vector_size * num_node_per_elem];
  double material_tangent[6 * 6 * num_int_pt_per_elem];  // correct for 3D,
  // overkill for 2D
  double element_tangent[num_node_per_elem * vector_size * num_node_per_elem * vector_size];

  for (int elem = 0; elem < num_elem; elem++) {
    for (int node = 0; node < num_node_per_elem; node++) {
      int node_id = elem_conn[elem * num_node_per_elem + node];
      for (int i = 0; i < vector_size; i++) {
        cur_coord[node * vector_size + i] =
            reference_coordinates[vector_size * node_id + i] + displacement[vector_size * node_id + i];
      }
    }

    material_->GetTangent(num_int_pt_per_elem, material_tangent);

    element_->ComputeTangent(cur_coord, material_tangent, element_tangent);

    for (int row = 0; row < num_node_per_elem; row++) {
      int global_row_node = global_node_ids[elem_conn[elem * num_node_per_elem + row]];
      for (int col = 0; col < num_node_per_elem; col++) {
        int global_col_node = global_node_ids[elem_conn[elem * num_node_per_elem + col]];
        for (int i = 0; i < vector_size; i++) {
          for (int j = 0; j < vector_size; j++) {
            int    local_row_index  = row * vector_size + i;
            int    local_col_index  = col * vector_size + j;
            int    global_row_index = global_row_node * vector_size + i;
            int    global_col_index = global_col_node * vector_size + j;
            double value = element_tangent[local_row_index * num_node_per_elem * vector_size + local_col_index];
            tangent_stiffness(global_row_index, global_col_index) += value;
            // tangent_stiffness.sumIntoValue(global_row_index,
            // global_col_index, value);
          }
        }
      }
    }
  }
}

// explicit template instantiation for MatrixContainer and CRSMatrixContainer
// matrix types
template void
BlockBase::ComputeTangentStiffnessMatrix<MatrixContainer>(
    int                 num_global_unknowns,
    const double* const reference_coordinates,
    const double* const displacement,
    int                 num_elem,
    const int* const    elem_conn,
    const int* const    global_node_ids,
    MatrixContainer&    tangent_stiffness) const;
template void
BlockBase::ComputeTangentStiffnessMatrix<CRSMatrixContainer>(
    int                 num_global_unknowns,
    const double* const reference_coordinates,
    const double* const displacement,
    int                 num_elem,
    const int* const    elem_conn,
    const int* const    global_node_ids,
    CRSMatrixContainer& tangent_stiffness) const;

}  // namespace nimble
