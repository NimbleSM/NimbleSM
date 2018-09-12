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

#include "nimble_kokkos_block.h"
#include "nimble_linear_solver.h"
//#include "nimble_data_manager.h"
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <algorithm>

#include "nimble_utils.h"

namespace nimble_kokkos {

  void Block::Initialize(std::string const & macro_material_parameters,
                         int num_elements) {
    macro_material_parameters_ = macro_material_parameters;
    InstantiateElement();
    int num_material_points = num_elements * element_host_->NumIntegrationPointsPerElement();
    InstantiateMaterialModel(num_material_points);
  }

  void Block::InstantiateMaterialModel(int num_material_points) {

    // the first entry in the material parameters string is the material model name
    size_t space_pos = macro_material_parameters_.find(" ");
    std::string name = macro_material_parameters_.substr(0, space_pos);

#ifdef NIMBLE_HAVE_EXTRAS
    // NGP LAME material models are designated with ngp_lame_
    bool is_ngp_lame_model = false;
    if (name.size() > 9 && name.substr(0,9) == "ngp_lame_") {
      is_ngp_lame_model = true;
    }
#endif

    char material_name[nimble::MaterialParameters::MAX_MAT_MODEL_STR_LEN];
    int num_material_parameters;
    char material_parameter_names[nimble::MaterialParameters::MAX_NUM_MAT_PARAM][nimble::MaterialParameters::MAX_MAT_MODEL_STR_LEN];
    double material_parameter_values[nimble::MaterialParameters::MAX_NUM_MAT_PARAM];
    nimble::ParseMaterialParametersString(macro_material_parameters_.c_str(), material_name, num_material_parameters, material_parameter_names, material_parameter_values);
    nimble::MaterialParameters material_parameters_struct(material_name, num_material_parameters, material_parameter_names, material_parameter_values, num_material_points);
    if (nimble::StringsAreEqual(material_name, "neohookean")) {
      material_host_ = std::make_shared<nimble::NeohookeanMaterial>(material_parameters_struct);
      material_device_ = static_cast<nimble::Material*>(Kokkos::kokkos_malloc<>("Material", sizeof(nimble::NeohookeanMaterial)));
      nimble::Material* pointer_that_lives_on_the_stack = material_device_;
      Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
          new (pointer_that_lives_on_the_stack) nimble::NeohookeanMaterial(material_parameters_struct);
        });
    }
#ifdef NIMBLE_HAVE_EXTRAS
    else if (is_ngp_lame_model) {
      std::cout << "DEBUGGING is_ngp_lame_model" << std::endl;
      material_host_ = std::make_shared<nimble::NGPLAMEMaterial>(material_parameters_struct);
      material_device_ = static_cast<nimble::Material*>(Kokkos::kokkos_malloc<>("Material", sizeof(nimble::NGPLAMEMaterial)));
      nimble::Material* pointer_that_lives_on_the_stack = material_device_;
      Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
          new (pointer_that_lives_on_the_stack) nimble::NGPLAMEMaterial(material_parameters_struct);
        });
    }
#endif
    else {
      throw std::logic_error("\nError in Block::InstantiateMaterialModel(), invalid material model name.\n");
    }
  }

  void Block::InstantiateElement() {

    // instantiate the element on the host (eventually we won't need this)
    element_host_ = std::make_shared<nimble::HexElement>();

    // instantiate the element on the device
    element_device_ = static_cast<nimble::Element*>(Kokkos::kokkos_malloc<>("Element", sizeof(nimble::HexElement)));
    nimble::Element* pointer_that_lives_on_the_stack = element_device_;
    Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
        new (pointer_that_lives_on_the_stack) nimble::HexElement();
      });
  }

  double Block::ComputeCriticalTimeStep(const double * const node_reference_coordinates,
                                        const double * const node_displacements,
                                        int num_elem,
                                        const int * const elem_conn) const {
    int dim = element_host_->Dim();
    int num_node_per_elem = element_host_->NumNodesPerElement();
    double sound_speed = std::sqrt( GetBulkModulus() / GetDensity() );
    double critical_time_step = std::numeric_limits<double>::max();

    int vector_size = 0;
    if (dim == 2) {
      vector_size = 2;
    }
    else if (dim == 3) {
      vector_size = 3;
    }
    double node_coord[vector_size*num_node_per_elem];

    for (int elem=0 ; elem<num_elem ; elem++) {

      for (int node=0 ; node<num_node_per_elem ; node++) {
        int node_id = elem_conn[elem*num_node_per_elem + node];
        for (int i=0 ; i<vector_size ; i++) {
          node_coord[node*vector_size + i] = node_reference_coordinates[vector_size*node_id + i] + node_displacements[vector_size*node_id + i];
        }
      }

      double elem_critical_time_step = element_host_->ComputeCharacteristicLength(node_coord) / sound_speed;
      if (elem_critical_time_step < critical_time_step) {
        critical_time_step = elem_critical_time_step;
      }
    }

    return critical_time_step;
  }


  template <typename MatT>
  void Block::ComputeTangentStiffnessMatrix(int num_global_unknowns,
                                            const double * const reference_coordinates,
                                            const double * const displacement,
                                            int num_elem,
                                            const int * const elem_conn,
                                            const int * const global_node_ids,
                                            MatT & tangent_stiffness) const {

    int dim = element_host_->Dim();
    int num_node_per_elem = element_host_->NumNodesPerElement();
    int num_int_pt_per_elem = element_host_->NumIntegrationPointsPerElement();

    int vector_size = LengthToInt(nimble::VECTOR, dim);
    int full_tensor_size = LengthToInt(nimble::FULL_TENSOR, dim);
    int sym_tensor_size = LengthToInt(nimble::SYMMETRIC_TENSOR, dim);

    double cur_coord[vector_size * num_node_per_elem];
    double material_tangent[6 * 6 * num_int_pt_per_elem]; // correct for 3D, overkill for 2D
    double element_tangent[num_node_per_elem * vector_size * num_node_per_elem * vector_size];

    for (int elem=0 ; elem<num_elem ; elem++) {

      for (int node=0 ; node<num_node_per_elem ; node++) {
        int node_id = elem_conn[elem*num_node_per_elem + node];
        for (int i=0 ; i<vector_size ; i++) {
          cur_coord[node*vector_size + i] = reference_coordinates[vector_size*node_id + i] + displacement[vector_size*node_id + i];
        }
      }

      material_host_->GetTangent(num_int_pt_per_elem, material_tangent);

      element_host_->ComputeTangent(cur_coord, material_tangent, element_tangent);

      for (int row=0 ; row<num_node_per_elem ; row++) {
        int global_row_node = global_node_ids[ elem_conn[elem*num_node_per_elem + row] ];
        for (int col=0 ; col<num_node_per_elem ; col++) {
          int global_col_node = global_node_ids[ elem_conn[elem*num_node_per_elem + col] ];
          for (int i=0 ; i<vector_size ; i++) {
            for (int j=0 ; j<vector_size ; j++) {
              int local_row_index = row * vector_size + i;
              int local_col_index = col * vector_size + j;
              int global_row_index = global_row_node * vector_size + i;
              int global_col_index = global_col_node * vector_size + j;
              double value = element_tangent[local_row_index*num_node_per_elem*vector_size + local_col_index];
              tangent_stiffness(global_row_index, global_col_index) += value;
              //tangent_stiffness.sumIntoValue(global_row_index, global_col_index, value);
            }
          }
        }
      }
    }
  }

  // explicit template instantiation for MatrixContainer and CRSMatrixContainer matrix types
  template
  void Block::ComputeTangentStiffnessMatrix<nimble::MatrixContainer>(int num_global_unknowns,
                                                                     const double * const reference_coordinates,
                                                                     const double * const displacement,
                                                                     int num_elem,
                                                                     const int * const elem_conn,
                                                                     const int * const global_node_ids,
                                                                     nimble::MatrixContainer & tangent_stiffness) const ;
  template
  void Block::ComputeTangentStiffnessMatrix<nimble::CRSMatrixContainer>(int num_global_unknowns,
                                                                        const double * const reference_coordinates,
                                                                        const double * const displacement,
                                                                        int num_elem,
                                                                        const int * const elem_conn,
                                                                        const int * const global_node_ids,
                                                                        nimble::CRSMatrixContainer & tangent_stiffness) const ;

} // namespace nimble
