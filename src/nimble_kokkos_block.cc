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

  void Block::Initialize(std::string const & macro_material_parameters) {
    macro_material_parameters_ = macro_material_parameters;
    InstantiateMaterialModel();
    InstantiateElement();
  }

  // void Block::Initialize(std::string const & macro_material_parameters,
  //                        std::map<int, std::string> const & rve_material_parameters,
  //                        GenesisMesh const & rve_mesh,
  //                        std::string rve_boundary_condition_strategy) {
  //   macro_material_parameters_ = macro_material_parameters;
  //   if (macro_material_parameters_ == "none") {
  //     rve_material_parameters_ = rve_material_parameters;
  //     rve_mesh_ = rve_mesh;
  //     rve_boundary_condition_strategy_ = rve_boundary_condition_strategy;
  //   }
  //   InstantiateMaterialModel();
  //   InstantiateElement();
  // }

  void Block::InstantiateMaterialModel() {

    if (macro_material_parameters_ != "none" && rve_material_parameters_.size() == 0) {

      // the first entry in the material parameters string is the material model name
      size_t space_pos = macro_material_parameters_.find(" ");
      std::string name = macro_material_parameters_.substr(0, space_pos);

      // LAME material models are designated with lame_
#ifdef NIMBLE_HAVE_EXTRAS
      bool is_lame_model = false;
      if (name.size() > 5 && name.substr(0,5) == "lame_") {
        is_lame_model = true;
      }

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
      nimble::MaterialParameters material_parameters_struct(material_name, num_material_parameters, material_parameter_names, material_parameter_values);

      if (material_name == "neohookean") {
        material_host_ = std::make_shared<nimble::NeohookeanMaterial>(material_parameters_struct);
        material_device_ = static_cast<nimble::Material*>(Kokkos::kokkos_malloc<>("Material", sizeof(nimble::NeohookeanMaterial)));
        nimble::Material* pointer_that_lives_on_the_stack = material_device_;
        Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
            new (pointer_that_lives_on_the_stack) nimble::NeohookeanMaterial(material_parameters_struct);
          });
      }
      // else if (is_lame_model) {
      //   material_host_ = std::make_shared<LAMEMaterial>(macro_material_parameters_);
      // }
      else if (is_ngp_lame_model) {
        material_host_ = std::make_shared<nimble::NGPLAMEMaterial>(material_parameters_struct);
        material_device_ = static_cast<nimble::Material*>(Kokkos::kokkos_malloc<>("Material", sizeof(nimble::NeohookeanMaterial)));
        nimble::Material* pointer_that_lives_on_the_stack = material_device_;
        Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
            new (pointer_that_lives_on_the_stack) nimble::NGPLAMEMaterial(material_parameters_struct);
          });
      }
      else {
        throw std::logic_error("\nError in Block::InstantiateMaterialModel(), invalid material model name.\n");
      }
    }
    // else if (macro_material_parameters_ == "none" && rve_material_parameters_.size() != 0) {
    //   material_host_ = std::make_shared<RVE>(rve_material_parameters_, rve_mesh_, rve_boundary_condition_strategy_);
    // }
    else if (macro_material_parameters_ != "none" && rve_material_parameters_.size() != 0) {
      throw std::logic_error("\nError:  Assigning both a macroscale material and an RVE material to the same block is currently not supported.\n");
    }
    else{
      throw std::logic_error("\nError in Block::InstantiateMaterialModel(), invalid material parameters\n");
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

  void Block::GetIntegrationPointDataLabelsAndLengths(std::vector< std::pair<std::string, nimble::Length> >& data_labels_and_lengths) {

    // kinematic data
    data_labels_and_lengths.push_back(std::pair<std::string, nimble::Length>("deformation_gradient", nimble::FULL_TENSOR));
    data_labels_and_lengths.push_back(std::pair<std::string, nimble::Length>("stress", nimble::SYMMETRIC_TENSOR));
    // material state data
    int num_data = material_host_->NumStateVariables();
    char label[nimble::MaterialParameters::MAX_MAT_MODEL_STR_LEN];
    for (int i=0 ; i<num_data ; ++i) {
      material_host_->GetStateVariableLabel(i, label);
      data_labels_and_lengths.push_back(std::pair<std::string, nimble::Length>(label, nimble::SCALAR));
    }
  }

  double Block::ComputeCriticalTimeStep(const double * const node_coordinates,
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
          node_coord[node*vector_size + i] = node_coordinates[vector_size*node_id + i];
        }
      }

      double elem_critical_time_step = element_host_->ComputeCharacteristicLength(node_coord) / sound_speed;
      if (elem_critical_time_step < critical_time_step) {
        critical_time_step = elem_critical_time_step;
      }
    }

    return critical_time_step;
  }

  void Block::InitializeElementData(int num_elem_in_block,
                                    std::vector<int> const & elem_global_ids_in_block,
                                    std::vector<int> const & rve_output_global_elem_ids,
                                    std::vector<std::string> const & elem_data_labels,
                                    std::vector<std::string> const & derived_elem_data_labels,
                                    std::vector<double>& elem_data_n,
                                    std::vector<double>& elem_data_np1,
                                    DataManager& data_manager) {

    int num_int_pts = element_host_->NumIntegrationPointsPerElement();
    for (int i_elem = 0 ; i_elem < num_elem_in_block; ++i_elem) {
      int elem_global_id = elem_global_ids_in_block.at(i_elem);
      // bool rve_exodus_output = false;
      // if (std::find(rve_output_global_elem_ids.begin(), rve_output_global_elem_ids.end(), elem_global_id) != rve_output_global_elem_ids.end()) {
      //   rve_exodus_output = true;
      // }
      // for (int i_ipt = 0 ; i_ipt < num_int_pts ; ++i_ipt) {
      //   material_host_->InitializeRVE(elem_global_id, i_ipt+1, data_manager, rve_exodus_output);
      // }
    }

    unsigned int stride = elem_data_labels.size();

    std::map<std::string, double> state_variable_initial_values;
    char label[nimble::MaterialParameters::MAX_MAT_MODEL_STR_LEN];
    for (int i=0 ; i<material_host_->NumStateVariables() ; ++i) {
      material_host_->GetStateVariableLabel(i, label);
      double state_variable_initial_value = material_host_->GetStateVariableInitialValue(i);
      std::string state_variable_label(label);
      state_variable_initial_values[state_variable_label] = state_variable_initial_value;
    }

    for (unsigned int i_variable = 0 ; i_variable < elem_data_labels.size() ; i_variable++) {

      std::string label = elem_data_labels[i_variable];
      std::string ipt_label = nimble::RemoveIntegrationPointPrefix(label);

      double initial_value = 0.0;
      auto it = state_variable_initial_values.find(label);
      if (it != state_variable_initial_values.end()) {
        initial_value = it->second;
      }

      for (int i_elem = 0 ; i_elem < num_elem_in_block; ++i_elem) {

        // override the material's intial values for the special cases
        // of global element id and integration point number
        if (ipt_label == "global_elem_id") {
          int elem_global_id = elem_global_ids_in_block.at(i_elem);
          initial_value = elem_global_id;
        }
        else if(ipt_label == "ipt_id") {
          int integration_point_id = static_cast<double>( nimble::LabelToIntegrationPointNumber(label) );
          initial_value = integration_point_id;
        }
        else if (ipt_label == "deformation_gradient_xx" || ipt_label == "deformation_gradient_yy" || ipt_label == "deformation_gradient_zz") {
          initial_value = 1.0;
        }

        int index = i_elem*stride + i_variable;
        elem_data_n[index] = initial_value;
        elem_data_np1[index] = initial_value;
      }
    }

    DetermineDataOffsets(elem_data_labels, derived_elem_data_labels);
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

  void Block::ComputeDerivedElementData(const double * const reference_coordinates,
                                        const double * const displacement,
                                        int num_elem,
                                        const int * const elem_conn,
                                        int num_elem_data,
                                        std::vector<double> const & elem_data_np1,
                                        int num_derived_elem_data,
                                        std::vector< std::vector<double> >& derived_elem_data) {

    // Currently supported derived element data are:
    // 1) Element volume
    // 2) Volume-averaged element data

    if (num_derived_elem_data == 0){
      return;
    }
    if (derived_elem_data.size() != num_derived_elem_data) {
      derived_elem_data.resize(num_derived_elem_data);
    }
    for (auto & vec : derived_elem_data) {
      if (vec.size() != num_elem) {
        vec.resize(num_elem);
      }
    }

    int dim = element_host_->Dim();
    int vector_size = LengthToInt(nimble::VECTOR, dim);
    int num_node_per_elem = element_host_->NumNodesPerElement();
    int num_int_pt_per_elem = element_host_->NumIntegrationPointsPerElement();
    int num_data_per_int_pt = num_elem_data / num_int_pt_per_elem;

    // Containers that will be passed to an individual element to compute element volume and volume-averaged quantities
    double ref_coord[vector_size*num_node_per_elem];
    double cur_coord[vector_size*num_node_per_elem];
    unsigned int num_vol_ave_data = vol_ave_offsets_.size() / num_int_pt_per_elem;
    double int_pt_quantities[num_vol_ave_data*num_int_pt_per_elem];
    double volume;
    double vol_ave_quantities[num_vol_ave_data];

    // Loop over each element and compute element volume and volume-averaged quantities
    for (int i_elem = 0 ; i_elem < num_elem ; i_elem++) {

      // Load per element data containers
      for (int node = 0 ; node < num_node_per_elem ; node++) {
        int node_id = elem_conn[i_elem*num_node_per_elem + node];
        for (int i = 0 ; i < vector_size ; i++) {
          ref_coord[node*vector_size + i] = reference_coordinates[vector_size*node_id + i];
          cur_coord[node*vector_size + i] = reference_coordinates[vector_size*node_id + i] + displacement[vector_size*node_id + i];
        }
      }
      for (int i = 0 ; i < num_vol_ave_data*num_int_pt_per_elem ; ++i) {
        int_pt_quantities[i] = elem_data_np1[i_elem*num_data_per_int_pt*num_int_pt_per_elem + vol_ave_offsets_.at(i)];
      }

      // Compute element volume and volume-averaged quantities
      element_host_->ComputeVolumeAverage(ref_coord,
                                          cur_coord,
                                          num_vol_ave_data,
                                          &int_pt_quantities[0],
                                          volume,
                                          vol_ave_quantities);

      // Copy the results into the global containers for derived data
      for (int i = 0 ; i < num_vol_ave_data ; i++) {
        derived_elem_data.at(vol_ave_index_to_derived_data_index_[i]).at(i_elem) = vol_ave_quantities[i];
      }
      if (vol_ave_volume_offset_ != -1) {
        derived_elem_data.at(vol_ave_volume_offset_).at(i_elem) = volume;
      }
    }
  }

  void Block::DetermineDataOffsets(std::vector<std::string> const & elem_data_labels,
                                   std::vector<std::string> const & derived_elem_data_labels)
  {
    int dim = element_host_->Dim();
    int num_int_pt_per_elem = element_host_->NumIntegrationPointsPerElement();
    int num_data_per_element = static_cast<int>( elem_data_labels.size() );
    int num_data_per_int_pt = num_data_per_element / num_int_pt_per_elem;
    int full_tensor_size = LengthToInt(nimble::FULL_TENSOR, dim);
    int sym_tensor_size = LengthToInt(nimble::SYMMETRIC_TENSOR, dim);

    // Determine the offset for the derformation gradient in the element data container
    def_grad_offset_.clear();
    for (int i_ipt = 0 ; i_ipt < num_int_pt_per_elem ; i_ipt++) {
      for (int i_component = 0 ; i_component < full_tensor_size ; i_component++) {
        std::string def_grad_label = GetComponentLabel("deformation_gradient", nimble::FULL_TENSOR, dim, i_component, i_ipt+1);
        for (unsigned int i_label = 0 ; i_label < num_data_per_element ; i_label++) {
          if (elem_data_labels.at(i_label) == def_grad_label) {
            def_grad_offset_.push_back(i_label);
          }
        }
      }
    }
    if (def_grad_offset_.size() != full_tensor_size*num_int_pt_per_elem) {
      throw std::logic_error("\nError in Block::DetermineDataOffsets(), failed to index def_grad into global data.\n");
    }

    // Determine the offset for the stress in the element data container
    stress_offset_.clear();
    for (int i_ipt = 0 ; i_ipt < num_int_pt_per_elem ; i_ipt++) {
      for (int i_component = 0 ; i_component < sym_tensor_size ; i_component++) {
        std::string stress_label = GetComponentLabel("stress", nimble::SYMMETRIC_TENSOR, dim, i_component, i_ipt+1);
        for (unsigned int i_label = 0 ; i_label < num_data_per_element ; i_label++) {
          if (elem_data_labels.at(i_label) == stress_label) {
            stress_offset_.push_back(i_label);
          }
        }
      }
    }
    if (stress_offset_.size() != sym_tensor_size*num_int_pt_per_elem) {
      throw std::logic_error("\nError in Block::DetermineDataOffsets(), failed to index stress into global data.\n");
    }

    // Determine the offset for the state data in the element data container
    int num_state_data = material_host_->NumStateVariables();
    std::vector<std::string> state_data_labels(num_state_data);
    for (int i=0 ; i<num_state_data ; ++i) {
      char label[nimble::MaterialParameters::MAX_MAT_MODEL_STR_LEN];
      material_host_->GetStateVariableLabel(i, label);
      state_data_labels[i] = label;
    }
    state_data_offset_.clear();
    for (int i_ipt = 0 ; i_ipt < num_int_pt_per_elem ; i_ipt++) {
      for (int i_state_data = 0 ; i_state_data < num_state_data ; i_state_data++) {
        std::string raw_label = state_data_labels.at(i_state_data);
        std::string state_data_label = GetComponentLabel(raw_label, nimble::SCALAR, dim, 0, i_ipt+1);
        for (unsigned int i_label = 0 ; i_label < num_data_per_element ; i_label++) {
          if (elem_data_labels.at(i_label) == state_data_label) {
            state_data_offset_.push_back(i_label);
          }
        }
      }
    }
    if (state_data_offset_.size() != num_state_data*num_int_pt_per_elem) {
      throw std::logic_error("\nError in Block::DetermineDataOffsets(), failed to index state data into global data.\n");
    }

    // Determine which requested derived data are volume-averaged data
    // and set up some bookkeeping
    std::vector<std::string> vol_ave_labels;
    vol_ave_volume_offset_ = -1;
    for (unsigned int i_derived_data = 0 ; i_derived_data < derived_elem_data_labels.size() ; i_derived_data++) {
      std::string const & label = derived_elem_data_labels[i_derived_data];
      if (label == "volume") {
        vol_ave_volume_offset_ = i_derived_data;
      }
      else {
        vol_ave_index_to_derived_data_index_[vol_ave_labels.size()] = i_derived_data;
        vol_ave_labels.push_back(label);
      }
    }
    unsigned int num_vol_ave_data = vol_ave_labels.size();
    vol_ave_offsets_.resize(num_vol_ave_data*num_int_pt_per_elem);
    // Bookkeeping to allow for loading the per element containers
    for (unsigned int i_vol_ave_label = 0 ; i_vol_ave_label < num_vol_ave_data ; ++i_vol_ave_label) {
      std::string const & vol_ave_label = vol_ave_labels[i_vol_ave_label];
      for (unsigned int i_elem_data_label = 0 ; i_elem_data_label < num_data_per_int_pt ; ++i_elem_data_label) {
        std::string const & elem_data_label = nimble::RemoveIntegrationPointPrefix(elem_data_labels[i_elem_data_label]);
        if (vol_ave_label == elem_data_label) {
          for (int i_int_pt = 0 ; i_int_pt < num_int_pt_per_elem ; ++i_int_pt) {
            vol_ave_offsets_[i_vol_ave_label + i_int_pt*num_vol_ave_data] = i_elem_data_label + i_int_pt*num_data_per_int_pt;
          }
        }
      }
    }
  }

} // namespace nimble
