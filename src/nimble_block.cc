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

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>
#include <nimble_block.h>
#include <nimble_data_utils.h>
#include <nimble_linear_solver.h>
#include <nimble_material.h>
#include <nimble_material_factory.h>

#ifndef NIMBLE_HAVE_KOKKOS
#include "nimble_rve.h"
#endif

#ifdef NIMBLE_HAVE_UQ
  #include "nimble_uq.h"
#endif

namespace nimble {

  void Block::Initialize(std::string const & macro_material_parameters,
                         MaterialFactory& factory) {
    macro_material_parameters_ = macro_material_parameters;
    InstantiateMaterialModel(factory);
    InstantiateElement();
  }

  void Block::Initialize(std::string const & macro_material_parameters,
                         std::map<int, std::string> const & rve_material_parameters,
                         GenesisMesh const & rve_mesh,
                         std::string rve_boundary_condition_strategy,
                         MaterialFactory& factory) {
    macro_material_parameters_ = macro_material_parameters;
    if (macro_material_parameters_ == "none") {
      rve_material_parameters_ = rve_material_parameters;
      rve_mesh_ = rve_mesh;
      rve_boundary_condition_strategy_ = rve_boundary_condition_strategy;
    }
    InstantiateMaterialModel(factory);
    InstantiateElement();
  }

  void Block::InstantiateMaterialModel(MaterialFactory& factory) {

    if (macro_material_parameters_ != "none" && rve_material_parameters_.size() == 0) {
      factory.parse_and_create(macro_material_parameters_);
      material_ = factory.get_material();
    }
#ifndef NIMBLE_HAVE_KOKKOS
    else if (macro_material_parameters_ == "none" && rve_material_parameters_.size() != 0) {
      material_ = std::make_shared<RVE>(rve_material_parameters_, rve_mesh_, rve_boundary_condition_strategy_);
    }
#endif
    else if (macro_material_parameters_ != "none" && rve_material_parameters_.size() != 0) {
      throw std::logic_error("\nError:  Assigning both a macroscale material and an RVE material to the same block is currently not supported.\n");
    }
    else{
      throw std::logic_error("\nError in Block::InstantiateMaterialModel(), invalid material parameters\n");
    }
  }

  void Block::InstantiateElement() {
    element_ = std::make_shared<HexElement>();
  }

  void Block::GetDataLabelsAndLengths(std::vector< std::pair<std::string, Length> >& data_labels_and_lengths) {
    int num_int_pts = element_->NumIntegrationPointsPerElement();
    std::vector< std::pair<std::string, Length> > mat_model_data_labels_and_lengths;
    // kinematic data
    mat_model_data_labels_and_lengths.push_back(std::pair<std::string, Length>("deformation_gradient", FULL_TENSOR));
    mat_model_data_labels_and_lengths.push_back(std::pair<std::string, Length>("stress", SYMMETRIC_TENSOR));
    // material state data
    int num_data = material_->NumStateVariables();
    char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN];
    for (int i=0 ; i<num_data ; ++i) {
      material_->GetStateVariableLabel(i, label);
      mat_model_data_labels_and_lengths.push_back(std::pair<std::string, Length>(label, SCALAR));
    }
    for (int i_ipt = 0 ; i_ipt < num_int_pts ; i_ipt++) {
      for (unsigned int i_label = 0 ; i_label < mat_model_data_labels_and_lengths.size() ; i_label++) {
        std::string material_label = mat_model_data_labels_and_lengths[i_label].first;
        Length length = mat_model_data_labels_and_lengths[i_label].second;
        std::string label = AddIntegrationPointPrefix(material_label, i_ipt+1);
        std::pair<std::string, Length> label_and_length(label, length);
        data_labels_and_lengths.push_back(label_and_length);
      }
    }
  }

  void Block::ComputeLumpedMassMatrix(const double * const reference_coordinates,
                                      int num_elem,
                                      const int * const elem_conn,
                                      double* lumped_mass) const {

    int dim = element_->Dim();
    int num_node_per_elem = element_->NumNodesPerElement();
    double density = GetDensity();

    int vector_size = 0;
    if (dim == 2) {
      vector_size = 2;
    }
    else if (dim == 3) {
      vector_size = 3;
    }

    double ref_coord[vector_size*num_node_per_elem];
    double mass[num_node_per_elem];

    for (int elem=0 ; elem<num_elem ; elem++) {

      for (int node=0 ; node<num_node_per_elem ; node++) {
        int node_id = elem_conn[elem*num_node_per_elem + node];
        for (int i=0 ; i<vector_size ; i++) {
          ref_coord[node*vector_size + i] = reference_coordinates[vector_size*node_id + i];
        }
      }

      element_->ComputeLumpedMass(density, ref_coord, mass);

      for (int node=0 ; node<num_node_per_elem ; node++) {
        int node_id = elem_conn[elem*num_node_per_elem + node];
        lumped_mass[node_id] += mass[node];
      }
    }
  }

  double Block::ComputeCriticalTimeStep(const double * const node_reference_coordinates,
                                        const double * const node_displacements,
                                        int num_elem,
                                        const int * const elem_conn) const {
    int dim = element_->Dim();
    int num_node_per_elem = element_->NumNodesPerElement();
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

      double elem_critical_time_step = element_->ComputeCharacteristicLength(node_coord) / sound_speed;
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
                                    MaterialFactory& material_factory,
                                    DataManager& data_manager) {

    int num_int_pts = element_->NumIntegrationPointsPerElement();
#ifndef NIMBLE_HAVE_KOKKOS
    for (int i_elem = 0 ; i_elem < num_elem_in_block; ++i_elem) {
      int elem_global_id = elem_global_ids_in_block.at(i_elem);
      bool rve_exodus_output = false;
      if (std::find(rve_output_global_elem_ids.begin(), rve_output_global_elem_ids.end(), elem_global_id) != rve_output_global_elem_ids.end()) {
        rve_exodus_output = true;
      }
      for (int i_ipt = 0 ; i_ipt < num_int_pts ; ++i_ipt) {
        material_->InitializeRVE(elem_global_id, i_ipt+1, data_manager, rve_exodus_output, material_factory);
      }
    }
#endif

    unsigned int stride = elem_data_labels.size();

    std::map<std::string, double> state_variable_initial_values;
    char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN];
    for (int i=0 ; i<material_->NumStateVariables() ; ++i) {
      material_->GetStateVariableLabel(i, label);
      double state_variable_initial_value = material_->GetStateVariableInitialValue(i);
      std::string state_variable_label(label);
      state_variable_initial_values[state_variable_label] = state_variable_initial_value;
    }

    for (unsigned int i_variable = 0 ; i_variable < elem_data_labels.size() ; i_variable++) {

      std::string label = elem_data_labels[i_variable];
      std::string ipt_label = RemoveIntegrationPointPrefix(label);

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
          int integration_point_id = static_cast<double>( LabelToIntegrationPointNumber(label) );
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

  void Block::ComputeInternalForce(const double * const reference_coordinates,
                                   const double * const displacement,
                                   const double * const velocity,
                                   const double * const rve_macroscale_deformation_gradient,
                                   double * const internal_force,
                                   double time_previous,
                                   double time_current,
                                   int num_elem,
                                   const int * const elem_conn,
                                   const int * const elem_global_ids,
                                   std::vector<std::string> const & elem_data_labels,
                                   std::vector<double> const & elem_data_n,
                                   std::vector<double> & elem_data_np1,
                                   DataManager& data_manager,
                                   bool is_output_step,
#ifdef NIMBLE_HAVE_UQ
                                   const bool & uq_enabled,
                                   const UqParameters* uq_params,
                                   const int & num_uq_samples,
                                   const std::vector<double*> offnominal_displacements,
                                   const std::vector<double*> offnominal_internal_forces,
                                   const std::vector<double*> displacement_sensitivities,
#endif
                                   bool compute_stress_only) const {

    int dim = element_->Dim();
    int num_node_per_elem = element_->NumNodesPerElement();
    int num_int_pt_per_elem = element_->NumIntegrationPointsPerElement();
    int num_element_data = static_cast<int>( elem_data_labels.size() );

    int vector_size = LengthToInt(VECTOR, dim);
    int full_tensor_size = LengthToInt(FULL_TENSOR, dim);
    int sym_tensor_size = LengthToInt(SYMMETRIC_TENSOR, dim);

    double ref_coord[vector_size*num_node_per_elem];
    double cur_coord[vector_size*num_node_per_elem];
    double def_grad_n[full_tensor_size*num_int_pt_per_elem];
    double def_grad_np1[full_tensor_size*num_int_pt_per_elem];
    double cauchy_stress_n[sym_tensor_size*num_int_pt_per_elem];
    double cauchy_stress_np1[sym_tensor_size*num_int_pt_per_elem];
    double force[vector_size*num_node_per_elem];

#ifdef NIMBLE_HAVE_UQ
    int num_uq_params = 0; std::vector<std::vector<double>> def_grad_sensitivities(0);
    std::vector<double> disp_sens_this_param(0), def_grad_this_sample(0);
    double cauchy_stress_this_sample[sym_tensor_size*num_int_pt_per_elem];
    if(uq_enabled) {
      num_uq_params = uq_params->GetNumParams();
      def_grad_sensitivities.resize(num_uq_params);
      for(int np=0; np<num_uq_params; np++) def_grad_sensitivities[np].resize(full_tensor_size*num_int_pt_per_elem, 0.0);
      disp_sens_this_param.resize(vector_size*num_node_per_elem, 0.0);
      def_grad_this_sample.resize(full_tensor_size*num_int_pt_per_elem, 0.0);
    }
#endif

    double* state_data_n(0);
    double* state_data_np1(0);
    std::vector<double> state_data_n_vec;
    std::vector<double> state_data_np1_vec;
    int num_state_data = material_->NumStateVariables();
    if (num_state_data > 0) {
      state_data_n_vec.resize(num_state_data*num_int_pt_per_elem);
      state_data_np1_vec.resize(num_state_data*num_int_pt_per_elem);
      state_data_n = state_data_n_vec.data();
      state_data_np1 = state_data_np1_vec.data();
    }

    std::vector<double> identity(full_tensor_size, 0.0);
    for (int i=0 ; i<dim ; i++) {
      identity[i] = 1.0;
    }

    for (int elem=0 ; elem<num_elem ; elem++) {

      for (int node=0 ; node<num_node_per_elem ; node++) {
        int node_id = elem_conn[elem*num_node_per_elem + node];
        for (int i=0 ; i<vector_size ; i++) {
          ref_coord[node*vector_size + i] = reference_coordinates[vector_size*node_id + i];
          cur_coord[node*vector_size + i] = reference_coordinates[vector_size*node_id + i] + displacement[vector_size*node_id + i];
        }
      }

      element_->ComputeDeformationGradients(ref_coord,
                                            cur_coord,
                                            def_grad_np1);

#ifdef NIMBLE_HAVE_UQ
      if(uq_enabled) {
        //Loop over num_uq_params, compute deformation_gradient_sensitivity
        //"sensitivity" implies derivative w.r.t. uncertain parameter \lambda
        //"gradient" implies derivative w.r.t reference coordinate \Xb
        //IMPORTANT DETAIL: deformation_gradient_sensitivity = displacement_sensitivity_gradient (Go Figure)
        for(int param=0; param<num_uq_params; param++){
          //Copy from global nodal vector to vector for nodes of this element
          for (int node=0 ; node<num_node_per_elem ; node++) {
            int node_id = elem_conn[elem*num_node_per_elem + node];
            for (int i=0 ; i<vector_size ; i++) {
              disp_sens_this_param[node*vector_size + i] = displacement_sensitivities[param][vector_size*node_id + i];
            }
          }
          element_->ComputeDeformationGradients(ref_coord,
                                                &(disp_sens_this_param[0]),
                                                &(def_grad_sensitivities[param][0]) );
        }
      }
#endif

      // Add in the macroscale displacement gradient (nonzero only for multiscale RVE problems)
      for (int i_ipt = 0 ; i_ipt < num_int_pt_per_elem ; i_ipt++) {
        for (int i_component = 0 ; i_component < full_tensor_size ; i_component++) {
          def_grad_np1[i_ipt*full_tensor_size + i_component] += rve_macroscale_deformation_gradient[i_component] - identity[i_component];
        }
      }

      // Copy data from the global data containers
      for (int i_ipt = 0 ; i_ipt < num_int_pt_per_elem ; i_ipt++) {
        for (int i_component = 0 ; i_component < full_tensor_size ; i_component++) {
          int def_grad_offset = def_grad_offset_.at(i_ipt*full_tensor_size + i_component);
          def_grad_n[i_ipt*full_tensor_size + i_component] = elem_data_n[elem*num_element_data + def_grad_offset];
        }
        for (int i_component = 0 ; i_component < sym_tensor_size ; i_component++) {
          int stress_offset = stress_offset_.at(i_ipt*sym_tensor_size + i_component);
          cauchy_stress_n[i_ipt*sym_tensor_size + i_component] = elem_data_n[elem*num_element_data + stress_offset];
        }
        for (int i_component = 0 ; i_component < num_state_data ; i_component++) {
          int state_data_offset = state_data_offset_.at(i_ipt*num_state_data + i_component);
          state_data_n[i_ipt*num_state_data + i_component] = elem_data_n[elem*num_element_data + state_data_offset];
        }
      }

      // DJL todo properly handle state data
      material_->GetStress(elem_global_ids[elem],
                           num_int_pt_per_elem,
                           time_previous,
                           time_current,
                           def_grad_n,
                           def_grad_np1,
                           cauchy_stress_n,
                           cauchy_stress_np1,
                           state_data_n,
                           state_data_np1,
                           data_manager,
                           is_output_step);

      // Copy data to the global containers
      for (int i_ipt = 0 ; i_ipt < num_int_pt_per_elem ; i_ipt++) {
        for (int i_component = 0 ; i_component < full_tensor_size ; i_component++) {
          int def_grad_offset = def_grad_offset_.at(i_ipt*full_tensor_size + i_component);
          elem_data_np1[elem*num_element_data + def_grad_offset] = def_grad_np1[i_ipt*full_tensor_size + i_component];
        }
        for (int i_component = 0 ; i_component < sym_tensor_size ; i_component++) {
          int stress_offset = stress_offset_.at(i_ipt*sym_tensor_size + i_component);
          elem_data_np1[elem*num_element_data + stress_offset] = cauchy_stress_np1[i_ipt*sym_tensor_size + i_component];
        }
        for (int i_component = 0 ; i_component < num_state_data ; i_component++) {
          int state_data_offset = state_data_offset_.at(i_ipt*num_state_data + i_component);
          elem_data_np1[elem*num_element_data + state_data_offset] = state_data_np1[i_ipt*num_state_data + i_component];
        }
      }

      if (!compute_stress_only) {

        // Compute element contribution to nodal force
        element_->ComputeNodalForces(cur_coord,
                                     cauchy_stress_np1,
                                     force);

        // Copy internal force to global containers
        for (int node=0 ; node<num_node_per_elem ; node++) {
          int node_id = elem_conn[elem*num_node_per_elem + node];
          for (int i=0 ; i<vector_size ; i++) {
            internal_force[vector_size*node_id + i] += force[node*vector_size + i];
          }
        }
      }

#ifdef NIMBLE_HAVE_UQ
      if(uq_enabled) {
        //---START MAIN LOOP OVER NUMBER OF UQ OFFNOMINAL SAMPLES--//
        for(int sample=0; sample < num_uq_samples; sample++) {

          //--STEP 1: Apply Taylor's series closure for deformation gradient of THIS sample
          uq_params->ApplyTaylorsSeriesClosure(sample,
                                               full_tensor_size*num_int_pt_per_elem,
                                               def_grad_np1,
                                               def_grad_sensitivities,
                                               def_grad_this_sample );


          //--STEP 2: Apply Material Model to Compute Stress of THIS sample
          material_->GetOffNominalStress(uq_params->GetParamsForSample(sample),
                                         bulk_modulus_uq_index_,
                                         shear_modulus_uq_index_,
                                         num_int_pt_per_elem,
                                         &(def_grad_this_sample[0]),
                                         cauchy_stress_this_sample);

          //--STEP 3: Perform assembly of stresses into nodal forces for THIS sample
          double * displacement_this_sample = offnominal_displacements[sample];
          for (int node=0 ; node<num_node_per_elem ; node++) {
            int node_id = elem_conn[elem*num_node_per_elem + node];
            for (int i=0 ; i<vector_size ; i++) {
              cur_coord[node*vector_size + i] = reference_coordinates[vector_size*node_id + i] + displacement_this_sample[vector_size*node_id + i];
            }
          }

          element_->ComputeNodalForces(cur_coord,
                                       cauchy_stress_this_sample,
                                       force);

          double * internal_force_this_sample = offnominal_internal_forces[sample];
          for (int node=0 ; node<num_node_per_elem ; node++) {
            int node_id = elem_conn[elem*num_node_per_elem + node];
            for (int i=0 ; i<vector_size ; i++) {
              internal_force_this_sample[vector_size*node_id + i] += force[node*vector_size + i];
            }
          }

        }
      }
#endif
    }
  }

  template <typename MatT>
  void Block::ComputeTangentStiffnessMatrix(int num_global_unknowns,
                                            const double * const reference_coordinates,
                                            const double * const displacement,
                                            int num_elem,
                                            const int * const elem_conn,
                                            const int * const global_node_ids,
                                            MatT & tangent_stiffness) const {

    int dim = element_->Dim();
    int num_node_per_elem = element_->NumNodesPerElement();
    int num_int_pt_per_elem = element_->NumIntegrationPointsPerElement();

    int vector_size = LengthToInt(VECTOR, dim);
    int full_tensor_size = LengthToInt(FULL_TENSOR, dim);
    int sym_tensor_size = LengthToInt(SYMMETRIC_TENSOR, dim);

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

      material_->GetTangent(num_int_pt_per_elem, material_tangent);

      element_->ComputeTangent(cur_coord, material_tangent, element_tangent);

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
  void Block::ComputeTangentStiffnessMatrix<MatrixContainer>(int num_global_unknowns,
                                                             const double * const reference_coordinates,
                                                             const double * const displacement,
                                                             int num_elem,
                                                             const int * const elem_conn,
                                                             const int * const global_node_ids,
                                                             MatrixContainer & tangent_stiffness) const ;
  template
  void Block::ComputeTangentStiffnessMatrix<CRSMatrixContainer>(int num_global_unknowns,
                                                                const double * const reference_coordinates,
                                                                const double * const displacement,
                                                                int num_elem,
                                                                const int * const elem_conn,
                                                                const int * const global_node_ids,
                                                                CRSMatrixContainer & tangent_stiffness) const ;

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

    int dim = element_->Dim();
    int vector_size = LengthToInt(VECTOR, dim);
    int num_node_per_elem = element_->NumNodesPerElement();
    int num_int_pt_per_elem = element_->NumIntegrationPointsPerElement();
    int num_data_per_int_pt = num_elem_data / num_int_pt_per_elem;

    // Containers that will be passed to an individual element to compute element volume and volume-averaged quantities
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
          cur_coord[node*vector_size + i] = reference_coordinates[vector_size*node_id + i] + displacement[vector_size*node_id + i];
        }
      }
      for (int i = 0 ; i < num_vol_ave_data*num_int_pt_per_elem ; ++i) {
        int_pt_quantities[i] = elem_data_np1[i_elem*num_data_per_int_pt*num_int_pt_per_elem + vol_ave_offsets_.at(i)];
      }

      // Compute element volume and volume-averaged quantities
      element_->ComputeVolumeAverage(cur_coord,
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
    int dim = element_->Dim();
    int num_int_pt_per_elem = element_->NumIntegrationPointsPerElement();
    int num_data_per_element = static_cast<int>( elem_data_labels.size() );
    int num_data_per_int_pt = num_data_per_element / num_int_pt_per_elem;
    int full_tensor_size = LengthToInt(FULL_TENSOR, dim);
    int sym_tensor_size = LengthToInt(SYMMETRIC_TENSOR, dim);

    // Determine the offset for the derformation gradient in the element data container
    def_grad_offset_.clear();
    for (int i_ipt = 0 ; i_ipt < num_int_pt_per_elem ; i_ipt++) {
      for (int i_component = 0 ; i_component < full_tensor_size ; i_component++) {
        std::string def_grad_label = GetComponentLabel("deformation_gradient", FULL_TENSOR, dim, i_component, i_ipt+1);
        for (unsigned int i_label = 0 ; i_label < num_data_per_element ; i_label++) {
          if (elem_data_labels.at(i_label) == def_grad_label) {
            def_grad_offset_.push_back(i_label);
          }
        }
      }
    }
    if (def_grad_offset_.size() != full_tensor_size*num_int_pt_per_elem) {
      throw std::logic_error("\nError in Block::ComputeInternalForce(), failed to index def_grad into global data.\n");
    }

    // Determine the offset for the stress in the element data container
    stress_offset_.clear();
    for (int i_ipt = 0 ; i_ipt < num_int_pt_per_elem ; i_ipt++) {
      for (int i_component = 0 ; i_component < sym_tensor_size ; i_component++) {
        std::string stress_label = GetComponentLabel("stress", SYMMETRIC_TENSOR, dim, i_component, i_ipt+1);
        for (unsigned int i_label = 0 ; i_label < num_data_per_element ; i_label++) {
          if (elem_data_labels.at(i_label) == stress_label) {
            stress_offset_.push_back(i_label);
          }
        }
      }
    }
    if (stress_offset_.size() != sym_tensor_size*num_int_pt_per_elem) {
      throw std::logic_error("\nError in Block::ComputeInternalForce(), failed to index stress into global data.\n");
    }

    // Determine the offset for the state data in the element data container
    int num_state_data = material_->NumStateVariables();
    std::vector<std::string> state_data_labels(num_state_data);
    for (int i=0 ; i<num_state_data ; ++i) {
      char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN];
      material_->GetStateVariableLabel(i, label);
      state_data_labels[i] = label;
    }
    state_data_offset_.clear();
    for (int i_ipt = 0 ; i_ipt < num_int_pt_per_elem ; i_ipt++) {
      for (int i_state_data = 0 ; i_state_data < num_state_data ; i_state_data++) {
        std::string raw_label = state_data_labels.at(i_state_data);
        std::string state_data_label = GetComponentLabel(raw_label, SCALAR, dim, 0, i_ipt+1);
        for (unsigned int i_label = 0 ; i_label < num_data_per_element ; i_label++) {
          if (elem_data_labels.at(i_label) == state_data_label) {
            state_data_offset_.push_back(i_label);
          }
        }
      }
    }
    if (state_data_offset_.size() != num_state_data*num_int_pt_per_elem) {
      throw std::logic_error("\nError in Block::ComputeInternalForce(), failed to index state data into global data.\n");
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
        std::string const & elem_data_label = RemoveIntegrationPointPrefix(elem_data_labels[i_elem_data_label]);
        if (vol_ave_label == elem_data_label) {
          for (int i_int_pt = 0 ; i_int_pt < num_int_pt_per_elem ; ++i_int_pt) {
            vol_ave_offsets_[i_vol_ave_label + i_int_pt*num_vol_ave_data] = i_elem_data_label + i_int_pt*num_data_per_int_pt;
          }
        }
      }
    }
  }

} // namespace nimble
