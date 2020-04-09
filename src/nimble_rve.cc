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

#include "nimble_rve.h"
#include "nimble_utils.h"
#include "nimble_mesh_utils.h"
#include "nimble_data_manager.h"
#include <sstream>
#include <limits>
#include <iostream>

namespace nimble {

  RVE::RVE(std::map<int, std::string> const & material_parameters_string,
           GenesisMesh const & rve_mesh,
           std::string rve_boundary_condition_strategy)
    : Material(MaterialParameters()), material_parameters_string_(material_parameters_string),
      rve_mesh_(rve_mesh),
      boundary_condition_strategy_(PERIODIC_BC),
      origin_x_(0.), origin_y_(0.), origin_z_(0.),
      num_global_data_(-1)
  {
    if (rve_boundary_condition_strategy == "periodic bc") {
      boundary_condition_strategy_ = PERIODIC_BC;
    }
    else if (rve_boundary_condition_strategy == "impose deformation gradient") {
      boundary_condition_strategy_ = IMPOSE_DEFORMATION_GRADIENT;
    }
    else {
      printf("\n****Error in RVE, unrecognized boundary condition strategy.\n");
      exit(1);
    }
  }

  std::map<std::string, double> RVE::ParseParametersString(std::string const & material_parameters_string) const {
    // The first string is the material name, followed by the material properties (key-value pairs)
    size_t space_pos = material_parameters_string.find(" ");
    std::string material_name = material_parameters_string.substr(0, space_pos);
    std::string material_props = material_parameters_string.substr(space_pos+1, material_parameters_string.size());
    std::map<std::string, double> material_parameters;
    std::stringstream ss(material_props);
    std::string key, value_string;
    double value;
    while (std::getline(ss, key, ' ') && std::getline(ss, value_string, ' ')) {
      value = std::atof(value_string.c_str());
      material_parameters[key] = value;
    }
    return material_parameters;
  }

  double RVE::GetDensity() const {
    // Return the minimum density of any material in the RVE
    // which leads to a conservative estimate of the critical time step
    double min_density = std::numeric_limits<double>::max();
    for (auto const & entry : material_parameters_string_) {
      std::map<std::string, double> material_parameters = ParseParametersString(entry.second);
      double density = material_parameters.at("density");
      if (density < min_density) {
        min_density = density;
      }
    }
    return min_density;
  }

  double RVE::GetBulkModulus() const {
    // Return the maximum bulk modulus of any material in the RVE
    // which leads to a conservative estimate of the critical time step
    double max_bulk_modulus = 0.0;
    for (auto const & entry : material_parameters_string_) {
      std::map<std::string, double> material_parameters = ParseParametersString(entry.second);
      double bulk_modulus = material_parameters.at("bulk_modulus");
      if (bulk_modulus > max_bulk_modulus) {
        max_bulk_modulus = bulk_modulus;
      }
    }
    return max_bulk_modulus;
  }

  void RVE::InitializeRVE(int elem_global_id,
                          int integration_point_id,
                          DataManager& data_manager,
                          bool write_exodus_output,
                          MaterialFactory& factory) {

    RVEData& rve_data = data_manager.AllocateRVEData(elem_global_id, integration_point_id);
    ModelData& model_data = rve_data.model_data_;
    std::map<int, std::vector< std::vector<double> > >& derived_elem_data = rve_data.derived_elem_data_;
    std::vector<double>& residual_vector = rve_data.residual_vector_;
    std::vector<double>& linear_solver_solution = rve_data.linear_solver_solution_;
    CRSMatrixContainer& tangent_stiffness = rve_data.tangent_stiffness_;
    rve_data.write_exodus_output_ = write_exodus_output;
    ExodusOutput& exodus_output = rve_data.exodus_output_;

    int dim = rve_mesh_.GetDim();
    int num_nodes = rve_mesh_.GetNumNodes();

    // Global data
    num_global_data_ = 0;
    global_data_labels_.resize(num_global_data_);
    global_data_.resize(num_global_data_);

    // Node data
    int reference_coordinate_field_id = model_data.AllocateNodeData(VECTOR, "reference_coordinate", rve_mesh_.GetNumNodes());
    int displacement_field_id = model_data.AllocateNodeData(VECTOR, "displacement", rve_mesh_.GetNumNodes());
    int displacement_fluctuation_field_id = model_data.AllocateNodeData(VECTOR, "displacement_fluctuation", rve_mesh_.GetNumNodes());
    int velocity_field_id = model_data.AllocateNodeData(VECTOR, "velocity", rve_mesh_.GetNumNodes());
    int internal_force_field_id = model_data.AllocateNodeData(VECTOR, "internal_force", rve_mesh_.GetNumNodes());

    // Blocks
    int num_blocks = rve_mesh_.GetNumBlocks();
    std::map<int, nimble::Block>& blocks = model_data.GetBlocks();
    std::vector<int> block_ids = rve_mesh_.GetBlockIds();
    for (int i=0 ; i<num_blocks ; i++){
      int block_id = block_ids[i];
      blocks[block_id] = nimble::Block();
      blocks[block_id].Initialize(material_parameters_string_.at(block_id), factory);
      std::vector< std::pair<std::string, nimble::Length> > data_labels_and_lengths;
      blocks[block_id].GetDataLabelsAndLengths(data_labels_and_lengths);
      model_data.DeclareElementData(block_id, data_labels_and_lengths);
    }
    model_data.AllocateElementData(rve_mesh_.GetNumElementsInBlock());
    model_data.SpecifyOutputFields("displacement displacement_fluctuation volume stress");
    std::map<int, std::vector<std::string> > const & elem_data_labels = model_data.GetElementDataLabels();
    std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();

    std::vector<int> rve_output_global_elem_ids;

    // Initialize the element data
    for (auto & block_it : blocks) {
      int block_id = block_it.first;
      int num_elem_in_block = rve_mesh_.GetNumElementsInBlock(block_id);
      std::vector<int> const & elem_global_ids = rve_mesh_.GetElementGlobalIdsInBlock(block_id);
      std::vector<double> & elem_data_n = model_data.GetElementDataOld(block_id);
      std::vector<double> & elem_data_np1 = model_data.GetElementDataNew(block_id);
      nimble::Block& block = block_it.second;
      block.InitializeElementData(num_elem_in_block,
                                  elem_global_ids,
                                  rve_output_global_elem_ids,
                                  elem_data_labels.at(block_id),
                                  derived_elem_data_labels.at(block_id),
                                  elem_data_n,
                                  elem_data_np1,
                                  factory,
                                  data_manager);
    }

     // Set up the global vectors
    unsigned int num_unknowns = num_nodes * rve_mesh_.GetDim();
    double* reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
    double* physical_displacement = model_data.GetNodeData(displacement_field_id);
    double* displacement_fluctuation = model_data.GetNodeData(displacement_fluctuation_field_id);
    double* velocity = model_data.GetNodeData(velocity_field_id);

    const double * const ref_coord_x = rve_mesh_.GetCoordinatesX();
    const double * const ref_coord_y = rve_mesh_.GetCoordinatesY();
    const double * const ref_coord_z = rve_mesh_.GetCoordinatesZ();
    for (int i=0 ; i<num_nodes ; i++) {
      reference_coordinate[3*i]   = ref_coord_x[i];
      reference_coordinate[3*i+1] = ref_coord_y[i];
      reference_coordinate[3*i+2] = ref_coord_z[i];
      physical_displacement[3*i]   = 0.0;
      physical_displacement[3*i+1] = 0.0;
      physical_displacement[3*i+2] = 0.0;
      displacement_fluctuation[3*i]   = 0.0;
      displacement_fluctuation[3*i+1] = 0.0;
      displacement_fluctuation[3*i+2] = 0.0;
      velocity[3*i]   = 0.0;
      velocity[3*i+1] = 0.0;
      velocity[3*i+2] = 0.0;
    }
    std::vector<double> center = rve_mesh_.BoundingBoxCenter();
    origin_x_ = center[0];
    origin_y_ = center[1];
    origin_z_ = center[2];

    std::map<int, std::string> const & node_set_names = rve_mesh_.GetNodeSetNames();
    std::map<int, std::vector<int> > const & node_sets = rve_mesh_.GetNodeSets();
    std::vector<std::string> bc_strings;
    bc_strings.push_back("periodic_rve");
    std::string time_integration_scheme = "quasistatic";
    bc_ = BoundaryConditionManager();
    linear_system_node_ids_ = std::vector<int>();
    map_from_linear_system_ = std::map<int, std::vector<int>>();
    bc_.Initialize(node_set_names, node_sets, bc_strings, dim, time_integration_scheme);

    const int * const global_node_ids = rve_mesh_.GetNodeGlobalIds();
    if (boundary_condition_strategy_ == PERIODIC_BC) {
      int rve_corner_node_id(-1);
      rve_mesh_.CreatePeriodicRVELinearSystemMap(global_node_ids,
                                                 linear_system_node_ids_,
                                                 map_from_linear_system_,
                                                 rve_corner_node_id);
      bc_.CreateRVEFixedCornersBoundaryConditions(rve_corner_node_id);
    }
    else {
      for (int n=0 ; n<num_nodes ; ++n) {
        int global_node_id = global_node_ids[n];
        linear_system_node_ids_.push_back(global_node_id);
        map_from_linear_system_[global_node_id] = std::vector<int>();
        map_from_linear_system_[global_node_id].push_back(n);
      }
    }
    int linear_system_num_nodes = static_cast<int>(map_from_linear_system_.size());
    int linear_system_num_unknowns = linear_system_num_nodes * dim;

    std::vector<int> i_index, j_index;
    nimble::DetermineTangentMatrixNonzeroStructure(rve_mesh_, linear_system_node_ids_, i_index, j_index);
    tangent_stiffness.AllocateNonzeros(i_index, j_index);
    residual_vector.resize(linear_system_num_unknowns, 0.0);
    linear_solver_solution.resize(linear_system_num_unknowns, 0.0);

    double* displacement = physical_displacement;
    if (bc_.IsPeriodicRVEProblem()) {
      displacement = displacement_fluctuation;
    }

    // Create and initialize the derived element data
    for (auto & block_it : blocks) {
      int block_id = block_it.first;
      derived_elem_data[block_id] = std::vector< std::vector<double> >();
      nimble::Block& block = block_it.second;
      int num_elem_in_block = rve_mesh_.GetNumElementsInBlock(block_id);
      int const * elem_conn = rve_mesh_.GetConnectivity(block_id);
      std::vector<double> const & elem_data_np1 = model_data.GetElementDataNew(block_id);
      block.ComputeDerivedElementData(reference_coordinate,
                                      displacement,
                                      num_elem_in_block,
                                      elem_conn,
                                      elem_data_labels.at(block_id).size(),
                                      elem_data_np1,
                                      derived_elem_data_labels.at(block_id).size(),
                                      derived_elem_data.at(block_id));
    }

    std::vector<double> global_data_for_output;
    std::vector< std::vector<double> > node_data_for_output;
    std::map<int, std::vector< std::vector<double> > > elem_data_for_output;
    if (rve_data.write_exodus_output_) {
      std::string filename = "rve_elem_" + std::to_string(elem_global_id) + "_ipt_1.e";
      exodus_output.Initialize(filename, rve_mesh_);
      std::vector<std::string> global_data_labels_for_output;
      std::vector<std::string> const & node_data_labels_for_output = model_data.GetNodeDataLabelsForOutput();
      std::map<int, std::vector<std::string> > const & elem_data_labels_for_output = model_data.GetElementDataLabelsForOutput();
      std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
      // Initialize the output file
      exodus_output.InitializeDatabase(rve_mesh_,
                                       global_data_labels_for_output,
                                       node_data_labels_for_output,
                                       elem_data_labels_for_output,
                                       derived_elem_data_labels);
      // Write output for time 0.0
      model_data.GetNodeDataForOutput(node_data_for_output);
      model_data.GetElementDataForOutput(elem_data_for_output);
      double time_current = 0.0;
      exodus_output.WriteStep(time_current,
                              global_data_for_output,
                              node_data_for_output,
                              elem_data_labels_for_output,
                              elem_data_for_output,
                              derived_elem_data_labels,
                              derived_elem_data);
    }
  }

  void RVE::GetStress(int elem_id,
                      int num_pts,
                      double time_previous,
                      double time_current,
                      const double * const deformation_gradient_n,
                      const double * const deformation_gradient_np1,
                      const double * const stress_n,
                      double* stress_np1,
                      const double * const state_data_n,
                      double* state_data_np1,
                      DataManager& data_manager,
                      bool is_output_step) {

    // Cauchy stress
    double* sig = stress_np1;
    double const * def_grad = deformation_gradient_np1;

    int dim = rve_mesh_.GetDim();
    int num_nodes = rve_mesh_.GetNumNodes();
    unsigned int num_unknowns = num_nodes * dim;

    std::vector<double> identity(dim*dim, 0.0);
    for (int i=0 ; i<dim ; i++) {
      identity[i] = 1.0;
    }

    for (int pt = 0 ; pt < num_pts ; pt++){

      RVEData& rve_data = data_manager.GetRVEData(elem_id, pt+1);
      ModelData& model_data = rve_data.model_data_;
      std::map<int, std::vector< std::vector<double> > >& derived_elem_data = rve_data.derived_elem_data_;
      std::vector<double>& residual_vector = rve_data.residual_vector_;
      std::vector<double>& linear_solver_solution = rve_data.linear_solver_solution_;
      CRSMatrixContainer& tangent_stiffness = rve_data.tangent_stiffness_;
      bool write_rve_exodus_output = rve_data.write_exodus_output_;
      ExodusOutput& rve_exodus_output = rve_data.exodus_output_;

      int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
      int displacement_field_id = model_data.GetFieldId("displacement");
      int displacement_fluctuation_field_id = model_data.GetFieldId("displacement_fluctuation");
      int velocity_field_id = model_data.GetFieldId("velocity");
      int internal_force_field_id = model_data.GetFieldId("internal_force");
      double* coord = model_data.GetNodeData(reference_coordinate_field_id);
      double* physical_displacement = model_data.GetNodeData(displacement_field_id);
      double* displacement_fluctuation = model_data.GetNodeData(displacement_fluctuation_field_id);
      double* velocity = model_data.GetNodeData(velocity_field_id);
      double* internal_force = model_data.GetNodeData(internal_force_field_id);
      int linear_system_num_nodes = static_cast<int>(map_from_linear_system_.size());
      int linear_system_num_unknowns = linear_system_num_nodes * dim;
      double* displacement = physical_displacement;
      if (bc_.IsPeriodicRVEProblem()) {
        displacement = displacement_fluctuation;
      }

#ifdef NIMBLE_HAVE_UQ
      std::vector<double*> offnominal_displacements(0), offnominal_internal_forces(0), displacement_sensitivities(0);
      std::vector<Viewify> bc_offnom_velocity_views(0);
#endif

      std::map<int, std::vector<std::string> > const & elem_data_labels = model_data.GetElementDataLabels();
      std::map<int, std::vector<std::string> > const & elem_data_labels_for_output = model_data.GetElementDataLabelsForOutput();
      std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();

      std::map<int, nimble::Block>& blocks = model_data.GetBlocks();

      std::vector<double> global_data_for_output;
      std::vector< std::vector<double> > node_data_for_output;
      std::map<int, std::vector< std::vector<double> > > elem_data_for_output;

      if (boundary_condition_strategy_ == IMPOSE_DEFORMATION_GRADIENT) {

        // Impose the deformation gradient on the entire RVE (too stiff, but a decent first step toward multiscale)
        int i_x, i_y, i_z;
        for (int i_node = 0 ; i_node < num_nodes ; ++i_node) {
          i_x = i_node*dim;
          i_y = i_node*dim + 1;
          i_z = i_node*dim + 2;
          displacement[i_x] = (def_grad[K_F_XX] - 1.0)*(coord[i_x] - origin_x_) +         def_grad[K_F_XY]*(coord[i_y] - origin_y_) +         def_grad[K_F_XZ]*(coord[i_z] - origin_z_);
          displacement[i_y] =         def_grad[K_F_YX]*(coord[i_x] - origin_x_) + (def_grad[K_F_YY] - 1.0)*(coord[i_y] - origin_y_) +         def_grad[K_F_YZ]*(coord[i_z] - origin_z_);
          displacement[i_z] =         def_grad[K_F_ZX]*(coord[i_x] - origin_x_) +         def_grad[K_F_ZY]*(coord[i_y] - origin_y_) + (def_grad[K_F_ZZ] - 1.0)*(coord[i_z] - origin_z_);
          internal_force[i_x] = 0.0;
          internal_force[i_y] = 0.0;
          internal_force[i_z] = 0.0;
        }

        // Compute the stress
        for (auto & id_block_pair : blocks) {
          int block_id = id_block_pair.first;
          int num_elem_in_block = rve_mesh_.GetNumElementsInBlock(block_id);
          int const * elem_conn = rve_mesh_.GetConnectivity(block_id);
          std::vector<int> const & elem_global_ids = rve_mesh_.GetElementGlobalIdsInBlock(block_id);
          std::vector<double> const & elem_data_n = model_data.GetElementDataOld(block_id);
          std::vector<double> & elem_data_np1 = model_data.GetElementDataNew(block_id);
          Block& block = id_block_pair.second;
          bool compute_stress_only = true;
          block.ComputeInternalForce(coord,
                                     displacement,
                                     velocity,
                                     identity.data(),
                                     internal_force,
                                     time_previous,
                                     time_current,
                                     num_elem_in_block,
                                     elem_conn,
                                     elem_global_ids.data(),
                                     elem_data_labels.at(block_id),
                                     elem_data_n,
                                     elem_data_np1,
                                     data_manager,
                                     compute_stress_only,
#ifdef NIMBLE_HAVE_UQ
                                     false,
                                     nullptr,
                                     0,
                                     offnominal_displacements,
                                     offnominal_internal_forces,
                                     displacement_sensitivities,
#endif
                                     false);
        }
        CheckVectorSanity(num_unknowns, internal_force, "RVE internal force");
      }
      else if (boundary_condition_strategy_ == PERIODIC_BC) {

        double delta_time = time_current - time_previous;
        // double final_time = parser.FinalTime();

        std::vector<double> rve_center = rve_mesh_.BoundingBoxCenter();

        bc_.ApplyKinematicBC(time_current, time_previous, Viewify(coord,3), Viewify(displacement,3), Viewify(velocity,3)
#ifdef NIMBLE_HAVE_UQ
                           , bc_offnom_velocity_views
#endif
                            );

        // Compute the residual, which is a norm of the (rearranged) nodal force vector, with the dof associated
        // with kinematic BC removed.
        for (int i=0 ; i<num_unknowns ; ++i) {
          internal_force[i] = 0.0;
        }
        for (auto & block_it : blocks) {
          int block_id = block_it.first;
          int num_elem_in_block = rve_mesh_.GetNumElementsInBlock(block_id);
          int const * elem_conn = rve_mesh_.GetConnectivity(block_id);
          std::vector<int> const & elem_global_ids = rve_mesh_.GetElementGlobalIdsInBlock(block_id);
          nimble::Block& block = block_it.second;
          std::vector<double> const & elem_data_n = model_data.GetElementDataOld(block_id);
          std::vector<double> & elem_data_np1 = model_data.GetElementDataNew(block_id);
          block.ComputeInternalForce(coord,
                                     displacement,
                                     velocity,
                                     def_grad,
                                     internal_force,
                                     time_previous,
                                     time_current,
                                     num_elem_in_block,
                                     elem_conn,
                                     elem_global_ids.data(),
                                     elem_data_labels.at(block_id),
                                     elem_data_n,
                                     elem_data_np1,
                                     data_manager,
                                     false
#ifdef NIMBLE_HAVE_UQ
                                    ,false
                                    ,nullptr
                                    ,0
                                    ,offnominal_displacements
                                    ,offnominal_internal_forces
                                    ,displacement_sensitivities
#endif
);
        }
        CheckVectorSanity(num_unknowns, internal_force, "RVE internal force, initial residual calculation");

        for (int i=0 ; i<linear_system_num_nodes*dim ; i++) {
          residual_vector[i] = 0.0;
        }

        for (int n=0 ; n<num_nodes ; n++) {
          int ls_id = linear_system_node_ids_[n];
          for (int dof=0 ; dof<dim ; dof++) {
            int local_index = n * dim + dof;
            int ls_index = ls_id * dim + dof;
            if (ls_index < 0 || ls_index > residual_vector.size() - 1) {
              throw std::logic_error("\nError:  Invalid index into residual vector in RVE::GetStress().\n");
            }
            residual_vector[ls_index] += -1.0 * internal_force[local_index];
          }
        }

        bc_.ModifyRHSForKinematicBC(linear_system_node_ids_.data(), residual_vector.data());
        double l2_norm(0.0), infinity_norm(0.0);
        for (int i=0 ; i<linear_system_num_unknowns ; i++) {
          double f = std::fabs(residual_vector[i]);
          l2_norm += f * f;
          if (f > infinity_norm) {
            infinity_norm = f;
          }
        }

        l2_norm = std::sqrt(l2_norm);
        double residual = l2_norm + 20.0 * infinity_norm;
        double absolute_convergence_tolerance = 1.0e-9;
        double relative_convergence_tolerance = 1.0e-6 * residual;
        double convergence_tolerance = absolute_convergence_tolerance > relative_convergence_tolerance ? absolute_convergence_tolerance : relative_convergence_tolerance;
        int iteration(0), max_nonlinear_iterations(50);

        bool verbose = false;

        if (verbose) {
          std::cout << "\n  RVE solve for element " << elem_id << ", integration point " << pt+1 << std::endl;
          std::cout << "    time " << time_current << ", delta_time " << delta_time << ", relative tolerance " << relative_convergence_tolerance << ", absolute tolerance " << absolute_convergence_tolerance << std::endl;
          std::cout << "    macroscale deformation gradient " << def_grad[0] << ", " << def_grad[1] << ", " << def_grad[2] << ", " << def_grad[3] << ", " << def_grad[4] << ", " << def_grad[5] << ", " << def_grad[6] << ", " << def_grad[7] << ", " << def_grad[8] << std::endl;
          std::cout << "      iteration " << iteration << " residual " << residual << std::endl;
        }

        while (residual > convergence_tolerance && iteration < max_nonlinear_iterations) {

          // compute the tangent stiffness matrix
          tangent_stiffness.SetAllValues(0.0);
          for (auto & block_it : blocks) {
            int block_id = block_it.first;
            int num_elem_in_block = rve_mesh_.GetNumElementsInBlock(block_id);
            int const * elem_conn = rve_mesh_.GetConnectivity(block_id);
            nimble::Block& block = block_it.second;
            block.ComputeTangentStiffnessMatrix(linear_system_num_unknowns,
                                                coord,
                                                displacement,
                                                num_elem_in_block,
                                                elem_conn,
                                                linear_system_node_ids_.data(),
                                                tangent_stiffness);
          }

          // For the dof with kinematic BC, zero out the rows and columns and put a non-zero on the diagonal
          double diagonal_entry(0.0);
          for (int i=0 ; i<linear_system_num_unknowns ; ++i) {
            diagonal_entry += std::fabs(tangent_stiffness(i,i));
          }
          diagonal_entry /= linear_system_num_unknowns;
          bc_.ModifyTangentStiffnessMatrixForKinematicBC(linear_system_num_unknowns, linear_system_node_ids_.data(), diagonal_entry, tangent_stiffness);

          bc_.ModifyRHSForKinematicBC(linear_system_node_ids_.data(), residual_vector.data());

          std::fill(linear_solver_solution.begin(), linear_solver_solution.end(), 0.0);
          int num_cg_iterations(0);
          bool success = nimble::CG_SolveSystem(tangent_stiffness,
                                                        residual_vector.data(),
                                                        cg_scratch_,
                                                        linear_solver_solution.data(),
                                                        num_cg_iterations);
          if (!success) {
            throw std::logic_error("\nError:  CG linear solver failed to converge.\n");
          }

          // update the trial solution
          for (int n=0 ; n<linear_system_num_nodes ; n++) {
            std::vector<int> const & node_ids = map_from_linear_system_[n];
            for (auto const & node_id : node_ids) {
              for (int dof=0 ; dof<dim ; dof++) {
                int local_index = node_id * dim + dof;
                int ls_index = n * dim + dof;
                if (ls_index < 0 || ls_index > residual_vector.size() - 1) {
                  throw std::logic_error("\nError:  Invalid index into residual vector in RVE::GetStress().\n");
                }
                displacement[local_index] -= linear_solver_solution[ls_index];
              }
              // if we are solving a standard problem, then displacement is the physical_displacement
              // if we are solving a periodic RVE problem, then displacement is the displacement_fluctuation,
              // and we need to set physical_displacement manually and include the macroscale deformation gradient
              if (bc_.IsPeriodicRVEProblem()) {
                int i_x = node_id * dim;
                int i_y = node_id * dim + 1;
                int i_z = node_id * dim + 2;
                // todo:  this is currently hard-coded to 3D
                physical_displacement[i_x] = displacement[i_x]
                  + (def_grad[K_F_XX] - identity[K_F_XX]) * (coord[i_x] - rve_center.at(0))
                  + (def_grad[K_F_XY] - identity[K_F_XY]) * (coord[i_y] - rve_center.at(1))
                  + (def_grad[K_F_XZ] - identity[K_F_XZ]) * (coord[i_z] - rve_center.at(2));
                physical_displacement[i_y] = displacement[i_y]
                  + (def_grad[K_F_YX] - identity[K_F_YX]) * (coord[i_x] - rve_center.at(0))
                  + (def_grad[K_F_YY] - identity[K_F_YY]) * (coord[i_y] - rve_center.at(1))
                  + (def_grad[K_F_YZ] - identity[K_F_YZ]) * (coord[i_z] - rve_center.at(2));
                physical_displacement[i_z] = displacement[i_z]
                  + (def_grad[K_F_ZX] - identity[K_F_ZX]) * (coord[i_x] - rve_center.at(0))
                  + (def_grad[K_F_ZY] - identity[K_F_ZY]) * (coord[i_y] - rve_center.at(1))
                  + (def_grad[K_F_ZZ] - identity[K_F_ZZ]) * (coord[i_z] - rve_center.at(2));
              }
            }
          }

          // check convergence
          for (int i=0 ; i<num_unknowns ; ++i) {
            internal_force[i] = 0.0;
          }
          for (auto & block_it : blocks) {
            int block_id = block_it.first;
            int num_elem_in_block = rve_mesh_.GetNumElementsInBlock(block_id);
            int const * elem_conn = rve_mesh_.GetConnectivity(block_id);
            std::vector<int> const & elem_global_ids = rve_mesh_.GetElementGlobalIdsInBlock(block_id);
            nimble::Block& block = block_it.second;
            std::vector<double> const & elem_data_n = model_data.GetElementDataOld(block_id);
            std::vector<double> & elem_data_np1 = model_data.GetElementDataNew(block_id);
            block.ComputeInternalForce(coord,
                                       displacement,
                                       velocity,
                                       def_grad,
                                       internal_force,
                                       time_previous,
                                       time_current,
                                       num_elem_in_block,
                                       elem_conn,
                                       elem_global_ids.data(),
                                       elem_data_labels.at(block_id),
                                       elem_data_n,
                                       elem_data_np1,
                                       data_manager,
                                       false
#ifdef NIMBLE_HAVE_UQ
                                      ,false
                                      ,nullptr
                                      ,0
                                      ,offnominal_displacements
                                      ,offnominal_internal_forces
                                      ,displacement_sensitivities
#endif
                                      );
          }
          CheckVectorSanity(num_unknowns, internal_force, "RVE internal force");
          for (int i=0 ; i<linear_system_num_nodes*dim ; i++) {
            residual_vector[i] = 0.0;
          }
          for (int n=0 ; n<num_nodes ; n++) {
            int ls_id = linear_system_node_ids_[n];
            for (int dof=0 ; dof<dim ; dof++) {
              int local_index = n * dim + dof;
              int ls_index = ls_id * dim + dof;
              if (ls_index < 0 || ls_index > residual_vector.size() - 1) {
                throw std::logic_error("\nError:  Invalid index into residual vector in RVE::GetStress().\n");
              }
              residual_vector[ls_index] += -1.0 * internal_force[local_index];
            }
          }
          bc_.ModifyRHSForKinematicBC(linear_system_node_ids_.data(), residual_vector.data());
          double l2_norm(0.0), infinity_norm(0.0);
          for (int i=0 ; i<num_nodes ; i++) {
            double f = std::fabs(residual_vector[i]);
            l2_norm += f * f;
            if (f > infinity_norm) {
              infinity_norm = f;
            }
          }
          l2_norm = std::sqrt(l2_norm);
          residual = l2_norm + 20.0 * infinity_norm;

          iteration += 1;

          if (verbose) {
            std::cout << "      iteration " << iteration << " residual " << residual << std::endl;
          }
        } // while (residual > convergence_tolerance ... )

        if (iteration == max_nonlinear_iterations) {
          std::cout << "\n**** RVE nonlinear solver failed to converge in " << max_nonlinear_iterations << " iterations.\n" << std::endl;
          throw std::logic_error("\nError:  RVE nonlinear solver failed to converge.n");
        }

        if (write_rve_exodus_output && is_output_step && pt == 0) {

          for (auto & block_it : blocks) {
            int block_id = block_it.first;
            nimble::Block& block = block_it.second;
            int num_elem_in_block = rve_mesh_.GetNumElementsInBlock(block_id);
            int const * elem_conn = rve_mesh_.GetConnectivity(block_id);
            std::vector<double> const & elem_data_np1 = model_data.GetElementDataNew(block_id);
            block.ComputeDerivedElementData(coord,
                                            displacement,
                                            num_elem_in_block,
                                            elem_conn,
                                            elem_data_labels.at(block_id).size(),
                                            elem_data_np1,
                                            derived_elem_data_labels.at(block_id).size(),
                                            derived_elem_data.at(block_id));
          }

          // Write output
          model_data.GetNodeDataForOutput(node_data_for_output);
          model_data.GetElementDataForOutput(elem_data_for_output);
          rve_exodus_output.WriteStep(time_current,
                                      global_data_for_output,
                                      node_data_for_output,
                                      elem_data_labels_for_output,
                                      elem_data_for_output,
                                      derived_elem_data_labels,
                                      derived_elem_data);

        }

        // todo: swap states
      }

      // Compute the volume-averaged stress
      for (auto & id_block_pair : blocks) {
        int block_id = id_block_pair.first;
        Block& block = id_block_pair.second;
        int num_elem_in_block = rve_mesh_.GetNumElementsInBlock(block_id);
        int const * elem_conn = rve_mesh_.GetConnectivity(block_id);
        std::vector<double> const & elem_data_np1 = model_data.GetElementDataNew(block_id);
        block.ComputeDerivedElementData(coord,
                                        displacement,
                                        num_elem_in_block,
                                        elem_conn,
                                        elem_data_labels.at(block_id).size(),
                                        elem_data_np1,
                                        derived_elem_data_labels.at(block_id).size(),
                                        derived_elem_data.at(block_id));
      }

      int volume_index(-1);
      int stress_xx_index(-1), stress_yy_index(-1), stress_zz_index(-1), stress_xy_index(-1), stress_yz_index(-1), stress_zx_index(-1);
      double total_volume(0.0);
      double vol_ave_stress_xx(0.0), vol_ave_stress_yy(0.0), vol_ave_stress_zz(0.0), vol_ave_stress_xy(0.0), vol_ave_stress_yz(0.0), vol_ave_stress_zx(0.0);
      for (auto const & id_block_pair : blocks) {
        int block_id = id_block_pair.first;
        std::vector<std::string> const & labels = derived_elem_data_labels.at(block_id);
        for (unsigned int i_label = 0 ; i_label < labels.size() ; i_label++) {
          if (labels[i_label] == "volume") {
            volume_index = i_label;
          }
          else if (labels[i_label] == "stress_xx") {
            stress_xx_index = i_label;
          }
          else if (labels[i_label] == "stress_yy") {
            stress_yy_index = i_label;
          }
          else if (labels[i_label] == "stress_zz") {
            stress_zz_index = i_label;
          }
          else if (labels[i_label] == "stress_xy") {
            stress_xy_index = i_label;
          }
          else if (labels[i_label] == "stress_yz") {
            stress_yz_index = i_label;
          }
          else if (labels[i_label] == "stress_zx") {
            stress_zx_index = i_label;
          }
        }
        if (volume_index == -1 ||
            stress_xx_index == -1 ||
            stress_yy_index == -1 ||
            stress_zz_index == -1 ||
            stress_xy_index == -1 ||
            stress_yz_index == -1 ||
            stress_zx_index == -1) {
          throw std::logic_error("**** Error in RVE::GetStress(), unable to process volume-average element data.\n");
        }
        int num_elem_in_block = rve_mesh_.GetNumElementsInBlock(block_id);
        for (int i_elem = 0 ; i_elem < num_elem_in_block ; i_elem++) {
          double element_volume = derived_elem_data.at(block_id).at(volume_index).at(i_elem);
          vol_ave_stress_xx += derived_elem_data.at(block_id).at(stress_xx_index).at(i_elem) * element_volume;
          vol_ave_stress_yy += derived_elem_data.at(block_id).at(stress_yy_index).at(i_elem) * element_volume;
          vol_ave_stress_zz += derived_elem_data.at(block_id).at(stress_zz_index).at(i_elem) * element_volume;
          vol_ave_stress_xy += derived_elem_data.at(block_id).at(stress_xy_index).at(i_elem) * element_volume;
          vol_ave_stress_yz += derived_elem_data.at(block_id).at(stress_yz_index).at(i_elem) * element_volume;
          vol_ave_stress_zx += derived_elem_data.at(block_id).at(stress_zx_index).at(i_elem) * element_volume;
          total_volume += element_volume;
        }
      }
      vol_ave_stress_xx /= total_volume;
      vol_ave_stress_yy /= total_volume;
      vol_ave_stress_zz /= total_volume;
      vol_ave_stress_xy /= total_volume;
      vol_ave_stress_yz /= total_volume;
      vol_ave_stress_zx /= total_volume;

      sig[K_S_XX] = vol_ave_stress_xx;
      sig[K_S_YY] = vol_ave_stress_yy;
      sig[K_S_ZZ] = vol_ave_stress_zz;
      sig[K_S_XY] = vol_ave_stress_xy;
      sig[K_S_YZ] = vol_ave_stress_yz;
      sig[K_S_ZX] = vol_ave_stress_zx;

      sig += 6;
      def_grad += 9;
    }
  }

  void RVE::GetTangent(int num_pts,
                       double* material_tangent) const {
  }
}
