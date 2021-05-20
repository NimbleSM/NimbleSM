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

#include "nimble_boundary_condition_manager.h"

#include "nimble_linear_solver.h"

namespace nimble {

void
BoundaryConditionManager::Initialize(
    std::map<int, std::string> const&      node_set_names,
    std::map<int, std::vector<int>> const& node_sets,
    std::map<int, std::string> const&      side_set_names,
    std::map<int, std::vector<int>> const& side_sets,
    std::vector<std::string> const&        bc_strings,
    int                                    dim,
    std::string                            time_integration_scheme)
{
  node_set_names_ = node_set_names;
  node_sets_      = node_sets;
  side_set_names_ = side_set_names;
  side_sets_      = side_sets;
  dim_            = dim;

  if (time_integration_scheme == "explicit") {
    time_integration_scheme_ = EXPLICIT;
  } else if (time_integration_scheme == "quasistatic") {
    time_integration_scheme_ = QUASISTATIC;
  }

  for (int i = 0; i < bc_strings.size(); i++) {
    BoundaryCondition bc;
    bool              is_valid = bc.Initialize(dim_, bc_strings[i], node_set_names, side_set_names);
    if (is_valid) { boundary_conditions_.push_back(bc); }
  }
}

bool
BoundaryConditionManager::IsPeriodicRVEProblem() const
{
  bool is_periodic_rve_problem(false);
  for (unsigned int i_bc = 0; i_bc < boundary_conditions_.size(); i_bc++) {
    if (boundary_conditions_.at(i_bc).bc_type_ == BoundaryCondition::PERIODIC_RVE) { is_periodic_rve_problem = true; }
  }
  return is_periodic_rve_problem;
}

template <typename MatT>
void
BoundaryConditionManager::ModifyTangentStiffnessMatrixForKinematicBC(
    int              num_unknowns,
    const int* const global_node_ids,
    double           diagonal_entry,
    MatT&            tangent_stiffness) const
{
  for (unsigned int i_bc = 0; i_bc < boundary_conditions_.size(); i_bc++) {
    BoundaryCondition const& bc          = boundary_conditions_[i_bc];
    int                      node_set_id = bc.node_set_id_;
    int                      coordinate  = bc.coordinate_;
    std::vector<int> const&  node_set    = node_sets_.at(node_set_id);

    if (bc.bc_type_ == BoundaryCondition::PRESCRIBED_DISPLACEMENT ||
        bc.bc_type_ == BoundaryCondition::PRESCRIBED_VELOCITY ||
        bc.bc_type_ == BoundaryCondition::RVE_FIXED_DISPLACEMENT) {
      for (unsigned int n = 0; n < node_set.size(); n++) {
        int index = global_node_ids[node_set[n]] * dim_ + coordinate;
        tangent_stiffness.SetRowValues(index, 0.0);
        tangent_stiffness.SetColumnValues(index, 0.0);
        tangent_stiffness(index, index) = diagonal_entry;
      }
    }
  }
}

// explicit template instantiation for MatrixContainer and CRSMatrixContainer
// matrix types
template void
BoundaryConditionManager::ModifyTangentStiffnessMatrixForKinematicBC<MatrixContainer>(
    int              num_unknowns,
    const int* const global_node_ids,
    double           diagonal_entry,
    MatrixContainer& tangent_stiffness) const;
template void
BoundaryConditionManager::ModifyTangentStiffnessMatrixForKinematicBC<CRSMatrixContainer>(
    int                 num_unknowns,
    const int* const    global_node_ids,
    double              diagonal_entry,
    CRSMatrixContainer& tangent_stiffness) const;

void
BoundaryConditionManager::ModifyRHSForKinematicBC(const int* const global_node_ids, double* rhs) const
{
  for (unsigned int i_bc = 0; i_bc < boundary_conditions_.size(); i_bc++) {
    BoundaryCondition const& bc          = boundary_conditions_[i_bc];
    int                      node_set_id = bc.node_set_id_;
    int                      coordinate  = bc.coordinate_;
    std::vector<int> const&  node_set    = node_sets_.at(node_set_id);

    if (bc.bc_type_ == BoundaryCondition::PRESCRIBED_DISPLACEMENT ||
        bc.bc_type_ == BoundaryCondition::PRESCRIBED_VELOCITY ||
        bc.bc_type_ == BoundaryCondition::RVE_FIXED_DISPLACEMENT) {
      for (unsigned int n = 0; n < node_set.size(); n++) {
        rhs[global_node_ids[node_set[n]] * dim_ + coordinate] = 0.0;
      }
    }
  }
}

void
BoundaryConditionManager::CreateRVEFixedCornersBoundaryConditions(int corner_node_id)
{
  // Create a node set for the RVE corner nodes
  int max_node_set_id(0);
  for (auto const& entry : node_sets_) {
    int node_set_id = entry.first;
    if (node_set_id > max_node_set_id) { max_node_set_id = node_set_id; }
  }
  int                        new_node_set_id   = max_node_set_id + 1;
  std::string                new_node_set_name = "nodelist_" + std::to_string(new_node_set_id);
  std::map<int, std::string> new_node_set_names;
  std::map<int, std::string> new_side_set_names;
  new_node_set_names[new_node_set_id] = new_node_set_name;
  std::vector<int> new_node_set_node_ids;
  new_node_set_node_ids.push_back(corner_node_id);
  node_set_names_[new_node_set_id] = new_node_set_name;
  node_sets_[new_node_set_id]      = new_node_set_node_ids;

  // Create fixed-displacement BCs
  std::string       bc_x_string = "rve_fixed_displacement " + new_node_set_name + " x 0.0";
  BoundaryCondition bc_x;
  bc_x.Initialize(dim_, bc_x_string, new_node_set_names, new_side_set_names);
  boundary_conditions_.push_back(bc_x);
  std::string       bc_y_string = "rve_fixed_displacement " + new_node_set_name + " y 0.0";
  BoundaryCondition bc_y;
  bc_y.Initialize(dim_, bc_y_string, new_node_set_names, new_side_set_names);
  boundary_conditions_.push_back(bc_y);
  std::string       bc_z_string = "rve_fixed_displacement " + new_node_set_name + " z 0.0";
  BoundaryCondition bc_z;
  bc_z.Initialize(dim_, bc_z_string, new_node_set_names, new_side_set_names);
  boundary_conditions_.push_back(bc_z);
}

void
BoundaryConditionManager::GetRVEMacroscaleDeformationGradient(
    double  time,
    double* deformation_gradient,
    double  x,
    double  y,
    double  z)
{
  for (int i = 0; i < dim_ * dim_; i++) { deformation_gradient[i] = 0.0; }
  for (int i = 0; i < dim_; i++) { deformation_gradient[i] = 1.0; }

  bool has_rve_macroscale_def_grad = false;
  for (unsigned int i_bc = 0; i_bc < boundary_conditions_.size(); i_bc++) {
    if (boundary_conditions_[i_bc].bc_type_ == BoundaryCondition::PERIODIC_RVE) { has_rve_macroscale_def_grad = true; }
  }

  // For RVE problems, there should only be the RVE boundary condition plus the
  // fixed-displacement boundary conditions at the RVE corner node
  if (has_rve_macroscale_def_grad && boundary_conditions_.size() != dim_ + 1) {
    throw std::logic_error(
        "\nError:  The periodic RVE boundary condition may not be used in "
        "combination with any other boundary condition.\n");
  }

  if (has_rve_macroscale_def_grad) {
    boundary_conditions_[0].GetRVEMacroscaleDeformationGradient(time, deformation_gradient, x, y, z);
  }
}

}  // namespace nimble
