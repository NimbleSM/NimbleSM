/*
// @HEADER
//  ************************************************************************
//
//                                 NimbleSM
//                              Copyright 2018
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
//  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
//  retains certain rights in this software.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//  1. Redistributions of source code must retain the above copyright
//  notice, this list of conditions and the following disclaimer.
//
//  2. Redistributions in binary form must reproduce the above copyright
//  notice, this list of conditions and the following disclaimer in the
//  documentation and/or other materials provided with the distribution.
//
//  3. Neither the name of the Corporation nor the names of the
//  contributors may be used to endorse or promote products derived from
//  this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
//  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
//  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
//  NO EVENT SHALL NTESS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Questions?  Contact David Littlewood (djlittl@sandia.gov)
//
//  ************************************************************************
// @HEADER
*/

#include "quasistatic_time_integrator.h"

#include <iostream>

#include "../nimble.h"
#include "../nimble_boundary_condition_manager.h"
#include "../nimble_data_manager.h"
#include "../nimble_mesh_utils.h"
#include "../nimble_parser.h"

namespace nimble {
namespace {

double
ComputeQuasistaticResidual(
    nimble::GenesisMesh&              mesh,
    nimble::DataManager&              data_manager,
    nimble::BoundaryConditionManager& bc,
    int                               linear_system_num_unknowns,
    std::vector<int>&                 linear_system_global_node_ids,
    double                            time_previous,
    double                            time_current,
    const nimble::Viewify<2>&         displacement,
    nimble::Viewify<2>&               internal_force,
    double*                           residual_vector,
    bool                              is_output_step)
{
  const int dim          = mesh.GetDim();
  const int num_nodes    = static_cast<int>(mesh.GetNumNodes());
  const int num_unknowns = num_nodes * mesh.GetDim();

  auto model_data = data_manager.GetModelData();

  model_data->ComputeInternalForce(
      data_manager, time_previous, time_current, is_output_step, displacement, internal_force);

  for (int i = 0; i < linear_system_num_unknowns; i++) residual_vector[i] = 0.0;

  for (int n = 0; n < num_nodes; n++) {
    int ls_id = linear_system_global_node_ids[n];
    for (int dof = 0; dof < dim; dof++) {
      int ls_index = ls_id * dim + dof;
#ifdef NIMBLE_DEBUG
      if (ls_index < 0 || ls_index > linear_system_num_unknowns - 1) {
        throw std::invalid_argument(
            "\nError:  Invalid index into residual vector in "
            "QuasistaticTimeIntegrator().\n");
      }
#endif
      residual_vector[ls_index] += -1.0 * internal_force(n, dof);
    }
  }
  bc.ModifyRHSForKinematicBC(linear_system_global_node_ids.data(), residual_vector);

  double l2_norm = InnerProduct(num_nodes, residual_vector, residual_vector);
  l2_norm        = sqrt(l2_norm);

  double infinity_norm(0.0);
  for (int i = 0; i < num_nodes; i++) { infinity_norm = std::max(infinity_norm, std::abs(residual_vector[i])); }
#ifdef NIMBLE_HAVE_MPI
  double restmp = infinity_norm;
  MPI_Allreduce(&restmp, &infinity_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  double residual = l2_norm + 20.0 * infinity_norm;
  return residual;
}
}  // namespace

QuasistaticTimeIntegrator::QuasistaticTimeIntegrator(
    NimbleApplication& app,
    GenesisMesh&       mesh,
    DataManager&       data_manager)
    : IntegratorBase(app, mesh, data_manager)
{
}

int
QuasistaticTimeIntegrator::Integrate()
{
  auto& parser       = App().GetParser();
  auto& mesh         = Mesh();
  auto& data_manager = GetDataManager();

  const int my_rank   = parser.GetRankID();
  const int num_ranks = parser.GetNumRanks();

  if (num_ranks > 1) {
    std::string msg =
        " Quasi-statics currently not implemented in parallel "
        "(work in progress).\n";
    throw std::invalid_argument(msg);
  }

  int status = 0;

  int        dim             = mesh.GetDim();
  int        num_nodes       = static_cast<int>(mesh.GetNumNodes());
  int        num_blocks      = static_cast<int>(mesh.GetNumBlocks());
  const int* global_node_ids = mesh.GetNodeGlobalIds();

  auto& bc = *(data_manager.GetBoundaryConditionManager());

  // Store various mappings from global to local node ids
  std::vector<int>                linear_system_global_node_ids;
  std::map<int, std::vector<int>> map_from_linear_system;
  for (int n = 0; n < num_nodes; ++n) {
    int global_node_id = global_node_ids[n];
    linear_system_global_node_ids.push_back(global_node_id);
    map_from_linear_system[global_node_id] = std::vector<int>();
    map_from_linear_system[global_node_id].push_back(n);
  }
  int linear_system_num_nodes    = static_cast<int>(map_from_linear_system.size());
  int linear_system_num_unknowns = linear_system_num_nodes * dim;

  std::vector<double> global_data;

  auto* model_data_ptr = dynamic_cast<nimble::ModelData*>(data_manager.GetModelData().get());
  if (model_data_ptr == nullptr) { throw std::invalid_argument(" Incompatible Model Data \n"); }
  nimble::ModelData& model_data = *model_data_ptr;

  // Set up the global vectors
  unsigned int num_unknowns = num_nodes * mesh.GetDim();

  auto reference_coordinate     = model_data.GetVectorNodeData("reference_coordinate");
  auto physical_displacement    = model_data.GetVectorNodeData("displacement");
  auto displacement_fluctuation = model_data.GetVectorNodeData("displacement_fluctuation");
  auto trial_displacement       = model_data.GetVectorNodeData("trial_displacement");

  auto velocity = model_data.GetVectorNodeData("velocity");

  auto internal_force       = model_data.GetVectorNodeData("internal_force");
  auto trial_internal_force = model_data.GetVectorNodeData("trial_internal_force");

  std::vector<double> residual_vector(linear_system_num_unknowns, 0.0);
  std::vector<double> trial_residual_vector(linear_system_num_unknowns, 0.0);
  std::vector<double> linear_solver_solution(linear_system_num_unknowns, 0.0);

  nimble::Viewify<2> displacement = physical_displacement;

  nimble::CRSMatrixContainer tangent_stiffness;
  nimble::CGScratchSpace     cg_scratch;
  std::vector<int>           i_index, j_index;
  nimble::DetermineTangentMatrixNonzeroStructure(mesh, linear_system_global_node_ids, i_index, j_index);
  tangent_stiffness.AllocateNonzeros(i_index, j_index);
  if (my_rank == 0) {
    std::cout << "Number of nonzeros in tangent stiffness matrix = " << tangent_stiffness.NumNonzeros() << "\n"
              << std::endl;
  }

#ifdef NIMBLE_HAVE_TRILINOS
  //  tpetra_container.AllocateTangentStiffnessMatrix(mesh);
  //  if (my_rank == 0) {
  //    std::cout << "Number of nonzeros in the crs tangent stiffness matrix = "
  //    << tpetra_container.TangentStiffnessMatrixNumNonzeros() << "\n" <<
  //    std::endl;
  //  }
#endif

  double initial_time = parser.InitialTime();
  double final_time   = parser.FinalTime();
  double time_current{initial_time};
  double time_previous{initial_time};
  double delta_time{0.0};
  int    num_load_steps           = parser.NumLoadSteps();
  int    output_frequency         = parser.OutputFrequency();
  double user_specified_time_step = (final_time - initial_time) / num_load_steps;

  if (my_rank == 0 && final_time < initial_time) {
    std::string msg = "Final time: " + std::to_string(final_time) +
                      " is less than initial time: " + std::to_string(initial_time) + "\n";
    throw std::invalid_argument(msg);
  }

#ifdef NIMBLE_HAVE_UQ
  if (parser.UseUQ()) NIMBLE_ABORT("\nError:  UQ enabled but not implemented for quasistatics.\n");
#endif

  model_data.ApplyKinematicConditions(data_manager, time_current, time_previous);

  std::vector<double> identity(dim * dim, 0.0);
  for (int i = 0; i < dim; i++) { identity[i] = 1.0; }

  auto& blocks = model_data.GetBlocks();

  data_manager.WriteOutput(time_current);

  if (my_rank == 0) { std::cout << "Beginning quasistatic time integration:" << std::endl; }

  for (int step = 0; step < num_load_steps; step++) {
    time_previous = time_current;
    time_current += user_specified_time_step;
    delta_time = time_current - time_previous;

    bool is_output_step = false;
    if (output_frequency != 0) {
      if (step % output_frequency == 0 || step == num_load_steps - 1) { is_output_step = true; }
    }

    model_data.ApplyKinematicConditions(data_manager, time_current, time_previous);

    // Compute the residual, which is a norm of the (rearranged) nodal force
    // vector, with the dof associated with kinematic BC removed.
    double residual = ComputeQuasistaticResidual(
        mesh,
        data_manager,
        bc,
        linear_system_num_unknowns,
        linear_system_global_node_ids,
        time_previous,
        time_current,
        displacement,
        internal_force,
        residual_vector.data(),
        is_output_step);

    int    iteration(0);
    int    max_nonlinear_iterations = parser.NonlinearSolverMaxIterations();
    double convergence_tolerance    = parser.NonlinearSolverRelativeTolerance() * residual;

    if (my_rank == 0) {
      std::cout << "\nStep " << step + 1 << std::scientific << std::setprecision(3) << ", time " << time_current
                << ", delta_time " << delta_time << ", convergence tolerance " << convergence_tolerance << std::endl;
      std::cout << "  iteration " << iteration << ": residual = " << residual << std::endl;
    }

    while (residual > convergence_tolerance && iteration < max_nonlinear_iterations) {
      tangent_stiffness.SetAllValues(0.0);
#ifdef NIMBLE_HAVE_TRILINOS
      //      tpetra_container.TangentStiffnessMatrixSetScalar(0.0);
#endif
      for (auto& block_it : blocks) {
        int        block_id          = block_it.first;
        int        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
        int const* elem_conn         = mesh.GetConnectivity(block_id);
        auto&      block             = block_it.second;
        block->ComputeTangentStiffnessMatrix(
            linear_system_num_unknowns,
            reference_coordinate.data(),
            displacement.data(),
            num_elem_in_block,
            elem_conn,
            linear_system_global_node_ids.data(),
            tangent_stiffness);
      }

      double diagonal_entry(0.0);
      for (int i = 0; i < linear_system_num_unknowns; ++i) { diagonal_entry += std::abs(tangent_stiffness(i, i)); }
      diagonal_entry /= linear_system_num_unknowns;

      // For the dof with kinematic BC, zero out the rows and columns and put a
      // non-zero on the diagonal
      bc.ModifyTangentStiffnessMatrixForKinematicBC(
          linear_system_num_unknowns, linear_system_global_node_ids.data(), diagonal_entry, tangent_stiffness);
      bc.ModifyRHSForKinematicBC(linear_system_global_node_ids.data(), residual_vector.data());

      // Solve the linear system with the tangent stiffness matrix
      int num_cg_iterations(0);
      std::fill(linear_solver_solution.begin(), linear_solver_solution.end(), 0.0);
      bool success = nimble::CG_SolveSystem(
          tangent_stiffness, residual_vector.data(), cg_scratch, linear_solver_solution.data(), num_cg_iterations);
      if (!success) {
        if (my_rank == 0) { std::cout << "\n**** CG solver failed to converge!\n" << std::endl; }
        status = 1;
        return status;
      }

      //
      // Apply a line search
      //

      // evaluate residual for alpha = 1.0
      for (int n = 0; n < linear_system_num_nodes; n++) {
        std::vector<int> const& node_ids = map_from_linear_system[n];
        for (auto const& node_id : node_ids) {
          for (int dof = 0; dof < dim; dof++) {
            int ls_index = n * dim + dof;
#ifdef NIMBLE_DEBUG
            if (ls_index < 0 || ls_index > linear_solver_solution.size() - 1) {
              throw std::invalid_argument(
                  "\nError:  Invalid index into linear solver solution in "
                  "QuasistaticTimeIntegrator().\n");
            }
#endif
            trial_displacement(node_id, dof) = displacement(node_id, dof) - linear_solver_solution[ls_index];
          }
        }
      }
      double trial_residual = ComputeQuasistaticResidual(
          mesh,
          data_manager,
          bc,
          linear_system_num_unknowns,
          linear_system_global_node_ids,
          time_previous,
          time_current,
          trial_displacement,
          trial_internal_force,
          trial_residual_vector.data(),
          is_output_step);

      //
      // secant line search
      //
      double sr        = nimble::InnerProduct(linear_solver_solution, residual_vector);
      double s_trial_r = nimble::InnerProduct(linear_solver_solution, trial_residual_vector);
      double alpha     = -1.0 * sr / (s_trial_r - sr);

      // evaluate residual for alpha computed with secant line search
      for (int n = 0; n < linear_system_num_nodes; n++) {
        std::vector<int> const& node_ids = map_from_linear_system[n];
        for (auto const& node_id : node_ids) {
          for (int dof = 0; dof < dim; dof++) {
            int ls_index = n * dim + dof;
#ifdef NIMBLE_DEBUG
            if (ls_index < 0 || ls_index > linear_solver_solution.size() - 1) {
              throw std::invalid_argument(
                  "\nError:  Invalid index into linear solver solution vector in "
                  "QuasistaticTimeIntegrator().\n");
            }
#endif
            displacement(node_id, dof) -= alpha * linear_solver_solution[ls_index];
          }
        }
      }
      residual = ComputeQuasistaticResidual(
          mesh,
          data_manager,
          bc,
          linear_system_num_unknowns,
          linear_system_global_node_ids,
          time_previous,
          time_current,
          displacement,
          internal_force,
          residual_vector.data(),
          is_output_step);

      // if the alpha = 1.0 result was better, use alpha = 1.0
      if (trial_residual < residual) {
        residual = trial_residual;
        displacement.copy(trial_displacement);
        internal_force.copy(trial_internal_force);
        for (int i = 0; i < linear_system_num_unknowns; i++) { residual_vector[i] = trial_residual_vector[i]; }
      }

      iteration += 1;

      if (my_rank == 0) {
        std::cout << "  iteration " << iteration << ": residual = " << residual
                  << ", linear cg iterations = " << num_cg_iterations << std::endl;
      }

    }  // while (residual > convergence_tolerance ... )

    if (iteration == max_nonlinear_iterations) {
      if (my_rank == 0) {
        std::cout << "\n**** Nonlinear solver failed to converge!\n" << std::endl;
        std::cout << "**** Relevant input deck parameters for the nonlinear solver:" << std::endl;
        std::cout << "****   nonlinear solver relative tolerance:  " << parser.NonlinearSolverRelativeTolerance()
                  << std::endl;
        std::cout << "****   nonlinear solver maximum iterations:  " << parser.NonlinearSolverMaxIterations()
                  << std::endl;
      }
      status = 1;
      return status;
    } else {
      if (is_output_step) data_manager.WriteOutput(time_current);
    }

    // swap states
    model_data.UpdateStates(data_manager);
  }

  if (my_rank == 0) { std::cout << "\nComplete.\n" << std::endl; }

  return status;
}
}  // namespace nimble