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

#include "explicit_time_integrator.h"
#include "../nimble.h"
#include "../nimble_contact_manager.h"
#include "../nimble_parser.h"
#include "../nimble_data_manager.h"
#include "../nimble_profiling_timer.h"
#include "../nimble_genesis_mesh.h"
#include "../nimble_timing_utils.h"
#include "../nimble.quanta.stopwatch.h"
#include <iostream>

namespace nimble {
ExplicitTimeIntegrator::ExplicitTimeIntegrator(NimbleApplication& app, GenesisMesh& mesh, DataManager& data_manager)
  : IntegratorBase(app, mesh, data_manager)
{}

int
ExplicitTimeIntegrator::Integrate()
{
  auto &mesh = Mesh();
  auto &parser = App().GetParser();
  auto &data_manager = GetDataManager();
  auto contact_interface = App().GetContactInterface();
  int        dim          = mesh.GetDim();
  int        num_nodes    = static_cast<int>(mesh.GetNumNodes());
  int        num_blocks   = static_cast<int>(mesh.GetNumBlocks());
  const int  my_rank      = App().Rank();
  const int  num_ranks    = App().NumRanks();
  const long num_unknowns = num_nodes * mesh.GetDim();

  auto contact_manager = nimble::GetContactManager(contact_interface, data_manager);

  bool contact_enabled       = parser.HasContact();
  bool contact_visualization = parser.ContactVisualization();

  auto                   myVectorCommunicator = data_manager.GetVectorCommunicator();
  nimble::ProfilingTimer watch_simulation;

  if (contact_enabled) {
    std::vector<std::string> contact_primary_block_names, contact_secondary_block_names;
    double                   penalty_parameter;
    nimble::ParseContactCommand(
        parser.ContactString(), contact_primary_block_names, contact_secondary_block_names, penalty_parameter);
    std::vector<int> contact_primary_block_ids, contact_secondary_block_ids;
    mesh.BlockNamesToOnProcessorBlockIds(contact_primary_block_names, contact_primary_block_ids);
    mesh.BlockNamesToOnProcessorBlockIds(contact_secondary_block_names, contact_secondary_block_ids);
    contact_manager->SetPenaltyParameter(penalty_parameter);
    contact_manager->CreateContactEntities(
        mesh, *myVectorCommunicator, contact_primary_block_ids, contact_secondary_block_ids);
    if (contact_visualization) {
      std::string tag = "out";
      std::string contact_visualization_exodus_file_name =
          nimble::IOFileName(parser.ContactVisualizationFileName(), "e", tag, my_rank, num_ranks);
      contact_manager->InitializeContactVisualization(contact_visualization_exodus_file_name);
    }
  }

  auto& model_data = *(data_manager.GetModelData());

  int status = 0;

  //
  // Extract view for global field vectors
  //
  auto reference_coordinate = model_data.GetVectorNodeData("reference_coordinate");
  auto displacement         = model_data.GetVectorNodeData("displacement");
  auto velocity             = model_data.GetVectorNodeData("velocity");
  auto acceleration         = model_data.GetVectorNodeData("acceleration");

  //
  // "View" objects for storing the internal and external forces
  // These forces will be filled and updated inside ModelData member routines.
  //
  auto internal_force = model_data.GetVectorNodeData("internal_force");
  auto external_force = model_data.GetVectorNodeData("external_force");

  nimble::Viewify<2> contact_force;
  if (contact_enabled) contact_force = model_data.GetVectorNodeData("contact_force");

  model_data.ComputeLumpedMass(data_manager);

  double     critical_time_step = model_data.GetCriticalTimeStep();
  auto const lumped_mass        = model_data.GetScalarNodeData("lumped_mass");

  double initial_time = parser.InitialTime();
  double final_time   = parser.FinalTime();
  double time_current{initial_time};
  double time_previous{initial_time};
  double delta_time{0.0};
  double half_delta_time{0.0};
  int    num_load_steps           = parser.NumLoadSteps();
  int    output_frequency         = parser.OutputFrequency();
  double user_specified_time_step = (final_time - initial_time) / num_load_steps;

  if (my_rank == 0 && final_time < initial_time) {
    std::string msg = "Final time: " + std::to_string(final_time)
        + " is less than initial time: " + std::to_string(initial_time)
        + "\n";
    throw std::invalid_argument(msg);
  }

  watch_simulation.push_region("BC enforcement");
  model_data.ApplyInitialConditions(data_manager);
  model_data.ApplyKinematicConditions(data_manager, 0.0, 0.0);
  watch_simulation.pop_region_and_report_time();

  data_manager.WriteOutput(time_current);

  if (contact_visualization) { contact_manager->ContactVisualizationWriteStep(time_current); }

  if (my_rank == 0) {
    std::cout << "\nUser specified time step:              " << user_specified_time_step << std::endl;
    std::cout << "Approximate maximum stable time step:  " << critical_time_step << "\n" << std::endl;
    if (user_specified_time_step > critical_time_step) {
      std::cout << "**** WARNING:  The user specified time step exceeds the "
                   "computed maximum stable time step.\n"
                   << std::endl;
    }
    std::cout << "Explicit time integration:\n    0% complete" << std::endl;
  }

  //
  // Timing records
  //

  double total_vector_reduction_time = 0.0;
  double total_dynamics_time = 0.0, total_exodus_write_time = 0.0;
  double total_force_time = 0.0, total_contact_time = 0.0;
  watch_simulation.push_region("Time stepping loop");

  nimble::ProfilingTimer     watch_internal;
  std::map<int, std::size_t> contactInfo;

  for (int step = 0; step < num_load_steps; step++) {
    if (my_rank == 0) {
      if (10 * (step + 1) % num_load_steps == 0 && step != num_load_steps - 1) {
        std::cout << "   " << static_cast<int>(100.0 * static_cast<double>(step + 1) / num_load_steps) << "% complete"
        << std::endl
        << std::flush;
      } else if (step == num_load_steps - 1) {
        std::cout << "  100% complete\n" << std::endl << std::flush;
      }
    }
    bool is_output_step = false;
    if (output_frequency != 0) {
      if (step % output_frequency == 0 || step == num_load_steps - 1) { is_output_step = true; }
    }

    time_previous = time_current;
    time_current += user_specified_time_step;
    delta_time      = time_current - time_previous;
    half_delta_time = 0.5 * delta_time;

    watch_internal.push_region("Time Integration Scheme");
    // V^{n+1/2} = V^{n} + (dt/2) * A^{n}
    velocity += half_delta_time * acceleration;

    model_data.UpdateWithNewVelocity(data_manager, half_delta_time);
    total_dynamics_time += watch_internal.pop_region_and_report_time();

    watch_internal.push_region("BC enforcement");
    model_data.ApplyKinematicConditions(data_manager, time_current, time_previous);
    watch_internal.pop_region_and_report_time();

    watch_internal.push_region("Time Integration Scheme");
    // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
    displacement += delta_time * velocity;

    model_data.UpdateWithNewDisplacement(data_manager, delta_time);
    total_dynamics_time += watch_internal.pop_region_and_report_time();

    watch_internal.push_region("BC enforcement");
    model_data.ApplyKinematicConditions(data_manager, time_current, time_previous);
    watch_internal.pop_region_and_report_time();

    //
    // Evaluate external body forces
    //
    watch_internal.push_region("Force calculation");
    model_data.ComputeExternalForce(data_manager, time_previous, time_current, is_output_step);

    //
    // Evaluate the internal force
    //
    model_data.ComputeInternalForce(
        data_manager, time_previous, time_current, is_output_step, displacement, internal_force);
    total_force_time += watch_internal.pop_region_and_report_time();

    // Evaluate the contact force
    if (contact_enabled) {
      watch_internal.push_region("Contact");
      contact_manager->ComputeContactForce(step + 1, contact_visualization && is_output_step, contact_force);
      total_contact_time += watch_internal.pop_region_and_report_time();
      auto tmpNum = contact_manager->numActiveContactFaces();
      if (tmpNum) contactInfo.insert(std::make_pair(step, tmpNum));
    }

    // fill acceleration vector A^{n+1} = M^{-1} ( F^{n} + b^{n} )
    watch_internal.push_region("Time Integration Scheme");
    if (contact_enabled) {
      for (int i = 0; i < num_nodes; ++i) {
        const double oneOverM = 1.0 / lumped_mass(i);
        acceleration(i, 0)    = oneOverM * (internal_force(i, 0) + external_force(i, 0) + contact_force(i, 0));
        acceleration(i, 1)    = oneOverM * (internal_force(i, 1) + external_force(i, 1) + contact_force(i, 1));
        acceleration(i, 2)    = oneOverM * (internal_force(i, 2) + external_force(i, 2) + contact_force(i, 2));
      }
    } else {
      for (int i = 0; i < num_nodes; ++i) {
        const double oneOverM = 1.0 / lumped_mass(i);
        acceleration(i, 0)    = oneOverM * (internal_force(i, 0) + external_force(i, 0));
        acceleration(i, 1)    = oneOverM * (internal_force(i, 1) + external_force(i, 1));
        acceleration(i, 2)    = oneOverM * (internal_force(i, 2) + external_force(i, 2));
      }
    }

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    velocity += half_delta_time * acceleration;

    model_data.UpdateWithNewVelocity(data_manager, half_delta_time);

    total_dynamics_time += watch_internal.pop_region_and_report_time();

    if (is_output_step) {
      watch_internal.push_region("Output");
      //
      model_data.ApplyKinematicConditions(data_manager, time_current, time_previous);
      //
      data_manager.WriteOutput(time_current);
      //
      if (contact_visualization) { contact_manager->ContactVisualizationWriteStep(time_current); }
      total_exodus_write_time += watch_internal.pop_region_and_report_time();
    }  // if (is_output_step)

    model_data.UpdateStates(data_manager);
  }
  double total_simulation_time = watch_simulation.pop_region_and_report_time();

  for (int irank = 0; irank < num_ranks; ++irank) {
#ifdef NIMBLE_HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if ((my_rank == irank) && (!contactInfo.empty())) {
      std::cout << " Rank " << irank << " has " << contactInfo.size() << " contact entries "
      << "(out of " << num_load_steps << " time steps)." << std::endl;
      std::cout.flush();
    }
#ifdef NIMBLE_HAVE_MPI
MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  if (my_rank == 0 && parser.WriteTimingDataFile()) {
    nimble::TimingInfo timing_writer{
      num_ranks,
      nimble::quanta::stopwatch::get_microsecond_timestamp(),
      total_simulation_time,
      total_force_time,
      total_contact_time,
      total_exodus_write_time,
      total_vector_reduction_time};
    timing_writer.BinaryWrite();
  }

  if (my_rank == 0) {
    std::cout << "======== Timing data: ========\n";
    std::cout << "Total step time = " << total_simulation_time << '\n';
    std::cout << " --- Update A, V, U: " << total_dynamics_time << '\n';
    std::cout << " --- Force: " << total_force_time << "\n";
    if ((contact_enabled) && (contact_manager)) {
      std::cout << " --- Contact time: " << total_contact_time << '\n';
      auto list_timers = contact_manager->getTimers();
      for (const auto& st_pair : list_timers)
        std::cout << " --- >>> >>> " << st_pair.first << " = " << st_pair.second << "\n";
    }
    if (num_ranks > 1) std::cout << " --- Vector Reduction = " << total_vector_reduction_time << "\n";
    //
    std::cout << " --- Exodus Write = " << total_exodus_write_time << "\n";
    //
  }

  return status;
}
}