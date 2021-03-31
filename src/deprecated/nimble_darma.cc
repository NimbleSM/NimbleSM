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

#include <cstdlib>
#include <iostream>

#include "nimble_boundary_condition_manager.h"
#include "nimble_darma_utils.h"
#include "nimble_data_manager.h"
#include "nimble_exodus_output.h"
#include "nimble_genesis_mesh.h"
#include "nimble_parser.h"
#include "nimble_version.h"

using namespace darma_runtime;
using namespace darma_runtime::keyword_arguments_for_access_handle_collection;
using namespace darma_runtime::keyword_arguments_for_task_creation;

using namespace nimble;

int
darma_main_task(std::vector<std::string> args)
{
  // question:  can the number of physical ranks be determined?
  // darma-build/repos/frontend/src/tests/frontend_validation/test_resource_count.cc
  int         num_ranks     = 1;
  std::string num_ranks_key = "num-virtual-ranks";
  std::string input_deck_name;
  for (unsigned int i = 0; i < args.size(); i++) {
    if (args[i] == num_ranks_key) {
      assert(i + 1 < args.size());
      num_ranks = std::atoi(args[i + 1].c_str());
    }
    if (i == args.size() - 1) { input_deck_name = args[i]; }
  }

  // Banner
  std::cout << "\n-- NimbleSM" << std::endl << std::flush;
  std::cout << "-- version " << nimble::NimbleVersion() << "\n" << std::endl;
  std::cout << "NimbleSM initialized on " << num_ranks << " virtual rank(s)."
            << std::endl
            << std::flush;

  auto macroscale_mesh_collection = initial_access_collection<GenesisMesh>(
      "macroscale mesh", index_range = Range1D<int>(num_ranks));
  auto rve_mesh_collection = initial_access_collection<GenesisMesh>(
      "rve mesh", index_range = Range1D<int>(num_ranks));
  auto data_manager_collection = initial_access_collection<DataManager>(
      "data manager", index_range = Range1D<int>(num_ranks));
  auto boundary_condition_manager_collection =
      initial_access_collection<BoundaryConditionManager>(
          "boundary condition manager", index_range = Range1D<int>(num_ranks));
  auto exodus_output_collection = initial_access_collection<ExodusOutput>(
      "exodus output", index_range = Range1D<int>(num_ranks));
  auto derived_element_data_collection =
      initial_access_collection<DerivedElementDataType>(
          "derived element data", index_range = Range1D<int>(num_ranks));
  auto boundary_node_global_ids_collection =
      initial_access_collection<std::vector<int>>(
          "boundary node global ids", index_range = Range1D<int>(num_ranks));
  auto boundary_ranks_collection = initial_access_collection<std::vector<int>>(
      "boundary ranks", index_range = Range1D<int>(num_ranks));

  std::cout << "\nReading input deck " << input_deck_name << std::endl
            << std::flush;
  Parser parser;
  parser.Initialize(input_deck_name);
  std::cout << "  complete." << std::endl << std::flush;

  double          final_time     = parser.FinalTime();
  int             num_load_steps = parser.NumLoadSteps();
  int             step(0);
  int             output_frequency  = parser.OutputFrequency();
  bool            is_output_step    = false;
  ProgressBarFlag progress_bar_flag = NONE;
  double          time_current(0.0), time_previous(0.0), delta_time;

  int rank_that_reads  = -1;
  int rank_that_writes = -1;

  //  {
  //  auto my_global_node_ids_collection =
  //  initial_access_collection<std::vector<int>>("my global node ids",
  //  index_range=Range1D<int>(num_ranks));
  //  create_concurrent_work<Initialize>(macroscale_mesh_collection,
  //                                     rve_mesh_collection,
  //                                     boundary_condition_manager_collection,
  //                                     data_manager_collection,
  //                                     boundary_node_global_ids_collection,
  //                                     boundary_ranks_collection,
  //                                     my_global_node_ids_collection,
  //                                     exodus_output_collection,
  //                                     parser,
  //                                     index_range = Range1D<int>(num_ranks));
  // }

  create_concurrent_work<ReadGenesisFiles>(
      macroscale_mesh_collection,
      rve_mesh_collection,
      parser,
      index_range = Range1D<int>(num_ranks));

  create_concurrent_work<InitializeDataManager>(
      macroscale_mesh_collection,
      rve_mesh_collection,
      data_manager_collection,
      parser,
      index_range = Range1D<int>(num_ranks));

  {
    auto my_global_node_ids_collection =
        initial_access_collection<std::vector<int>>(
            "my global node ids", index_range = Range1D<int>(num_ranks));
    create_concurrent_work<IdentifyGloballySharedNodes>(
        macroscale_mesh_collection,
        data_manager_collection,
        boundary_node_global_ids_collection,
        boundary_ranks_collection,
        my_global_node_ids_collection,
        index_range = Range1D<int>(num_ranks));
  }

  create_concurrent_work<InitializeBoundaryConditionManager>(
      macroscale_mesh_collection,
      boundary_condition_manager_collection,
      parser,
      index_range = Range1D<int>(num_ranks));

  create_concurrent_work<InitializeExodusOutput>(
      macroscale_mesh_collection,
      data_manager_collection,
      exodus_output_collection,
      parser,
      index_range = Range1D<int>(num_ranks));

  create_concurrent_work<ComputeLumpedMass>(
      macroscale_mesh_collection,
      data_manager_collection,
      index_range = Range1D<int>(num_ranks));

  {
    auto lumped_mass_reduction_buffer_collection =
        initial_access_collection<std::vector<double>>(
            "lumped mass reduction buffer",
            index_range = Range1D<int>(num_ranks));
    create_concurrent_work<DistributedVectorReduction>(
        data_manager_collection,
        boundary_node_global_ids_collection,
        boundary_ranks_collection,
        lumped_mass_reduction_buffer_collection,
        LUMPED_MASS,
        step,
        index_range = Range1D<int>(num_ranks));
  }

  create_concurrent_work<ComputeCriticalTimeStep>(
      macroscale_mesh_collection,
      data_manager_collection,
      parser,
      index_range = Range1D<int>(num_ranks));

  create_concurrent_work<ApplyInitialConditions>(
      boundary_condition_manager_collection,
      data_manager_collection,
      index_range = Range1D<int>(num_ranks));

  create_concurrent_work<ApplyKinematicBC>(
      boundary_condition_manager_collection,
      data_manager_collection,
      time_current,
      time_previous,
      index_range = Range1D<int>(num_ranks));

  create_concurrent_work<ComputeDerivedElementData>(
      macroscale_mesh_collection,
      data_manager_collection,
      derived_element_data_collection,
      index_range = Range1D<int>(num_ranks));

  create_concurrent_work<ExodusWriteStep>(
      data_manager_collection,
      exodus_output_collection,
      derived_element_data_collection,
      time_current,
      index_range = Range1D<int>(num_ranks));

  auto internal_force_reduction_buffer_collection =
      initial_access_collection<std::vector<double>>(
          "internal force reduction buffer",
          index_range = Range1D<int>(num_ranks));

  for (int step = 0; step < num_load_steps; step++) {
    is_output_step = false;
    if (step % output_frequency == 0 || step == num_load_steps - 1) {
      is_output_step = true;
    }

    progress_bar_flag = NONE;
    if (step == 0) {
      progress_bar_flag = FIRST_STEP;
    } else if (
        10 * (step + 1) % num_load_steps == 0 && step != num_load_steps - 1) {
      progress_bar_flag = PRINT_PROGRESS;
    } else if (step == num_load_steps - 1) {
      progress_bar_flag = LAST_STEP;
    }

    time_previous = time_current;
    time_current += final_time / num_load_steps;

    create_concurrent_work<ExplicitTimeStep>(
        macroscale_mesh_collection,
        boundary_condition_manager_collection,
        data_manager_collection,
        boundary_node_global_ids_collection,
        boundary_ranks_collection,
        internal_force_reduction_buffer_collection,
        step,
        num_load_steps,
        time_current,
        time_previous,
        progress_bar_flag,
        is_output_step,
        index_range = Range1D<int>(num_ranks),
        name        = darma_runtime::experimental::backend_hint::LBIteration);

    if (is_output_step) {
      create_concurrent_work<ComputeDerivedElementData>(
          macroscale_mesh_collection,
          data_manager_collection,
          derived_element_data_collection,
          index_range = Range1D<int>(num_ranks));

      create_concurrent_work<ApplyKinematicBC>(
          boundary_condition_manager_collection,
          data_manager_collection,
          time_current,
          time_previous,
          index_range = Range1D<int>(num_ranks));

      create_concurrent_work<ExodusWriteStep>(
          data_manager_collection,
          exodus_output_collection,
          derived_element_data_collection,
          time_current,
          index_range = Range1D<int>(num_ranks));
    }
  }

  return 0;
}

DARMA_REGISTER_TOP_LEVEL_FUNCTION(darma_main_task);
