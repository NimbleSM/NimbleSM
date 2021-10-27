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

#include "nimble_main.h"

#include "nimble.quanta.stopwatch.h"
#include "nimble_block_material_interface_factory_base.h"
#include "nimble_boundary_condition_manager.h"
#include "nimble_contact_interface.h"
#include "nimble_contact_manager.h"
#include "nimble_data_manager.h"
#include "nimble_genesis_mesh.h"
#include "nimble_linear_solver.h"
#include "nimble_material_factory.h"
#include "nimble_mesh_utils.h"
#include "nimble_model_data.h"
#include "nimble_parser.h"
#include "nimble_profiling_timer.h"
#include "nimble_timing_utils.h"
#include "nimble_utils.h"
#include "nimble_vector_communicator.h"
#include "nimble_version.h"
#include "nimble_view.h"

#ifdef NIMBLE_HAVE_MPI
#include <mpi.h>
#endif

#ifdef NIMBLE_HAVE_TRILINOS
#include "nimble_tpetra_utils.h"
#endif

#ifdef NIMBLE_HAVE_VT
#include <vt/transport.h>
#endif

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>

#include "nimble_timer.h"

namespace nimble {

namespace details {

int
ExplicitTimeIntegrator(
    const nimble::Parser&                     parser,
    nimble::GenesisMesh&                      mesh,
    nimble::DataManager&                      data_manager,
    std::shared_ptr<nimble::ContactInterface> contact_interface);

int
QuasistaticTimeIntegrator(const nimble::Parser& parser, nimble::GenesisMesh& mesh, nimble::DataManager& data_manager);

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
    bool                              is_output_step);

int
parseCommandLine(int argc, char** argv, nimble::Parser& parser)
{
  for (int ia = 1; ia < argc; ++ia) {
    std::string my_arg = std::string(argv[ia]);
    //
    // Skip flags for Kokkos, Tpetra, and VT as they have already been processed
    //
    if ((my_arg == "--use_kokkos") || (my_arg == "--use_kokkos") || (my_arg == "--use_vt")) continue;
    //
    if (my_arg == "--use_uq") {
#ifdef NIMBLE_HAVE_UQ
      parser.SetToUseUQ();
#else
      std::cerr << "\n Flag '--use_uq' ignored \n\n";
#endif
      continue;
    }

    if (my_arg.substr(0, 4) == "--vt") continue;
    //
    parser.SetInputFilename(std::string(my_arg));
  }
  return 0;
}

}  // namespace details

#ifdef NIMBLE_HAVE_VT
static ::vt::RuntimePtrType vt_rt = nullptr;
#endif

void
NimbleInitializeAndGetInput(int argc, char** argv, nimble::Parser& parser)
{
  int my_rank = 0, num_ranks = 1;

  // --- Parse the command line
  details::parseCommandLine(argc, argv, parser);

#ifdef NIMBLE_HAVE_TRILINOS
  if (parser.UseTpetra()) {
    auto sguard = new Tpetra::ScopeGuard(&argc, &argv);
    parser.ResetTpetraScope(sguard);
    auto comm = Tpetra::getDefaultComm();
    num_ranks = comm->getSize();
    my_rank   = comm->getRank();
  } else
#endif
#ifdef NIMBLE_HAVE_MPI
  {
    MPI_Init(&argc, &argv);
    int mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (mpi_err != MPI_SUCCESS) { NIMBLE_ABORT("\nError:  MPI_Comm_rank() returned nonzero error code.\n"); }
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    if (mpi_err != MPI_SUCCESS) { NIMBLE_ABORT("\nError:  MPI_Comm_size() returned nonzero error code.\n"); }
  }
#endif

#ifdef NIMBLE_HAVE_KOKKOS
  Kokkos::initialize(argc, argv);
#endif

  if (argc < 2) {
    if (my_rank == 0) {
#ifdef NIMBLE_HAVE_MPI
      std::cout << "Usage:  mpirun -np NP NimbleSM <input_deck.in>\n" << std::endl;
#else
      std::cout << "Usage:  NimbleSM <input_deck.in>\n" << std::endl;
#endif
    }
    throw std::invalid_argument("\nError: Inappropriate set of parameters.\n");
  }

  // Banner
  if (my_rank == 0) {
    std::cout << "\n-- NimbleSM" << std::endl;
    std::cout << "-- version " << nimble::NimbleVersion() << "\n";
    if (parser.UseKokkos()) {
      std::cout << "-- Using Kokkos interface \n";
    } else if (parser.UseTpetra()) {
      std::cout << "-- Using Tpetra interface \n";
    } else if (parser.UseVT()) {
      std::cout << "-- Using VT runtime \n";
    }
    std::cout << "-- Number of rank";
    if (num_ranks > 1) std::cout << "(s)";
    std::cout << " = " << num_ranks << "\n";
    std::cout << std::endl;
  }

  // Initialize VT if needed
#ifdef NIMBLE_HAVE_VT
  if (parser.UseVT() == true) {
    MPI_Comm vt_comm = MPI_COMM_WORLD;
    vt_rt            = ::vt::CollectiveOps::initialize(argc, argv, ::vt::no_workers, true, &vt_comm);
    //
    // Check whether we need the next line
    //
    // bvh::vt::context vt_ctx{argc, argv, &comm_mpi};
  }
#endif

  parser.SetRankID(my_rank);
  parser.SetNumRanks(num_ranks);

  parser.Initialize();
}

int
NimbleMain(
    const std::shared_ptr<nimble::MaterialFactoryBase>&               material_factory_base,
    std::shared_ptr<nimble::ContactInterface>                         contact_interface,
    const std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase>& block_material,
    const nimble::Parser&                                             parser)
{
  const int my_rank   = parser.GetRankID();
  const int num_ranks = parser.GetNumRanks();

  // Read the mesh
  nimble::GenesisMesh mesh;
  {
    std::string genesis_file_name = nimble::IOFileName(parser.GenesisFileName(), "g", "", my_rank, num_ranks);
    mesh.ReadFile(genesis_file_name);
  }

  std::string tag                = "out";
  std::string output_exodus_name = nimble::IOFileName(parser.ExodusFileName(), "e", tag, my_rank, num_ranks);

  int dim       = mesh.GetDim();
  int num_nodes = static_cast<int>(mesh.GetNumNodes());

  nimble::DataManager data_manager(parser, mesh);
  data_manager.SetBlockMaterialInterfaceFactory(block_material);

  if (my_rank == 0) {
    std::cout << "\n";
    if (num_ranks == 1) {
      std::cout << " Number of Nodes = " << num_nodes << "\n";
      std::cout << " Number of Elements = " << mesh.GetNumElements() << "\n";
    }
    std::cout << " Number of Global Blocks = " << mesh.GetNumGlobalBlocks() << "\n";
    std::cout << "\n";
    std::cout << " Number of Ranks         = " << num_ranks << "\n";
#ifdef _OPENMP
    std::cout << " Number of Threads       = " << omp_get_max_threads() << "\n";
#endif
    std::cout << "\n";
  }

  auto model_data = data_manager.GetModelData();
  model_data->InitializeBlocks(data_manager, material_factory_base);

  //
  // Initialize the output file
  //

  data_manager.InitializeOutput(output_exodus_name);

  int        status                  = 0;
  const auto time_integration_scheme = parser.TimeIntegrationScheme();
  if (time_integration_scheme == "explicit") {
    status = details::ExplicitTimeIntegrator(parser, mesh, data_manager, contact_interface);
  } else if (time_integration_scheme == "quasistatic") {
    status = details::QuasistaticTimeIntegrator(parser, mesh, data_manager);
  }

  return status;
}

void
NimbleFinalize(const nimble::EnvironmentFlags& env_flags)
{
#ifdef NIMBLE_HAVE_VT
  if (env_flags.use_vt_) ::vt::finalize(std::move(vt_rt));
#endif

#ifdef NIMBLE_HAVE_KOKKOS
  if (env_flags.use_kokkos_) Kokkos::finalize();
#endif

#ifdef NIMBLE_HAVE_TRILINOS
  if (!env_flags.use_tpetra_) {
#ifdef NIMBLE_HAVE_MPI
    MPI_Finalize();
#endif
  }
#else
#ifdef NIMBLE_HAVE_MPI
  MPI_Finalize();
#endif
#endif
}

namespace details {

int
ExplicitTimeIntegrator(
    const nimble::Parser&                     parser,
    nimble::GenesisMesh&                      mesh,
    nimble::DataManager&                      data_manager,
    std::shared_ptr<nimble::ContactInterface> contact_interface)
{
  int        dim          = mesh.GetDim();
  int        num_nodes    = static_cast<int>(mesh.GetNumNodes());
  int        num_blocks   = static_cast<int>(mesh.GetNumBlocks());
  const int  my_rank      = parser.GetRankID();
  const int  num_ranks    = parser.GetNumRanks();
  const long num_unknowns = num_nodes * mesh.GetDim();

  auto contact_manager = nimble::GetContactManager(contact_interface, data_manager);

  bool contact_enabled       = parser.HasContact();
  bool contact_visualization = parser.ContactVisualization();

#ifdef NIMBLE_HAVE_TRILINOS
  //
  // -- This code does not seem to be active
  // -- The lines are commented out for the moment.
  //
/*
  int contact_block_id_1 = 0;
  int contact_block_id_2 = 0;
  double* contact_data_block_1 = 0;
  double* contact_data_block_2 = 0;
  std::vector<int> contact_results_block_1;
  std::vector<int> contact_results_block_2;
  bool has_contact = parser.HasContact();
  if (has_contact) {
    std::stringstream ss(parser.ContactString());
    std::string block_name;
    ss >> block_name;
    contact_block_id_1 = mesh.GetBlockId(block_name);
    contact_results_block_1.resize(mesh.GetNumElementsInBlock(contact_block_id_1));
    ss >> block_name;
    contact_block_id_2 = mesh.GetBlockId(block_name);
    contact_results_block_2.resize(mesh.GetNumElementsInBlock(contact_block_id_2));
  }
*/
#endif

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

int
QuasistaticTimeIntegrator(const nimble::Parser& parser, nimble::GenesisMesh& mesh, nimble::DataManager& data_manager)
{
  const int my_rank   = parser.GetRankID();
  const int num_ranks = parser.GetNumRanks();

  if (num_ranks > 1) {
    std::string msg = " Quasi-statics currently not implemented in parallel "
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
    std::string msg = "Final time: " + std::to_string(final_time)
                      + " is less than initial time: " + std::to_string(initial_time)
                      + "\n";
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
        if (my_rank == 0) {
          std::cout << "\n**** CG solver failed to converge!\n" << std::endl;
        }
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

}  // namespace details
}  // namespace nimble
