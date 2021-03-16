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

#include "nimble_version.h"
#include "nimble_mpi.h"

#include "nimble_parser.h"
#include "nimble_genesis_mesh.h"
#include "nimble_boundary_condition_manager.h"
#include "nimble_data_manager.h"
#include "nimble_linear_solver.h"
#include "nimble_utils.h"
#include "nimble_mesh_utils.h"
#include "nimble_contact_interface.h"
#include "nimble_contact_manager.h"
#include "nimble_view.h"
#include "nimble_block_material_interface_factory_base.h"
#include "nimble_material_factory.h"
#include "nimble_model_data.h"
#include "nimble_model_data_utils.h"
#include "nimble_vector_communicator.h"
#include "nimble.quanta.stopwatch.h"
#include "nimble_timing_utils.h"

#ifdef NIMBLE_HAVE_BVH
#include "contact/parallel/bvh_contact_manager.h"
#endif

#ifdef NIMBLE_HAVE_MPI
#include <mpi.h>
#endif

#ifdef NIMBLE_HAVE_TRILINOS
#include "nimble_tpetra_utils.h"
#endif

#ifdef NIMBLE_HAVE_VT
#include <vt/transport.h>
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <random>
#include <limits>
#include <fstream>
#include <sstream>

#include "nimble_timer.h"

namespace nimble {

namespace details {

int ExplicitTimeIntegrator(
    const nimble::Parser &parser,
    nimble::GenesisMesh &mesh,
    nimble::DataManager &data_manager,
    std::shared_ptr<nimble::ContactInterface> contact_interface
);

int QuasistaticTimeIntegrator(
    const nimble::Parser &parser,
    nimble::GenesisMesh &mesh,
    nimble::DataManager &data_manager
);

double ComputeQuasistaticResidual(nimble::GenesisMesh & mesh,
                                  nimble::DataManager & data_manager,
                                  nimble::BoundaryConditionManager & bc,
                                  int linear_system_num_unknowns,
                                  std::vector<int> & linear_system_global_node_ids,
                                  double time_previous, double time_current,
                                  const nimble::Viewify<2> &displacement,
                                  nimble::Viewify<2> &internal_force,
                                  double* residual_vector,
                                  double const * rve_macroscale_deformation_gradient,
                                  bool is_output_step);

int parseCommandLine(int argc, char **argv, nimble::Parser &parser)
{
  for (int ia = 1; ia < argc; ++ia) {
    std::string my_arg = std::string(argv[ia]);
    if (my_arg == "--use_vt") {
#ifdef NIMBLE_HAVE_VT
      parser.SetToUseVT();
#else
      std::cerr << "\n Flag '--use_vt' ignored \n\n";
#endif
      continue;
    }
    if (my_arg == "--use_kokkos") {
#ifdef NIMBLE_HAVE_KOKKOS
      parser.SetToUseKokkos();
#else
      std::cerr << "\n Flag '--use_kokkos' ignored \n\n";
#endif
      continue;
    }
    if (my_arg == "--use_tpetra") {
#ifdef NIMBLE_HAVE_TRILINOS
      parser.SetToUseTpetra();
#else
      std::cerr << "\n Flag '--use_tpetra' ignored \n\n";
#endif
      continue;
    }
    //
    if (my_arg == "--use_uq") {
#ifdef NIMBLE_HAVE_UQ
      parser.SetToUseUQ();
#else
      std::cerr << "\n Flag '--use_uq' ignored \n\n";
#endif
      continue;
    }
    //
    parser.SetInputFilename(std::string(my_arg));
  }
  return 0;
}


}


void NimbleInitializeAndGetInput(int argc, char **argv, nimble::Parser &parser) {

  int my_rank = 0, num_ranks = 1;

  // --- Parse the command line
  details::parseCommandLine(argc, argv, parser);

#ifdef NIMBLE_HAVE_TRILINOS
  if (parser.UseTpetra()) {
    auto sguard = new Tpetra::ScopeGuard(&argc,&argv);
    parser.ResetTpetraScope(sguard);
    auto comm = Tpetra::getDefaultComm();
    num_ranks = comm->getSize();
    my_rank = comm->getRank();
  }
  else
#endif
#ifdef NIMBLE_HAVE_MPI
  {
    MPI_Init(&argc, &argv);
    int mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if (mpi_err != MPI_SUCCESS) {
      throw std::logic_error(
          "\nError:  MPI_Comm_rank() returned nonzero error code.\n");
    }
    mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    if (mpi_err != MPI_SUCCESS) {
      throw std::logic_error(
          "\nError:  MPI_Comm_size() returned nonzero error code.\n");
    }
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
    throw std::runtime_error("\nError: Inappropriate set of parameters.\n");
  }

  // Banner
  if (my_rank == 0) {
    std::cout << "\n-- NimbleSM" << std::endl;
    std::cout << "-- version " << nimble::NimbleVersion() << "\n";
    if (parser.UseKokkos()) {
      std::cout << "-- Using Kokkos interface \n";
    }
    else if (parser.UseTpetra()) {
      std::cout << "-- Using Tpetra interface \n";
    }
    else if (parser.UseVT()) {
      std::cout << "-- Using VT runtime \n";
    }
    std::cout << "-- Number of rank";
    if (num_ranks > 1)
      std::cout << "(s)";
    std::cout << " = " << num_ranks << "\n";
    std::cout << std::endl;
  }

  // Initialize VT if needed
#ifdef NIMBLE_HAVE_VT
  if (parser.UseVT() == true) {
    MPI_Comm vt_comm = MPI_COMM_WORLD;
    ::vt::CollectiveOps::initialize(argc, argv, ::vt::no_workers, true, &vt_comm );
  }
#endif

  parser.SetRankID(my_rank);
  parser.SetNumRanks(num_ranks);

  parser.Initialize();

}

int NimbleMain(const std::shared_ptr<MaterialFactoryType> &material_factory_base,
               std::shared_ptr<nimble::ContactInterface> contact_interface,
               const std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase>& block_material,
               const nimble::Parser &parser)
{

  const int my_rank = parser.GetRankID();
  const int num_ranks = parser.GetNumRanks();

#ifdef NIMBLE_HAVE_BVH
  auto comm_mpi = MPI_COMM_WORLD;
  //bvh::vt::context vt_ctx{ argc, argv, &comm_mpi };
#endif

  // Read the mesh
  nimble::GenesisMesh mesh;
  nimble::GenesisMesh rve_mesh;
  {
    std::string genesis_file_name = nimble::IOFileName(parser.GenesisFileName(), "g", "", my_rank, num_ranks);
    std::string rve_genesis_file_name = nimble::IOFileName(parser.RVEGenesisFileName(), "g");
    mesh.ReadFile(genesis_file_name);
    if (rve_genesis_file_name != "none") {
      rve_mesh.ReadFile(rve_genesis_file_name);
    }
  }

  std::string tag = parser.GetOutputTag();
  std::string output_exodus_name = nimble::IOFileName(parser.ExodusFileName(), "e", tag, my_rank, num_ranks);

  int dim = mesh.GetDim();
  int num_nodes = static_cast<int>(mesh.GetNumNodes());

  nimble::DataManager data_manager(parser, mesh, rve_mesh);
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

  auto macroscale_data = data_manager.GetMacroScaleData();
  macroscale_data->InitializeBlocks(data_manager, material_factory_base);

  //
  // Initialize the output file
  //

  data_manager.InitializeOutput(output_exodus_name);

  int status = 0;
  const auto time_integration_scheme = parser.TimeIntegrationScheme();
  if (time_integration_scheme == "explicit") {
    status = details::ExplicitTimeIntegrator(parser, mesh, data_manager,
                                             contact_interface);
  }
  else if (time_integration_scheme == "quasistatic") {
    status = details::QuasistaticTimeIntegrator(parser, mesh, data_manager);
  }

  return status;
}


void NimbleFinalize(const nimble::Parser &parser) {

#ifdef NIMBLE_HAVE_VT
  while ( !::vt::curRT->isTerminated() )
      ::vt::runScheduler();
#endif

#ifdef NIMBLE_HAVE_KOKKOS
  Kokkos::finalize();
#endif

#ifdef NIMBLE_HAVE_TRILINOS
  if (!parser.UseTpetra()) {
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

int ExplicitTimeIntegrator(
    const nimble::Parser &parser,
    nimble::GenesisMesh &mesh,
    nimble::DataManager &data_manager,
    std::shared_ptr<nimble::ContactInterface> contact_interface
)
{

  int dim = mesh.GetDim();
  int num_nodes = static_cast<int>(mesh.GetNumNodes());
  int num_blocks = static_cast<int>(mesh.GetNumBlocks());
  const int my_rank = parser.GetRankID();
  const int num_ranks = parser.GetNumRanks();
  const long num_unknowns = num_nodes * mesh.GetDim();

#ifdef NIMBLE_HAVE_BVH
  nimble::BvhContactManager contact_manager(contact_interface, parser.ContactDicing());
#else
  nimble::ContactManager contact_manager(contact_interface);
#endif

  bool contact_enabled = parser.HasContact();
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

  auto myVectorCommunicator = data_manager.GetVectorCommunicator();

  if (contact_enabled) {
    std::vector<std::string> contact_primary_block_names, contact_secondary_block_names;
    double penalty_parameter;
    nimble::ParseContactCommand(parser.ContactString(),
                                contact_primary_block_names,
                                contact_secondary_block_names,
                                penalty_parameter);
    std::vector<int> contact_primary_block_ids, contact_secondary_block_ids;
    mesh.BlockNamesToOnProcessorBlockIds(contact_primary_block_names,
                                         contact_primary_block_ids);
    mesh.BlockNamesToOnProcessorBlockIds(contact_secondary_block_names,
                                         contact_secondary_block_ids);
    contact_manager.SetPenaltyParameter(penalty_parameter);
    contact_manager.CreateContactEntities(mesh, *myVectorCommunicator,
                                          contact_primary_block_ids,
                                          contact_secondary_block_ids);
    if (contact_visualization) {
      std::string tag = parser.GetOutputTag();
      std::string contact_visualization_exodus_file_name = nimble::IOFileName(parser.ContactVisualizationFileName(), "e", tag, my_rank, num_ranks);
      contact_manager.InitializeContactVisualization(contact_visualization_exodus_file_name);
    }
  }

  nimble::ModelData &model_data = to_ModelData(data_manager.GetMacroScaleData());

  int status = 0;

  //
  // Extract view for global field vectors
  //
  auto reference_coordinate = model_data.GetVectorNodeData("reference_coordinate");
  auto displacement = model_data.GetVectorNodeData("displacement");
  auto velocity = model_data.GetVectorNodeData("velocity");
  auto acceleration = model_data.GetVectorNodeData("acceleration");

  //
  // "View" objects for storing the internal and external forces
  // These forces will be filled and updated inside ModelData member routines.
  //
  auto internal_force = model_data.GetVectorNodeData("internal_force");
  auto external_force = model_data.GetVectorNodeData("external_force");

  nimble::Viewify<2> contact_force;
  if (contact_enabled)
    contact_force = model_data.GetVectorNodeData("contact_force");

  std::map<int, nimble::Block>& blocks = model_data.GetBlocks();
  std::map<int, nimble::Block>::iterator block_it;

  model_data.ComputeLumpedMass(data_manager);

  double critical_time_step = model_data.GetCriticalTimeStep();
  const auto lumped_mass = model_data.GetScalarNodeData("lumped_mass");

  double time_current(0.0), time_previous(0.0);
  double final_time = parser.FinalTime();
  double delta_time, half_delta_time;
  int num_load_steps = parser.NumLoadSteps();
  int output_frequency = parser.OutputFrequency();

  model_data.ApplyInitialConditions(data_manager);

  model_data.ApplyKinematicConditions(data_manager, 0.0, 0.0);

  // For explicit dynamics, the macroscale model is never treated as an RVE
  // so rve_macroscale_deformation_gradient will always be the identity matrix
  std::vector<double> rve_macroscale_deformation_gradient(dim*dim, 0.0);
  for (int i=0 ; i<dim ; i++) {
    rve_macroscale_deformation_gradient[i] = 1.0;
  }

  data_manager.WriteOutput(time_current);

  if (contact_visualization) {
    contact_manager.ContactVisualizationWriteStep(time_current);
  }

  double user_specified_time_step = final_time/num_load_steps;

  if (my_rank == 0) {
    std::cout << "\nUser specified time step:              " << user_specified_time_step << std::endl;
    std::cout << "Approximate maximum stable time step:  " << critical_time_step << "\n" << std::endl;
    if (user_specified_time_step > critical_time_step) {
      std::cout << "**** WARNING:  The user specified time step exceeds the computed maximum stable time step.\n" << std::endl;
    }
    std::cout << "Explicit time integration:\n    0% complete" << std::endl;
  }

  // Timing occurs for this portion of the code
  nimble::quanta::stopwatch main_simulation_timer;
  main_simulation_timer.reset();

  nimble::TimeKeeper total_step_time;
  nimble::TimeKeeper total_contact_time;
  nimble::TimeKeeper total_dynamics_time;
  double total_exodus_write_time = 0.0, total_vector_reduction_time = 0.0;

  for (int step=0 ; step<num_load_steps ; step++) {
    total_step_time.Start();

    if (my_rank == 0) {
      if (10*(step+1) % num_load_steps == 0 && step != num_load_steps - 1) {
        std::cout << "   " << static_cast<int>( 100.0 * static_cast<double>(step+1)/num_load_steps ) << "% complete" << std::endl << std::flush;
      }
      else if (step == num_load_steps - 1) {
        std::cout << "  100% complete\n" << std::endl << std::flush;
      }
    }
    bool is_output_step = false;
    if (output_frequency != 0) {
      if (step%output_frequency == 0 || step == num_load_steps - 1) {
        is_output_step = true;
      }
    }

    total_dynamics_time.Start();

    time_previous = time_current;
    time_current += final_time/num_load_steps;
    delta_time = time_current - time_previous;
    half_delta_time = 0.5*delta_time;

    // V^{n+1/2} = V^{n} + (dt/2) * A^{n}
    velocity += half_delta_time * acceleration;

    model_data.UpdateWithNewVelocity(data_manager, half_delta_time);

    model_data.ApplyKinematicConditions(data_manager, time_current, time_previous);

    // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
    displacement += delta_time * velocity;

    model_data.UpdateWithNewDisplacement(data_manager, delta_time);

    //
    // Evaluate external body forces
    //
    model_data.ComputeExternalForce(data_manager, time_previous, time_current,
                                    is_output_step);

    //
    // Evaluate the internal force
    //
    model_data.ComputeInternalForce(data_manager, time_previous, time_current,
                                    is_output_step, displacement, internal_force);

    total_dynamics_time.Stop();

    // Evaluate the contact force
    if (contact_enabled) {
      total_contact_time.Start();
      contact_manager.ApplyDisplacements(displacement.data());
#ifdef NIMBLE_HAVE_BVH
      contact_manager.ComputeParallelContactForce(step + 1, is_output_step);
#else
      contact_manager.ComputeContactForce(step+1, contact_visualization && is_output_step);
#endif
      contact_manager.GetForces(contact_force.data());
      const int vector_dimension = 3;
      myVectorCommunicator->VectorReduction(vector_dimension, contact_force.data());
      total_contact_time.Stop();
    }

    total_dynamics_time.Start();

    // fill acceleration vector A^{n+1} = M^{-1} ( F^{n} + b^{n} )
    if (contact_enabled) {
      for (int i = 0; i < num_nodes; ++i) {
        const double oneOverM = 1.0 / lumped_mass(i);
        acceleration(i, 0) = oneOverM * (internal_force(i, 0)
                                         + external_force(i, 0)
                                         + contact_force(i, 0));
        acceleration(i, 1) = oneOverM * (internal_force(i, 1)
                                         + external_force(i, 1)
                                           + contact_force(i, 1));
        acceleration(i, 2) = oneOverM * (internal_force(i, 2)
                                         + external_force(i, 2)
                                           + contact_force(i, 2));
      }
    }
    else {
      for (int i = 0; i < num_nodes; ++i) {
        const double oneOverM = 1.0 / lumped_mass(i);
        acceleration(i, 0) = oneOverM * (internal_force(i, 0)
                                         + external_force(i, 0));
        acceleration(i, 1) = oneOverM * (internal_force(i, 1)
                                         + external_force(i, 1));
        acceleration(i, 2) = oneOverM * (internal_force(i, 2)
                                         + external_force(i, 2));
      }
    }

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    velocity += half_delta_time * acceleration;

    model_data.UpdateWithNewVelocity(data_manager, half_delta_time);

    total_dynamics_time.Stop();

    if (is_output_step) {
      nimble::quanta::stopwatch exodus_write_timer;
      exodus_write_timer.reset();

      model_data.ApplyKinematicConditions(data_manager, time_current, time_previous);

      data_manager.WriteOutput(time_current);

      if (contact_visualization) {
        contact_manager.ContactVisualizationWriteStep(time_current);
      }
      total_exodus_write_time += exodus_write_timer.age();

    } // if (is_output_step)

    model_data.UpdateStates(data_manager);
    total_step_time.Stop();
  }

  double total_simulation_time = main_simulation_timer.age();
  if (my_rank == 0) {
    std::cout << "======== Timing data: ========\n";
    std::cout << "Total step time: " << total_step_time.GetElapsedTime() << '\n';
    std::cout << "Total contact time: " << total_contact_time.GetElapsedTime() << '\n';
    std::cout << "Total dynamics time: " << total_dynamics_time.GetElapsedTime() << '\n';
    std::cout << "Total search time: " << contact_manager.total_search_time.GetElapsedTime() << '\n';
    std::cout << "Total enforcement time: " << contact_manager.total_enforcement_time.GetElapsedTime() << '\n';
    std::cout << "Total num contacts: " << contact_manager.total_num_contacts << '\n';
  }
  if (my_rank == 0 && parser.WriteTimingDataFile()) {
    nimble::TimingInfo timing_writer{
        num_ranks,
        nimble::quanta::stopwatch::get_microsecond_timestamp(),
        total_simulation_time,
        0.0, 0.0,
        total_exodus_write_time,
        total_vector_reduction_time
    };
    timing_writer.BinaryWrite();
  }
  return status;
}


int QuasistaticTimeIntegrator(const nimble::Parser &parser,
                              nimble::GenesisMesh & mesh,
                              nimble::DataManager & data_manager
)
{

  const int my_rank = parser.GetRankID();
  const int num_ranks = parser.GetNumRanks();

  if (num_ranks > 1) {
    std::cerr << "Error:  Quasi-statics currently not implemented (work in progress).\n" << std::endl;
    throw std::runtime_error("No Quasistatics in parallel");
  }

  int status = 0;

  int dim = mesh.GetDim();
  int num_nodes = static_cast<int>(mesh.GetNumNodes());
  int num_blocks = static_cast<int>(mesh.GetNumBlocks());
  const int *global_node_ids = mesh.GetNodeGlobalIds();

  auto &bc = *( data_manager.GetBoundaryConditionManager() );

  // Store various mappings from global to local node ids
  // Things get a bit tricky for periodic boundary conditions,
  // where dof are condensed out of the linear system
  std::vector<int> linear_system_global_node_ids;
  std::map<int, std::vector<int>> map_from_linear_system;
  if (bc.IsPeriodicRVEProblem()) {
    int rve_corner_node_id(-1);
    mesh.CreatePeriodicRVELinearSystemMap(global_node_ids,
                                          linear_system_global_node_ids,
                                          map_from_linear_system,
                                          rve_corner_node_id);
    bc.CreateRVEFixedCornersBoundaryConditions(rve_corner_node_id);
  }
  else {
    for (int n=0 ; n<num_nodes ; ++n) {
      int global_node_id = global_node_ids[n];
      linear_system_global_node_ids.push_back(global_node_id);
      map_from_linear_system[global_node_id] = std::vector<int>();
      map_from_linear_system[global_node_id].push_back(n);
    }
  }
  int linear_system_num_nodes = static_cast<int>(map_from_linear_system.size());
  int linear_system_num_unknowns = linear_system_num_nodes * dim;

  std::vector<double> global_data;

  nimble::ModelData &model_data = details::to_ModelData(data_manager.GetMacroScaleData());

  // Set up the global vectors
  unsigned int num_unknowns = num_nodes * mesh.GetDim();

  auto reference_coordinate = model_data.GetVectorNodeData("reference_coordinate");
  auto physical_displacement = model_data.GetVectorNodeData("displacement");
  auto displacement_fluctuation = model_data.GetVectorNodeData("displacement_fluctuation");
  auto trial_displacement = model_data.GetVectorNodeData("trial_displacement");

  auto velocity = model_data.GetVectorNodeData("velocity");

  auto internal_force = model_data.GetVectorNodeData("internal_force");
  auto trial_internal_force = model_data.GetVectorNodeData("trial_internal_force");

  std::vector<double> residual_vector(linear_system_num_unknowns, 0.0);
  std::vector<double> trial_residual_vector(linear_system_num_unknowns, 0.0);
  std::vector<double> linear_solver_solution(linear_system_num_unknowns, 0.0);

  nimble::Viewify<2> displacement = physical_displacement;
  if (bc.IsPeriodicRVEProblem()) {
    std::cout << " periodic RVE problem ... displacement is set to displacement_fluctuation \n";
    displacement = displacement_fluctuation;
    model_data.SetUseDisplacementFluctuations();
  }

  nimble::CRSMatrixContainer tangent_stiffness;
  nimble::CGScratchSpace cg_scratch;
  std::vector<int> i_index, j_index;
  nimble::DetermineTangentMatrixNonzeroStructure(mesh, linear_system_global_node_ids, i_index, j_index);
  tangent_stiffness.AllocateNonzeros(i_index, j_index);
  if (my_rank == 0) {
    std::cout << "Number of nonzeros in tangent stiffness matrix = " << tangent_stiffness.NumNonzeros() << "\n" << std::endl;
  }

#ifdef NIMBLE_HAVE_TRILINOS
//  tpetra_container.AllocateTangentStiffnessMatrix(mesh);
//  if (my_rank == 0) {
//    std::cout << "Number of nonzeros in the crs tangent stiffness matrix = " << tpetra_container.TangentStiffnessMatrixNumNonzeros() << "\n" << std::endl;
//  }
#endif

  double time_current(0.0), time_previous(0.0), delta_time(0.0);
  double final_time = parser.FinalTime();
  int num_load_steps = parser.NumLoadSteps();
  int output_frequency = parser.OutputFrequency();

#ifdef NIMBLE_HAVE_UQ
  if (parser.UseUQ())
    throw std::logic_error("\nError:  UQ enabled but not implemented for quasistatics.\n");
#endif

  model_data.ApplyKinematicConditions(data_manager, time_current, time_previous);

  std::vector<double> identity(dim*dim, 0.0);
  for (int i=0 ; i<dim ; i++) {
    identity[i] = 1.0;
  }
  std::vector<double> rve_center;
  if (bc.IsPeriodicRVEProblem()) {
    rve_center = mesh.BoundingBoxCenter();
  }

  std::map<int, nimble::Block>& blocks = model_data.GetBlocks();

  data_manager.WriteOutput(time_current);

  if (my_rank == 0) {
    std::cout << "Beginning quasistatic time integration:" << std::endl;
  }

  auto &rve_macroscale_deformation_gradient = data_manager.GetRVEDeformationGradient();

  for (int step=0 ; step<num_load_steps ; step++) {

    time_previous = time_current;
    time_current += final_time/num_load_steps;
    delta_time = time_current - time_previous;

    bool is_output_step = false;
    if (output_frequency != 0) {
      if (step%output_frequency == 0 || step == num_load_steps - 1) {
        is_output_step = true;
      }
    }

    model_data.ApplyKinematicConditions(data_manager, time_current, time_previous);

    bc.GetRVEMacroscaleDeformationGradient(time_current, rve_macroscale_deformation_gradient.data());

    // Compute the residual, which is a norm of the (rearranged) nodal force vector, with the dof associated
    // with kinematic BC removed.
    double residual = ComputeQuasistaticResidual(mesh,
                                                 data_manager,
                                                 bc,
                                                 linear_system_num_unknowns,
                                                 linear_system_global_node_ids,
                                                 time_previous,
                                                 time_current,
                                                 displacement,
                                                 internal_force,
                                                 residual_vector.data(),
                                                 rve_macroscale_deformation_gradient.data(),
                                                 is_output_step);

    int iteration(0);
    int max_nonlinear_iterations = parser.NonlinearSolverMaxIterations();
    double convergence_tolerance = parser.NonlinearSolverRelativeTolerance() * residual;

    if (my_rank == 0) {
      std::cout << "\nStep " << step+1 << std::scientific << std::setprecision(3) << ", time " << time_current << ", delta_time " << delta_time << ", convergence tolerance " << convergence_tolerance << std::endl;
      std::cout << "  iteration " << iteration << ": residual = " << residual << std::endl;
    }

    while (residual > convergence_tolerance && iteration < max_nonlinear_iterations) {

      tangent_stiffness.SetAllValues(0.0);
#ifdef NIMBLE_HAVE_TRILINOS
//      tpetra_container.TangentStiffnessMatrixSetScalar(0.0);
#endif
      for (auto& block_it : blocks) {
        int block_id = block_it.first;
        int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
        int const * elem_conn = mesh.GetConnectivity(block_id);
        nimble::Block& block = block_it.second;
        block.ComputeTangentStiffnessMatrix(linear_system_num_unknowns,
                                            reference_coordinate.data(),
                                            displacement.data(),
                                            num_elem_in_block,
                                            elem_conn,
                                            linear_system_global_node_ids.data(),
                                            tangent_stiffness);
      }

      double diagonal_entry(0.0);
      for (int i=0 ; i<linear_system_num_unknowns ; ++i) {
        diagonal_entry += std::abs(tangent_stiffness(i,i));
      }
      diagonal_entry /= linear_system_num_unknowns;

      // For the dof with kinematic BC, zero out the rows and columns and put a non-zero on the diagonal
      bc.ModifyTangentStiffnessMatrixForKinematicBC(linear_system_num_unknowns, linear_system_global_node_ids.data(), diagonal_entry, tangent_stiffness);
      bc.ModifyRHSForKinematicBC(linear_system_global_node_ids.data(), residual_vector.data());

      // Solve the linear system with the tangent stiffness matrix
      int num_cg_iterations(0);
      std::fill(linear_solver_solution.begin(), linear_solver_solution.end(), 0.0);
      bool success = nimble::CG_SolveSystem(tangent_stiffness,
                                            residual_vector.data(),
                                            cg_scratch,
                                            linear_solver_solution.data(),
                                            num_cg_iterations);
      if (!success) {
        throw std::logic_error("\nError:  CG linear solver failed to converge.\n");
      }

      //
      // Apply a line search
      //

      // evaluate residual for alpha = 1.0
      for (int n=0 ; n<linear_system_num_nodes ; n++) {
        std::vector<int> const & node_ids = map_from_linear_system[n];
        for (auto const & node_id : node_ids) {
          for (int dof=0 ; dof<dim ; dof++) {
            int ls_index = n * dim + dof;
            if (ls_index < 0 || ls_index > residual_vector.size() - 1) {
              throw std::logic_error("\nError:  Invalid index into residual vector in QuasistaticTimeIntegrator().\n");
            }
            trial_displacement(node_id, dof) = displacement(node_id, dof)
                                               - linear_solver_solution[ls_index];
          }
        }
      }
      double trial_residual = ComputeQuasistaticResidual(mesh,
                                                         data_manager,
                                                         bc,
                                                         linear_system_num_unknowns,
                                                         linear_system_global_node_ids,
                                                         time_previous,
                                                         time_current,
                                                         trial_displacement,
                                                         trial_internal_force,
                                                         trial_residual_vector.data(),
                                                         rve_macroscale_deformation_gradient.data(),
                                                         is_output_step);

      //
      // secant line search
      //
      double sr = nimble::InnerProduct(linear_solver_solution, residual_vector);
      double s_trial_r = nimble::InnerProduct(linear_solver_solution, trial_residual_vector);
      double alpha = -1.0 * sr / (s_trial_r - sr);

      // evaluate residual for alpha computed with secant line search
      for (int n=0 ; n<linear_system_num_nodes ; n++) {
        std::vector<int> const & node_ids = map_from_linear_system[n];
        for (auto const & node_id : node_ids) {
          for (int dof=0 ; dof<dim ; dof++) {
            int ls_index = n * dim + dof;
            if (ls_index < 0 || ls_index > residual_vector.size() - 1) {
              throw std::logic_error("\nError:  Invalid index into residual vector in QuasistaticTimeIntegrator().\n");
            }
            displacement(node_id, dof) -= alpha*linear_solver_solution[ls_index];
          }
        }
      }
      residual = ComputeQuasistaticResidual(mesh,
                                            data_manager,
                                            bc,
                                            linear_system_num_unknowns,
                                            linear_system_global_node_ids,
                                            time_previous,
                                            time_current,
                                            displacement,
                                            internal_force,
                                            residual_vector.data(),
                                            rve_macroscale_deformation_gradient.data(),
                                            is_output_step);

      // if the alpha = 1.0 result was better, use alpha = 1.0
      if (trial_residual < residual) {
        residual = trial_residual;
        displacement.copy(trial_displacement);
        internal_force.copy(trial_internal_force);
        for (int i=0 ; i<linear_system_num_unknowns ; i++) {
          residual_vector[i] = trial_residual_vector[i];
        }
      }

      // if we are solving a standard problem, then displacement is the physical_displacement
      // if we are solving a periodic RVE problem, then displacement is the displacement_fluctuation,
      // and we need to set physical_displacement manually and include the macroscale deformation gradient
      if (bc.IsPeriodicRVEProblem()) {
        for (int n=0 ; n<linear_system_num_nodes ; n++) {
          std::vector<int> const & node_ids = map_from_linear_system[n];
          for (auto const & node_id : node_ids) {
            std::vector<double> const & F = rve_macroscale_deformation_gradient;
            physical_displacement(node_id, 0) = displacement(node_id, 0)
                                         + (F[K_F_XX] - identity[K_F_XX]) * (reference_coordinate(node_id, 0) - rve_center.at(0))
                                         + (F[K_F_XY] - identity[K_F_XY]) * (reference_coordinate(node_id, 1) - rve_center.at(1))
                                         + (F[K_F_XZ] - identity[K_F_XZ]) * (reference_coordinate(node_id, 2) - rve_center.at(2));
            physical_displacement(node_id, 1) = displacement(node_id, 1)
                                         + (F[K_F_YX] - identity[K_F_YX]) * (reference_coordinate(node_id, 0) - rve_center.at(0))
                                         + (F[K_F_YY] - identity[K_F_YY]) * (reference_coordinate(node_id, 1) - rve_center.at(1))
                                         + (F[K_F_YZ] - identity[K_F_YZ]) * (reference_coordinate(node_id, 2) - rve_center.at(2));
            physical_displacement(node_id, 2) = displacement(node_id, 2)
                                         + (F[K_F_ZX] - identity[K_F_ZX]) * (reference_coordinate(node_id, 0) - rve_center.at(0))
                                         + (F[K_F_ZY] - identity[K_F_ZY]) * (reference_coordinate(node_id, 1) - rve_center.at(1))
                                         + (F[K_F_ZZ] - identity[K_F_ZZ]) * (reference_coordinate(node_id, 2) - rve_center.at(2));
          }
        }
      }

      iteration += 1;

      if (my_rank == 0) {
        std::cout << "  iteration " << iteration << ": residual = " << residual << ", linear cg iterations = " << num_cg_iterations << std::endl;
      }

    } // while (residual > convergence_tolerance ... )

    if (iteration == max_nonlinear_iterations) {
      if (my_rank == 0) {
        std::cout << "\n**** Nonlinear solver failed to converge!\n"
                  << std::endl;
        std::cout
            << "**** Relevant input deck parameters for the nonlinear solver:"
            << std::endl;
        std::cout << "****   nonlinear solver relative tolerance:  "
                  << parser.NonlinearSolverRelativeTolerance() << std::endl;
        std::cout << "****   nonlinear solver maximum iterations:  "
                  << parser.NonlinearSolverMaxIterations() << std::endl;
      }
      status = 1;
      return status;
    }
    else {
      if (is_output_step)
        data_manager.WriteOutput(time_current);
    }

    // swap states
    model_data.UpdateStates(data_manager);

  }

  if (my_rank == 0) {
    std::cout << "\nComplete.\n" << std::endl;
  }

  return status;

}


double ComputeQuasistaticResidual(nimble::GenesisMesh & mesh,
                                  nimble::DataManager & data_manager,
                                  nimble::BoundaryConditionManager & bc,
                                  int linear_system_num_unknowns,
                                  std::vector<int> & linear_system_global_node_ids,
                                  double time_previous, double time_current,
                                  const nimble::Viewify<2> &displacement,
                                  nimble::Viewify<2> &internal_force,
                                  double* residual_vector,
                                  double const * rve_macroscale_deformation_gradient,
                                  bool is_output_step) {

  const int dim = mesh.GetDim();
  const int num_nodes = static_cast<int>(mesh.GetNumNodes());
  const int num_unknowns = num_nodes * mesh.GetDim();

  auto model_data = data_manager.GetMacroScaleData();

  model_data->ComputeInternalForce(data_manager, time_previous, time_current,
                                   is_output_step, displacement, internal_force);

  for (int i=0 ; i<linear_system_num_unknowns ; i++)
    residual_vector[i] = 0.0;

  for (int n=0 ; n<num_nodes ; n++) {
    int ls_id = linear_system_global_node_ids[n];
    for (int dof=0 ; dof<dim ; dof++) {
      int ls_index = ls_id * dim + dof;
      if (ls_index < 0 || ls_index > linear_system_num_unknowns - 1) {
        throw std::logic_error("\nError:  Invalid index into residual vector in QuasistaticTimeIntegrator().\n");
      }
      residual_vector[ls_index] += -1.0 * internal_force(n, dof);
    }
  }
  bc.ModifyRHSForKinematicBC(linear_system_global_node_ids.data(), residual_vector);

  double l2_norm = InnerProduct(num_nodes, residual_vector, residual_vector);
  l2_norm = sqrt(l2_norm);

  double infinity_norm(0.0);
  for (int i=0 ; i<num_nodes ; i++) {
    infinity_norm = std::max(infinity_norm, std::abs(residual_vector[i]));
  }
#ifdef NIMBLE_HAVE_MPI
  double restmp = infinity_norm;
  MPI_Allreduce(&restmp, &infinity_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  double residual = l2_norm + 20.0 * infinity_norm;
  return residual;
}

} }
