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

#include "nimble_kokkos.h"
#include "nimble.quanta.stopwatch.h"
#include "nimble_boundary_condition_manager.h"
#include "nimble_contact_manager.h"
#include "nimble_data_manager.h"
#include "nimble_exodus_output.h"
#include "nimble_kokkos_block.h"
#include "nimble_kokkos_defs.h"
#include "nimble_kokkos_material_factory.h"
#include "nimble_kokkos_model_data.h"
#include "nimble_kokkos_profiling.h"
#include "nimble_model_data_utils.h"
#include "nimble_parser.h"
#include "nimble_timer.h"
#include "nimble_timing_utils.h"
#include "nimble_utils.h"
#include "nimble_version.h"
#include "nimble_view.h"
#include <cassert>
#include <nimble_contact_interface.h>
#include <nimble_kokkos_block_material_interface.h>
#include <nimble_kokkos_block_material_interface_factory.h>

#ifdef NIMBLE_HAVE_ARBORX
#ifdef NIMBLE_HAVE_MPI
    #include "contact/parallel/arborx_parallel_contact_manager.h"
  #endif
    #include "contact/serial/arborx_serial_contact_manager.h"
#endif


#include <iostream>

namespace nimble {

namespace details_kokkos {

int ExplicitTimeIntegrator(const nimble::Parser & parser,
                           nimble::GenesisMesh & mesh,
                           nimble::DataManager & data_manager,
                           nimble::BoundaryConditionManager & boundary_condition_manager,
                           nimble::ExodusOutput & exodus_output,
                           std::shared_ptr<nimble::ContactInterface> contact_interface,
                           std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase> block_material_interface_factory
);

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
    parser.SetInputFilename(std::string(my_arg));
  }
  return 0;
}

}


void NimbleKokkosInitializeAndGetInput(int argc, char **argv, nimble::Parser &parser)
{

  int my_rank = 0, num_ranks = 1;

  // --- Parse the command line
  details_kokkos::parseCommandLine(argc, argv, parser);

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


int NimbleKokkosFinalize(const nimble::Parser &parser) {

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

  return 0;

}

void NimbleKokkosMain(const std::shared_ptr<MaterialFactoryType>& material_factory,
                      std::shared_ptr<nimble::ContactInterface> contact_interface,
                      const std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase> &block_material,
                      const nimble::Parser &parser)
{

  const int my_rank = parser.GetRankID();
  const int num_ranks = parser.GetNumRanks();

  //--- Define timers
  nimble_kokkos::ProfilingTimer watch_simulation;
  watch_simulation.push_region("Parse and read mesh");

  // Read the mesh
  nimble::GenesisMesh mesh;
  nimble::GenesisMesh rve_mesh;
  {
    //--- This part is independent of Kokkos
    std::string genesis_file_name = nimble::IOFileName(parser.GenesisFileName(), "g", "", my_rank, num_ranks);
    std::string rve_genesis_file_name = nimble::IOFileName(parser.RVEGenesisFileName(), "g");
    mesh.ReadFile(genesis_file_name);
    if (rve_genesis_file_name != "none") {
      rve_mesh.ReadFile(rve_genesis_file_name);
    }
  }

  nimble::DataManager data_manager(parser, mesh, rve_mesh);

  watch_simulation.pop_region_and_report_time();

#ifdef NIMBLE_HAVE_ARBORX
  std::string tag = "arborx";
#else
  std::string tag = "kokkos";
#endif
  std::string output_exodus_name = nimble::IOFileName(parser.ExodusFileName(), "e", tag, my_rank, num_ranks);
  int dim = mesh.GetDim();
  int num_nodes = static_cast<int>(mesh.GetNumNodes());
  int num_blocks = static_cast<int>(mesh.GetNumBlocks());

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
  watch_simulation.push_region("Model data and field allocation");

  auto macroscale_data = data_manager.GetMacroScaleData();
  macroscale_data->InitializeBlocks(data_manager, material_factory);

  //
  // Initialize the initial- and boundary-condition manager
  //
  std::map<int, std::string> const & node_set_names = mesh.GetNodeSetNames();
  std::map<int, std::vector<int> > const & node_sets = mesh.GetNodeSets();
  std::vector<std::string> const & bc_strings = parser.GetBoundaryConditionStrings();
  std::string const & time_integration_scheme = parser.TimeIntegrationScheme();
  nimble::BoundaryConditionManager boundary_condition_manager;
  boundary_condition_manager.Initialize(node_set_names, node_sets, bc_strings, dim, time_integration_scheme);

  //
  // Initialize the output file
  //

  std::vector<std::string> global_data_labels;

  nimble::ExodusOutput exodus_output;
  exodus_output.Initialize(output_exodus_name, mesh);

  auto model_data = data_manager.GetMacroScaleData();
  auto &node_data_labels_for_output = model_data->GetNodeDataLabelsForOutput();
  auto &elem_data_labels_for_output = model_data->GetElementDataLabelsForOutput();
  auto &derived_elem_data_labels = model_data->GetDerivedElementDataLabelsForOutput();

  exodus_output.InitializeDatabase(mesh,
                                   global_data_labels,
                                   node_data_labels_for_output,
                                   elem_data_labels_for_output,
                                   derived_elem_data_labels);

  watch_simulation.pop_region_and_report_time();

  if (time_integration_scheme == "explicit") {
    details_kokkos::ExplicitTimeIntegrator(parser, mesh, data_manager,
                                           boundary_condition_manager,
                                           exodus_output,
                                           contact_interface, block_material);
  }
  else {
    throw std::runtime_error("\n Time Integration Scheme Not Implemented \n");
  }

}


namespace details_kokkos {

int ExplicitTimeIntegrator(const nimble::Parser & parser,
                           nimble::GenesisMesh & mesh,
                           nimble::DataManager & data_manager,
                           nimble::BoundaryConditionManager & boundary_condition_manager,
                           nimble::ExodusOutput & exodus_output,
                           std::shared_ptr<nimble::ContactInterface> contact_interface,
                           std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase> block_material_interface_factory
)
{
  const int my_rank = parser.GetRankID();
  const int num_ranks = parser.GetNumRanks();

  int dim = mesh.GetDim();
  int num_nodes = static_cast<int>(mesh.GetNumNodes());
  int num_blocks = static_cast<int>(mesh.GetNumBlocks());

  /// Temporary Solution while refactoring
  auto &field_ids = data_manager.GetFieldIDs();
  nimble_kokkos::ModelData &model_data = ::nimble_kokkos::details_kokkos::to_ModelData(data_manager.GetMacroScaleData());
  auto &exodus_output_manager = model_data.GetExodusOutputManager();
  auto &blocks = model_data.GetBlocks();
  /////////////////

  // Build up block data for stress computation
  std::vector<nimble::BlockData> block_data;
  for (auto &&block_it : blocks)
  {
    int block_id = block_it.first;
    nimble_kokkos::Block &block = block_it.second;
    nimble::Material *material_d = block.GetDeviceMaterialModel();
    int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    int num_integration_points_per_element = block.GetHostElement()->NumIntegrationPointsPerElement();
    block_data.emplace_back(&block, material_d, block_id, num_elem_in_block, num_integration_points_per_element);
  }

  nimble_kokkos::ProfilingTimer watch_simulation;
  watch_simulation.push_region("Lumped mass gather and compute");

  nimble_kokkos::HostVectorNodeView reference_coordinate_h = model_data.GetHostVectorNodeData(field_ids.reference_coordinates);
  nimble_kokkos::DeviceVectorNodeView reference_coordinate_d = model_data.GetDeviceVectorNodeData(field_ids.reference_coordinates);

  nimble_kokkos::HostVectorNodeView displacement_h = model_data.GetHostVectorNodeData(field_ids.displacement);
  nimble_kokkos::DeviceVectorNodeView displacement_d = model_data.GetDeviceVectorNodeData(field_ids.displacement);
  Kokkos::deep_copy(displacement_h, (double)(0.0));

  nimble_kokkos::HostVectorNodeView velocity_h = model_data.GetHostVectorNodeData(field_ids.velocity);
  nimble_kokkos::DeviceVectorNodeView velocity_d = model_data.GetDeviceVectorNodeData(field_ids.velocity);
  Kokkos::deep_copy(velocity_h, (double)(0.0));

  nimble_kokkos::HostVectorNodeView acceleration_h = model_data.GetHostVectorNodeData(field_ids.acceleration);
  Kokkos::deep_copy(acceleration_h, (double)(0.0));

  nimble_kokkos::HostVectorNodeView internal_force_h = model_data.GetHostVectorNodeData(field_ids.internal_force);
  nimble_kokkos::DeviceVectorNodeView internal_force_d = model_data.GetDeviceVectorNodeData(field_ids.internal_force);
  Kokkos::deep_copy(internal_force_h, (double)(0.0));

  nimble_kokkos::HostVectorNodeView contact_force_h = model_data.GetHostVectorNodeData(field_ids.contact_force);
  nimble_kokkos::DeviceVectorNodeView contact_force_d = model_data.GetDeviceVectorNodeData(field_ids.contact_force);
  Kokkos::deep_copy(internal_force_h, (double)(0.0));

  model_data.ComputeLumpedMass(data_manager);
  Viewify lumped_mass = model_data.GetScalarNodeData("lumped_mass");

  auto gathered_reference_coordinate_d = model_data.GetGatheredRefCoord();
  auto &gathered_displacement_d = model_data.GetGatheredDisp();
  auto &gathered_internal_force_d = model_data.GetGatheredInternalForce();
  auto &gathered_contact_force_d = model_data.GetGatheredContactForce();

  double critical_time_step = model_data.GetCriticalTimeStep();

  watch_simulation.push_region("Contact setup");

#ifdef NIMBLE_HAVE_ARBORX
  #ifdef NIMBLE_HAVE_MPI
    nimble::ArborXParallelContactManager contact_manager(contact_interface);
  #else
    nimble::ArborXSerialContactManager contact_manager(contact_interface);
  #endif // NIMBLE_HAVE_MPI
#else
  nimble::ContactManager contact_manager(contact_interface);
#endif // NIMBLE_HAVE_ARBORX

  auto myVectorCommunicator = data_manager.GetVectorCommunicator();

  bool contact_enabled = parser.HasContact();
  bool contact_visualization = parser.ContactVisualization();
  if (contact_enabled) {
    std::vector<std::string> contact_master_block_names, contact_slave_block_names;
    double penalty_parameter;
    nimble::ParseContactCommand(parser.ContactString(),
                                contact_master_block_names,
                                contact_slave_block_names,
                                penalty_parameter);
    std::vector<int> contact_master_block_ids, contact_slave_block_ids;
    mesh.BlockNamesToOnProcessorBlockIds(contact_master_block_names,
                                         contact_master_block_ids);
    mesh.BlockNamesToOnProcessorBlockIds(contact_slave_block_names,
                                         contact_slave_block_ids);
    contact_manager.SetPenaltyParameter(penalty_parameter);
    contact_manager.CreateContactEntities(mesh,
                                          *myVectorCommunicator,
                                          contact_master_block_ids,
                                          contact_slave_block_ids);
    if (contact_visualization) {
#ifdef NIMBLE_HAVE_ARBORX
      std::string tag = "arborx";
#else
      std::string tag = "kokkos";
#endif
      std::string contact_visualization_exodus_file_name = nimble::IOFileName(parser.ContactVisualizationFileName(), "e", tag, my_rank, num_ranks);
      contact_manager.InitializeContactVisualization(contact_visualization_exodus_file_name);
    }
  }

  watch_simulation.pop_region_and_report_time();

  watch_simulation.push_region("BC enforcement");

  double time_current(0.0), time_previous(0.0);
  double final_time = parser.FinalTime();
  double delta_time(0.0), half_delta_time(0.0);
  const int num_load_steps = parser.NumLoadSteps();
  int output_frequency = parser.OutputFrequency();

  boundary_condition_manager.ApplyInitialConditions(reference_coordinate_h, velocity_h);
  boundary_condition_manager.ApplyKinematicBC(time_current, time_previous,
                                              reference_coordinate_h, displacement_h,
                                              velocity_h);
  Kokkos::deep_copy(displacement_d, displacement_h);
  Kokkos::deep_copy(velocity_d, velocity_h);

  watch_simulation.pop_region_and_report_time();

  // Output to Exodus file
  watch_simulation.push_region("Output");

  exodus_output_manager.ComputeElementData(mesh, model_data, blocks,
                                           gathered_reference_coordinate_d,
                                           gathered_displacement_d);
  std::vector<double> global_data;
  std::vector< std::vector<double> > const & node_data_for_output = exodus_output_manager.GetNodeDataForOutput(model_data);

  std::map<int, std::vector< std::vector<double> > > const & elem_data_for_output = exodus_output_manager.GetElementDataForOutput(model_data);
  auto elem_data_labels_for_output = exodus_output_manager.GetElementDataLabelsForOutput();

  std::map<int, std::vector< std::vector<double> > > derived_elem_data;
  auto const &derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();

  exodus_output.WriteStep(time_current,
                          global_data,
                          node_data_for_output,
                          elem_data_labels_for_output,
                          elem_data_for_output,
                          derived_elem_data_labels,
                          derived_elem_data);
  if (contact_visualization) {
    contact_manager.ContactVisualizationWriteStep(time_current);
  }

  watch_simulation.pop_region_and_report_time();
  //
  double total_internal_force_time = 0.0, total_contact_time = 0.0;
  double total_arborx_time = 0.0,
      total_contact_applyd = 0.0, total_contact_getf = 0.0;
  double total_vector_reduction_time = 0.0;
  double total_update_avu_time = 0.0;
  double total_exodus_write_time = 0.0;
  //
  std::map<int, std::size_t> contactInfo;
  //
  watch_simulation.push_region("Time stepping loop");
  nimble_kokkos::ProfilingTimer watch_internal;

  int block_index = 0;
  std::map<int, nimble_kokkos::Block>::iterator block_it;

  for (int step = 0 ; step < num_load_steps ; ++step) {

    if (my_rank == 0) {
      if (10*(step+1) % num_load_steps == 0 && step != num_load_steps - 1) {
        std::cout << "   " << static_cast<int>( 100.0 * static_cast<double>(step+1)/num_load_steps )
                  << "% complete" << std::endl << std::flush;
      }
      else if (step == num_load_steps - 1) {
        std::cout << "  100% complete\n" << std::endl << std::flush;
      }
    }
    bool is_output_step = (step%output_frequency == 0 || step == num_load_steps - 1);

    watch_internal.push_region("Central difference");
    time_previous = time_current;
    time_current += final_time/num_load_steps;
    delta_time = time_current - time_previous;
    half_delta_time = 0.5*delta_time;

    // V^{n+1/2} = V^{n} + (dt/2) * A^{n}
    for (int i=0 ; i<num_nodes ; ++i) {
      velocity_h(i,0) += half_delta_time * acceleration_h(i,0);
      velocity_h(i,1) += half_delta_time * acceleration_h(i,1);
      velocity_h(i,2) += half_delta_time * acceleration_h(i,2);
    }
    total_update_avu_time += watch_internal.pop_region_and_report_time();

    // Apply kinematic boundary conditions
    watch_internal.push_region("BC enforcement");
    boundary_condition_manager.ApplyKinematicBC(time_current, time_previous, reference_coordinate_h, displacement_h, velocity_h);
    watch_internal.pop_region_and_report_time();

    // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
    watch_internal.push_region("Central difference");
    for (int i=0 ; i<num_nodes ; ++i) {
      displacement_h(i,0) += delta_time * velocity_h(i,0);
      displacement_h(i,1) += delta_time * velocity_h(i,1);
      displacement_h(i,2) += delta_time * velocity_h(i,2);
    }

    // Copy the current displacement and velocity value to device memory
    Kokkos::deep_copy(displacement_d, displacement_h);
    Kokkos::deep_copy(velocity_d, velocity_h);

    total_update_avu_time += watch_internal.pop_region_and_report_time();

    watch_internal.push_region("Force calculation");

    nimble_kokkos::ProfilingTimer watch_internal_details;
    Kokkos::deep_copy(internal_force_d, (double)(0.0));

    // Compute element-level kinematics
    constexpr int mpi_vector_dim = 3;

    watch_internal_details.push_region("Element kinematics");
    for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
      //
      int block_id = block_it->first;
      nimble_kokkos::Block& block = block_it->second;
      nimble::Element* element_d = block.GetDeviceElement();
      int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      int num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);

      nimble_kokkos::DeviceElementConnectivityView elem_conn_d = block.GetDeviceElementConnectivityView();
      nimble_kokkos::DeviceVectorNodeGatheredView gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
      nimble_kokkos::DeviceVectorNodeGatheredView gathered_displacement_block_d = gathered_displacement_d.at(block_index);
      nimble_kokkos::DeviceVectorNodeGatheredView gathered_internal_force_block_d = gathered_internal_force_d.at(block_index);

      model_data.GatherVectorNodeData(field_ids.reference_coordinates, // TODO SHOULD JUST PASS IN VIEW?
                                      num_elem_in_block, // TODO SHOULD BE ABLE TO GET THIS OFF VIEW "EXTENT"
                                      num_nodes_per_elem,
                                      elem_conn_d,
                                      gathered_reference_coordinate_block_d);

      model_data.GatherVectorNodeData(field_ids.displacement,
                                      num_elem_in_block,
                                      num_nodes_per_elem,
                                      elem_conn_d,
                                      gathered_displacement_block_d);

      auto deformation_gradient_step_np1_d = model_data
          .GetDeviceFullTensorIntegrationPointData(block_id, field_ids.deformation_gradient, nimble::STEP_NP1);

      // COMPUTE DEFORMATION GRADIENTS
      Kokkos::parallel_for("Deformation Gradient", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
        nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL(), Kokkos::ALL());
        nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d = Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL(), Kokkos::ALL());
        nimble_kokkos::DeviceFullTensorIntPtSubView element_deformation_gradient_step_np1_d = Kokkos::subview(deformation_gradient_step_np1_d, i_elem, Kokkos::ALL(), Kokkos::ALL());
        element_d->ComputeDeformationGradients(element_reference_coordinate_d,
                                               element_displacement_d,
                                               element_deformation_gradient_step_np1_d);
      });
    }
    watch_internal_details.pop_region_and_report_time();

    watch_internal_details.push_region("Material stress calculation");
    auto block_material_interface = block_material_interface_factory->create(time_previous, time_current, field_ids, block_data, &model_data);
    block_material_interface->ComputeStress();
    watch_internal_details.pop_region_and_report_time();

    // Stress divergence
    watch_internal_details.push_region("Stress divergence calculation");

    for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
      int block_id = block_it->first;
      nimble_kokkos::Block& block = block_it->second;
      nimble::Element* element_d = block.GetDeviceElement();
      int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      int num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);

      nimble_kokkos::DeviceElementConnectivityView elem_conn_d = block.GetDeviceElementConnectivityView();
      nimble_kokkos::DeviceVectorNodeGatheredView gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
      nimble_kokkos::DeviceVectorNodeGatheredView gathered_displacement_block_d = gathered_displacement_d.at(block_index);
      nimble_kokkos::DeviceVectorNodeGatheredView gathered_internal_force_block_d = gathered_internal_force_d.at(block_index);

      nimble_kokkos::DeviceSymTensorIntPtView stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id,
                                                                                                                    field_ids.stress,
                                                                                                                    nimble::STEP_NP1);

      // COMPUTE NODAL FORCES
      Kokkos::parallel_for("Force", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
        nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
        nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d = Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
        nimble_kokkos::DeviceSymTensorIntPtSubView element_stress_step_np1_d = Kokkos::subview(stress_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
        nimble_kokkos::DeviceVectorNodeGatheredSubView element_internal_force_d = Kokkos::subview(gathered_internal_force_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
        element_d->ComputeNodalForces(element_reference_coordinate_d,
                                      element_displacement_d,
                                      element_stress_step_np1_d,
                                      element_internal_force_d);
      });

      model_data.ScatterVectorNodeData(field_ids.internal_force,
                                       num_elem_in_block,
                                       num_nodes_per_elem,
                                       elem_conn_d,
                                       gathered_internal_force_block_d);
    } // loop over blocks
    watch_internal_details.pop_region_and_report_time();

    Kokkos::deep_copy(internal_force_h, internal_force_d);

    total_internal_force_time += watch_internal.pop_region_and_report_time();

    // Perform a reduction to obtain correct values on MPI boundaries
    watch_internal.push_region("MPI reduction");
    myVectorCommunicator->VectorReduction(mpi_vector_dim, internal_force_h);
    total_vector_reduction_time += watch_internal.pop_region_and_report_time();
    //

    // Evaluate the contact force
    if (contact_enabled) {
      watch_internal_details.push_region("Contact");
      //
      Kokkos::deep_copy(contact_force_d, (double)(0.0));
      contact_manager.ApplyDisplacements(displacement_d);
      contact_manager.ComputeContactForce(step+1, is_output_step);
      //
      contact_manager.GetForces(contact_force_d);
      Kokkos::deep_copy(contact_force_h, contact_force_d);
      total_contact_time += watch_internal_details.pop_region_and_report_time();
      // Perform a reduction to obtain correct values on MPI boundaries
      {
        watch_internal_details.push_region("MPI reduction");
        myVectorCommunicator->VectorReduction(mpi_vector_dim, contact_force_h);
        total_vector_reduction_time +=
            watch_internal_details.pop_region_and_report_time();
      }
      //
      auto tmpNum = contact_manager.numActiveContactFaces();
      if (tmpNum)
        contactInfo.insert(std::make_pair(step, tmpNum));
      //
//      if (contact_visualization && is_output_step)
//        contact_manager.ContactVisualizationWriteStep(time_current);
    }

    // fill acceleration vector A^{n+1} = M^{-1} ( F^{n} + b^{n} )
    watch_internal.push_region("Central difference");
    if (contact_enabled) {
      for (int i = 0; i < num_nodes; ++i) {
        acceleration_h(i, 0) = (1.0 / lumped_mass(i, 0)) *
                               (internal_force_h(i, 0) + contact_force_h(i, 0));
        acceleration_h(i, 1) = (1.0 / lumped_mass(i, 0)) *
                               (internal_force_h(i, 1) + contact_force_h(i, 1));
        acceleration_h(i, 2) = (1.0 / lumped_mass(i, 0)) *
                               (internal_force_h(i, 2) + contact_force_h(i, 2));
      }
    }
    else {
      for (int i = 0; i < num_nodes; ++i) {
        acceleration_h(i, 0) = (1.0 / lumped_mass(i, 0)) * internal_force_h(i, 0);
        acceleration_h(i, 1) = (1.0 / lumped_mass(i, 0)) * internal_force_h(i, 1);
        acceleration_h(i, 2) = (1.0 / lumped_mass(i, 0)) * internal_force_h(i, 2);
      }
    }

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    for (int i=0 ; i<num_nodes ; ++i) {
      velocity_h(i,0) += half_delta_time * acceleration_h(i,0);
      velocity_h(i,1) += half_delta_time * acceleration_h(i,1);
      velocity_h(i,2) += half_delta_time * acceleration_h(i,2);
    }
    total_update_avu_time += watch_internal.pop_region_and_report_time();

    if (is_output_step) {
      //
      watch_internal.push_region("Output");

      boundary_condition_manager.ApplyKinematicBC(time_current, time_previous, reference_coordinate_h, displacement_h, velocity_h);
      Kokkos::deep_copy(displacement_d, displacement_h);
      Kokkos::deep_copy(velocity_d, velocity_h);
      exodus_output_manager.ComputeElementData(mesh, model_data, blocks, gathered_reference_coordinate_d, gathered_displacement_d);
      std::vector<double> glbl_data;
      std::vector< std::vector<double> > const &node_data_output = exodus_output_manager.GetNodeDataForOutput(model_data);
      std::map<int, std::vector< std::vector<double> > > const &elem_data_output = exodus_output_manager.GetElementDataForOutput(model_data);
      std::map<int, std::vector< std::vector<double> > > drvd_elem_data;
      exodus_output.WriteStep(time_current,
                              glbl_data,
                              node_data_output,
                              elem_data_labels_for_output,
                              elem_data_output,
                              derived_elem_data_labels,
                              drvd_elem_data);
      //
      //
      if (contact_visualization) {
        contact_manager.ContactVisualizationWriteStep(time_current);
      }

      total_exodus_write_time += watch_internal.pop_region_and_report_time();
    }

    watch_internal.push_region("Copy field data new to old");
    // Copy STEP_NP1 data to STEP_N
    for (block_index = 0, block_it = blocks.begin(); block_it != blocks.end(); block_index++, block_it++) {
      int block_id = block_it->first;
      auto deformation_gradient_step_n_d = model_data.GetDeviceFullTensorIntegrationPointData(
          block_id, field_ids.deformation_gradient, nimble::STEP_N);
      auto unrotated_stress_step_n_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id,
                                                                                         field_ids.unrotated_stress,
                                                                                         nimble::STEP_N);
      auto stress_step_n_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id, field_ids.stress,
                                                                               nimble::STEP_N);
      auto deformation_gradient_step_np1_d = model_data.GetDeviceFullTensorIntegrationPointData(
          block_id, field_ids.deformation_gradient, nimble::STEP_NP1);
      auto unrotated_stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id,
                                                                                           field_ids.unrotated_stress,
                                                                                           nimble::STEP_NP1);
      auto stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id, field_ids.stress,
                                                                                 nimble::STEP_NP1);
      Kokkos::deep_copy(deformation_gradient_step_n_d, deformation_gradient_step_np1_d);
      Kokkos::deep_copy(unrotated_stress_step_n_d, unrotated_stress_step_np1_d);
      Kokkos::deep_copy(stress_step_n_d, stress_step_np1_d);
    }
    watch_internal.pop_region_and_report_time();

  } // loop over time steps
  double total_simulation_time = watch_simulation.pop_region_and_report_time();

  for (int irank = 0; irank < num_ranks; ++irank) {
#ifdef NIMBLE_HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if ((my_rank == irank) && (!contactInfo.empty())) {
      std::cout << " Rank " << irank << " has " << contactInfo.size()
                << " contact entries "
                << "(out of " << num_load_steps << " time steps)."<< std::endl;
      std::cout.flush();
    }
#ifdef NIMBLE_HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  if (my_rank == 0 && parser.WriteTimingDataFile()) {
//    double tcontact = total_contact_applyd + total_arborx_time + total_contact_getf;
    nimble::TimingInfo timing_writer{
        num_ranks,
        nimble::quanta::stopwatch::get_microsecond_timestamp(),
        total_simulation_time,
        total_internal_force_time,
        total_contact_time,
        total_exodus_write_time,
        total_vector_reduction_time
    };
    timing_writer.BinaryWrite();
  }

  if (my_rank == 0) {
    std::cout << " Total Time Loop = " << total_simulation_time << "\n";
    std::cout << " --- Internal Forces = " << total_internal_force_time << "\n";
    if (contact_enabled) {
//      double tcontact = total_contact_applyd + total_arborx_time + total_contact_getf;
      std::cout << " --- Contact = " << total_contact_time << "\n";
//      std::cout << " --- >>> Apply displ. = " << total_contact_applyd << "\n";
//      std::cout << " --- >>> Search / Project / Enforce = " << total_arborx_time << "\n";
      auto list_timers = contact_manager.getTimers();
      for (const auto& st_pair : list_timers)
        std::cout << " --- >>> >>> " << st_pair.first << " = " << st_pair.second << "\n";
      std::cout << " --- >>> Get Forces = " << total_contact_getf << "\n";
    }
    std::cout << " --- Exodus Write = " << total_exodus_write_time << "\n";
    std::cout << " --- Update AVU = " << total_update_avu_time << "\n";
    std::cout << " --- Vector Reduction = " << total_vector_reduction_time << "\n";
  }

  return 0;
}

}


}
