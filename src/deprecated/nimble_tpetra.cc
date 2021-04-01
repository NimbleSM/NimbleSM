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
#include <iomanip>
#include <iostream>
#include <limits>

#include "nimble_boundary_condition_manager.h"
#include "nimble_contact_manager.h"
#include "nimble_data_manager.h"
#include "nimble_exodus_output.h"
#include "nimble_genesis_mesh.h"
#include "nimble_linear_solver.h"
#include "nimble_material_factory.h"
#include "nimble_mesh_utils.h"
#include "nimble_parser.h"
#include "nimble_tpetra_utils.h"
#include "nimble_utils.h"
#include "nimble_version.h"
#include "nimble_view.h"

#ifdef NIMBLE_HAVE_EXTRAS
#include "nimble_extras_material_factory.h"
#endif

int
ExplicitTimeIntegrator(
    nimble::Parser&                   parser,
    nimble::GenesisMesh&              mesh,
    nimble::DataManager&              data_manager,
    nimble::BoundaryConditionManager& boundary_condition_manager,
    nimble::ExodusOutput&             exodus_output,
    int                               num_mpi_ranks,
    int                               my_mpi_rank,
    nimble::comm_type const&          comm);

int
QuasistaticTimeIntegrator(
    nimble::Parser&                   parser,
    nimble::GenesisMesh&              mesh,
    nimble::DataManager&              data_manager,
    nimble::BoundaryConditionManager& boundary_condition_manager,
    nimble::ExodusOutput&             exodus_output,
    int                               num_mpi_ranks,
    int                               my_mpi_rank,
    nimble::comm_type const&          comm);

int
main(int argc, char* argv[])
{
  bool using_mpi     = false;
  int  num_mpi_ranks = 1;
  int  my_mpi_rank   = 0;

#ifdef NIMBLE_HAVE_TRILINOS
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  auto               comm = Tpetra::getDefaultComm();
  using_mpi               = true;
  num_mpi_ranks           = comm->getSize();
  my_mpi_rank             = comm->getRank();
#endif

#ifdef NIMBLE_HAVE_KOKKOS
  Kokkos::initialize(argc, argv);
#endif

#ifdef NIMBLE_HAVE_EXTRAS
  using MaterialFactoryType = nimble::ExtrasMaterialFactory;
#else
  using MaterialFactoryType = nimble::MaterialFactory;
#endif

  int status = 0;

  // Banner
  if (!using_mpi) {
    std::cout << "\n-- NimbleSM" << std::endl;
    std::cout << "-- version " << nimble::NimbleVersion() << "\n" << std::endl;
    if (argc != 2) {
      std::cout << "Usage:  NimbleSM_Tpetra <input_deck.in>\n" << std::endl;
      exit(1);
    }
  } else if (using_mpi && my_mpi_rank == 0) {
    std::cout << "\n-- NimbleSM\n" << std::endl;
    if (argc != 2) {
      std::cout << "Usage:  mpirun -np NP NimbleSM_Tpetra <input_deck.in>\n" << std::endl;
      exit(1);
    }
    std::cout << "NimbleSM_Tpetra initialized on " << num_mpi_ranks << " mpi rank(s).\n" << std::endl;
  }

  std::string    input_deck_name = argv[1];
  nimble::Parser parser;
  parser.Initialize(input_deck_name);

  // Read the mesh
  std::string genesis_file_name     = nimble::IOFileName(parser.GenesisFileName(), "g", "", my_mpi_rank, num_mpi_ranks);
  std::string rve_genesis_file_name = nimble::IOFileName(parser.RVEGenesisFileName(), "g");
  nimble::GenesisMesh mesh;
  mesh.ReadFile(genesis_file_name);
  nimble::GenesisMesh rve_mesh;
  if (rve_genesis_file_name != "none") { rve_mesh.ReadFile(rve_genesis_file_name); }
  std::string tag = "serial";
  if (using_mpi) { tag = "tpetra"; }
  std::string output_exodus_name = nimble::IOFileName(parser.ExodusFileName(), "e", tag, my_mpi_rank, num_mpi_ranks);
  int         dim                = mesh.GetDim();
  int         num_nodes          = mesh.GetNumNodes();
  int         num_blocks         = mesh.GetNumBlocks();

  nimble::DataManager data_manager;
  nimble::ModelData&  macroscale_data = data_manager.GetMacroScaleData();
  macroscale_data.SetDimension(dim);

  // Global data
  int                      num_global_data = 0;
  std::vector<std::string> global_data_labels(num_global_data);
  std::vector<double>      global_data(num_global_data);

  int lumped_mass_field_id = macroscale_data.AllocateNodeData(nimble::SCALAR, "lumped_mass", num_nodes);
  int reference_coordinate_field_id =
      macroscale_data.AllocateNodeData(nimble::VECTOR, "reference_coordinate", num_nodes);
  int displacement_field_id = macroscale_data.AllocateNodeData(nimble::VECTOR, "displacement", num_nodes);
  int displacement_fluctuation_field_id =
      macroscale_data.AllocateNodeData(nimble::VECTOR, "displacement_fluctuation", num_nodes);
  int velocity_field_id       = macroscale_data.AllocateNodeData(nimble::VECTOR, "velocity", num_nodes);
  int acceleration_field_id   = macroscale_data.AllocateNodeData(nimble::VECTOR, "acceleration", num_nodes);
  int internal_force_field_id = macroscale_data.AllocateNodeData(nimble::VECTOR, "internal_force", num_nodes);
  int external_force_field_id = macroscale_data.AllocateNodeData(nimble::VECTOR, "external_force", num_nodes);
  // Blocks
  std::map<int, nimble::Block>&          blocks = macroscale_data.GetBlocks();
  std::map<int, nimble::Block>::iterator block_it;
  std::vector<int>                       block_ids = mesh.GetBlockIds();
  for (int i = 0; i < num_blocks; i++) {
    int                               block_id                  = block_ids[i];
    std::string const&                macro_material_parameters = parser.GetMacroscaleMaterialParameters(block_id);
    std::map<int, std::string> const& rve_material_parameters   = parser.GetMicroscaleMaterialParameters();
    std::string                       rve_bc_strategy           = parser.GetMicroscaleBoundaryConditionStrategy();
    MaterialFactoryType               factory;
    blocks[block_id] = nimble::Block();
    blocks[block_id].Initialize(macro_material_parameters, rve_material_parameters, rve_mesh, rve_bc_strategy, factory);
    std::vector<std::pair<std::string, nimble::Length>> data_labels_and_lengths;
    blocks[block_id].GetDataLabelsAndLengths(data_labels_and_lengths);
    macroscale_data.DeclareElementData(block_id, data_labels_and_lengths);
  }
  std::map<int, int> num_elem_in_each_block = mesh.GetNumElementsInBlock();
  macroscale_data.AllocateElementData(num_elem_in_each_block);
  macroscale_data.SpecifyOutputFields(parser.GetOutputFieldString());
  std::map<int, std::vector<std::string>> const& elem_data_labels = macroscale_data.GetElementDataLabels();
  std::map<int, std::vector<std::string>> const& elem_data_labels_for_output =
      macroscale_data.GetElementDataLabelsForOutput();
  std::map<int, std::vector<std::string>> const& derived_elem_data_labels =
      macroscale_data.GetDerivedElementDataLabelsForOutput();

  // Initialize the element data
  std::vector<int> rve_output_elem_ids = parser.MicroscaleOutputElementIds();
  for (block_it = blocks.begin(); block_it != blocks.end(); block_it++) {
    int                     block_id          = block_it->first;
    int                     num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    std::vector<int> const& elem_global_ids   = mesh.GetElementGlobalIdsInBlock(block_id);
    nimble::Block&          block             = block_it->second;
    std::vector<double>&    elem_data_n       = macroscale_data.GetElementDataOld(block_id);
    std::vector<double>&    elem_data_np1     = macroscale_data.GetElementDataNew(block_id);
    MaterialFactoryType     factory;
    block.InitializeElementData(
        num_elem_in_block,
        elem_global_ids,
        rve_output_elem_ids,
        elem_data_labels.at(block_id),
        derived_elem_data_labels.at(block_id),
        elem_data_n,
        elem_data_np1,
        factory,
        data_manager);
  }

  // Initialize the initial- and boundary-condition manager
  std::map<int, std::string> const&      node_set_names          = mesh.GetNodeSetNames();
  std::map<int, std::vector<int>> const& node_sets               = mesh.GetNodeSets();
  std::vector<std::string> const&        bc_strings              = parser.GetBoundaryConditionStrings();
  std::string                            time_integration_scheme = parser.TimeIntegrationScheme();
  nimble::BoundaryConditionManager       bc;
  bc.Initialize(node_set_names, node_sets, bc_strings, dim, time_integration_scheme);
  // Initialize the output file
  nimble::ExodusOutput exodus_output;
  exodus_output.Initialize(output_exodus_name, mesh);
  std::vector<std::string> const& node_data_labels_for_output = macroscale_data.GetNodeDataLabelsForOutput();
  exodus_output.InitializeDatabase(
      mesh, global_data_labels, node_data_labels_for_output, elem_data_labels_for_output, derived_elem_data_labels);

  const double* const ref_coord_x          = mesh.GetCoordinatesX();
  const double* const ref_coord_y          = mesh.GetCoordinatesY();
  const double* const ref_coord_z          = mesh.GetCoordinatesZ();
  double*             reference_coordinate = macroscale_data.GetNodeData(reference_coordinate_field_id);
  for (int i = 0; i < num_nodes; i++) {
    reference_coordinate[3 * i]     = ref_coord_x[i];
    reference_coordinate[3 * i + 1] = ref_coord_y[i];
    reference_coordinate[3 * i + 2] = ref_coord_z[i];
  }

  if (time_integration_scheme == "explicit") {
    ExplicitTimeIntegrator(parser, mesh, data_manager, bc, exodus_output, num_mpi_ranks, my_mpi_rank, comm);
  } else if (time_integration_scheme == "quasistatic") {
    QuasistaticTimeIntegrator(parser, mesh, data_manager, bc, exodus_output, num_mpi_ranks, my_mpi_rank, comm);
  }

#ifdef NIMBLE_HAVE_KOKKOS
  Kokkos::finalize();
#endif

  return status;
}

int
ExplicitTimeIntegrator(
    nimble::Parser&                   parser,
    nimble::GenesisMesh&              mesh,
    nimble::DataManager&              data_manager,
    nimble::BoundaryConditionManager& bc,
    nimble::ExodusOutput&             exodus_output,
    int                               num_mpi_ranks,
    int                               my_mpi_rank,
    nimble::comm_type const&          comm)
{
  int dim        = mesh.GetDim();
  int num_nodes  = mesh.GetNumNodes();
  int num_blocks = mesh.GetNumBlocks();

  int              contact_block_id_1   = 0;
  int              contact_block_id_2   = 0;
  double*          contact_data_block_1 = 0;
  double*          contact_data_block_2 = 0;
  std::vector<int> contact_results_block_1;
  std::vector<int> contact_results_block_2;
  bool             has_contact = parser.HasContact();
  if (has_contact) {
    std::stringstream ss(parser.ContactString());
    std::string       block_name;
    ss >> block_name;
    contact_block_id_1 = mesh.GetBlockId(block_name);
    contact_results_block_1.resize(mesh.GetNumElementsInBlock(contact_block_id_1));
    ss >> block_name;
    contact_block_id_2 = mesh.GetBlockId(block_name);
    contact_results_block_2.resize(mesh.GetNumElementsInBlock(contact_block_id_2));
  }

  std::vector<int> global_node_ids(num_nodes);
  int const* const global_node_ids_ptr = mesh.GetNodeGlobalIds();
  for (int n = 0; n < num_nodes; ++n) { global_node_ids[n] = global_node_ids_ptr[n]; }

  nimble::TpetraContainer tpetra_container;
  tpetra_container.Initialize(mesh, global_node_ids, comm);

  int                 num_global_data = 0;
  std::vector<double> global_data(num_global_data);

  nimble::ModelData& macroscale_data = data_manager.GetMacroScaleData();

  int lumped_mass_field_id          = macroscale_data.GetFieldId("lumped_mass");
  int reference_coordinate_field_id = macroscale_data.GetFieldId("reference_coordinate");
  int displacement_field_id         = macroscale_data.GetFieldId("displacement");
  int velocity_field_id             = macroscale_data.GetFieldId("velocity");
  int acceleration_field_id         = macroscale_data.GetFieldId("acceleration");
  int internal_force_field_id       = macroscale_data.GetFieldId("internal_force");
  int external_force_field_id       = macroscale_data.GetFieldId("external_force");

  int status = 0;
  // Set up the global vectors
  unsigned int num_unknowns         = num_nodes * mesh.GetDim();
  double*      lumped_mass          = macroscale_data.GetNodeData(lumped_mass_field_id);
  double*      reference_coordinate = macroscale_data.GetNodeData(reference_coordinate_field_id);
  double*      displacement         = macroscale_data.GetNodeData(displacement_field_id);
  double*      velocity             = macroscale_data.GetNodeData(velocity_field_id);
  double*      acceleration         = macroscale_data.GetNodeData(acceleration_field_id);
  double*      internal_force       = macroscale_data.GetNodeData(internal_force_field_id);
  double*      external_force       = macroscale_data.GetNodeData(external_force_field_id);

  std::map<int, std::vector<std::string>> const& elem_data_labels = macroscale_data.GetElementDataLabels();
  std::map<int, std::vector<std::string>> const& elem_data_labels_for_output =
      macroscale_data.GetElementDataLabelsForOutput();
  std::map<int, std::vector<std::string>> const& derived_elem_data_labels =
      macroscale_data.GetDerivedElementDataLabelsForOutput();

  // Computed the lumped mass matrix (diagonal matrix) and the critical time
  // step
  double                                 critical_time_step = std::numeric_limits<double>::max();
  std::map<int, nimble::Block>&          blocks             = macroscale_data.GetBlocks();
  std::map<int, nimble::Block>::iterator block_it;
  for (block_it = blocks.begin(); block_it != blocks.end(); block_it++) {
    int            block_id          = block_it->first;
    int            num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    int const*     elem_conn         = mesh.GetConnectivity(block_id);
    nimble::Block& block             = block_it->second;
    block.ComputeLumpedMassMatrix(reference_coordinate, num_elem_in_block, elem_conn, lumped_mass);
    double block_critical_time_step =
        block.ComputeCriticalTimeStep(reference_coordinate, displacement, num_elem_in_block, elem_conn);
    if (block_critical_time_step < critical_time_step) { critical_time_step = block_critical_time_step; }
  }
  tpetra_container.VectorReduction(macroscale_data, "lumped_mass");

  double time_current(0.0), time_previous(0.0);
  double final_time = parser.FinalTime();
  double delta_time, half_delta_time;
  int    num_load_steps   = parser.NumLoadSteps();
  int    output_frequency = parser.OutputFrequency();

  bc.ApplyInitialConditions(Viewify(reference_coordinate, 3), Viewify(velocity, 3));
  // Allow prescribed velocities to trump initial velocities
  bc.ApplyKinematicBC(0.0, 0.0, Viewify(reference_coordinate, 3), Viewify(displacement, 3), Viewify(velocity, 3));

  // For explicit dynamics, the macroscale model is never treated as an RVE
  // so rve_macroscale_deformation_gradient will always be the identity matrix
  std::vector<double> rve_macroscale_deformation_gradient(dim * dim, 0.0);
  for (int i = 0; i < dim; i++) { rve_macroscale_deformation_gradient[i] = 1.0; }

  std::map<int, std::vector<std::vector<double>>> derived_elem_data;
  for (block_it = blocks.begin(); block_it != blocks.end(); block_it++) {
    int                        block_id          = block_it->first;
    nimble::Block&             block             = block_it->second;
    int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    int const*                 elem_conn         = mesh.GetConnectivity(block_id);
    std::vector<double> const& elem_data_np1     = macroscale_data.GetElementDataNew(block_id);
    derived_elem_data[block_id]                  = std::vector<std::vector<double>>();
    block.ComputeDerivedElementData(
        reference_coordinate,
        displacement,
        num_elem_in_block,
        elem_conn,
        elem_data_labels.at(block_id).size(),
        elem_data_np1,
        derived_elem_data_labels.at(block_id).size(),
        derived_elem_data.at(block_id));
  }
  std::vector<std::vector<double>> node_data_for_output;
  macroscale_data.GetNodeDataForOutput(node_data_for_output);
  std::map<int, std::vector<std::vector<double>>> elem_data_for_output;
  macroscale_data.GetElementDataForOutput(elem_data_for_output);
  exodus_output.WriteStep(
      time_current,
      global_data,
      node_data_for_output,
      elem_data_labels_for_output,
      elem_data_for_output,
      derived_elem_data_labels,
      derived_elem_data);

  double user_specified_time_step = final_time / num_load_steps;
  if (my_mpi_rank == 0) {
    std::cout << "User specified time step:              " << user_specified_time_step << std::endl;
    std::cout << "Approximate maximum stable time step:  " << critical_time_step << "\n" << std::endl;
    if (user_specified_time_step > critical_time_step) {
      std::cout << "**** WARNING:  The user specified time step exceeds the "
                   "computed maximum stable time step.\n"
                << std::endl;
    }
    std::cout << "Explicit time integration:\n    0% complete" << std::endl;
  }

  for (int step = 0; step < num_load_steps; step++) {
    if (my_mpi_rank == 0) {
      if (10 * (step + 1) % num_load_steps == 0 && step != num_load_steps - 1) {
        std::cout << "   " << static_cast<int>(100.0 * static_cast<double>(step + 1) / num_load_steps) << "% complete"
                  << std::endl
                  << std::flush;
      } else if (step == num_load_steps - 1) {
        std::cout << "  100% complete\n" << std::endl << std::flush;
      }
    }
    bool is_output_step = false;
    if (step % output_frequency == 0 || step == num_load_steps - 1) { is_output_step = true; }

    time_previous = time_current;
    time_current += final_time / num_load_steps;
    delta_time      = time_current - time_previous;
    half_delta_time = 0.5 * delta_time;

    // V^{n+1/2} = V^{n} + (dt/2) * A^{n}
    for (int i = 0; i < num_unknowns; ++i) { velocity[i] += half_delta_time * acceleration[i]; }

    bc.ApplyKinematicBC(
        time_current, time_previous, Viewify(reference_coordinate, 3), Viewify(displacement, 3), Viewify(velocity, 3));

    // Evaluate external body forces
    for (int i = 0; i < num_unknowns; ++i) { external_force[i] = 0.0; }

    // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
    for (int i = 0; i < num_unknowns; ++i) { displacement[i] += delta_time * velocity[i]; }

    // Evaluate the internal force
    for (int i = 0; i < num_unknowns; ++i) { internal_force[i] = 0.0; }
    for (block_it = blocks.begin(); block_it != blocks.end(); block_it++) {
      int                        block_id          = block_it->first;
      int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      int const*                 elem_conn         = mesh.GetConnectivity(block_id);
      std::vector<int> const&    elem_global_ids   = mesh.GetElementGlobalIdsInBlock(block_id);
      nimble::Block&             block             = block_it->second;
      std::vector<double> const& elem_data_n       = macroscale_data.GetElementDataOld(block_id);
      std::vector<double>&       elem_data_np1     = macroscale_data.GetElementDataNew(block_id);
      block.ComputeInternalForce(
          reference_coordinate,
          displacement,
          velocity,
          rve_macroscale_deformation_gradient.data(),
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
          is_output_step);
    }
    tpetra_container.VectorReduction(macroscale_data, "internal_force");

    // fill acceleration vector A^{n+1} = M^{-1} ( F^{n} + b^{n} )
    for (int i = 0; i < num_unknowns; ++i) {
      acceleration[i] = (1.0 / lumped_mass[i / 3]) * (internal_force[i] + external_force[i]);
    }

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    for (int i = 0; i < num_unknowns; ++i) { velocity[i] += half_delta_time * acceleration[i]; }

    if (is_output_step) {
      for (block_it = blocks.begin(); block_it != blocks.end(); block_it++) {
        int                        block_id          = block_it->first;
        nimble::Block&             block             = block_it->second;
        int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
        int const*                 elem_conn         = mesh.GetConnectivity(block_id);
        std::vector<double> const& elem_data_np1     = macroscale_data.GetElementDataNew(block_id);
        block.ComputeDerivedElementData(
            reference_coordinate,
            displacement,
            num_elem_in_block,
            elem_conn,
            elem_data_labels.at(block_id).size(),
            elem_data_np1,
            derived_elem_data_labels.at(block_id).size(),
            derived_elem_data.at(block_id));
      }

      bc.ApplyKinematicBC(
          time_current,
          time_previous,
          Viewify(reference_coordinate, 3),
          Viewify(displacement, 3),
          Viewify(velocity, 3));

      // Write output
      macroscale_data.GetNodeDataForOutput(node_data_for_output);
      macroscale_data.GetElementDataForOutput(elem_data_for_output);
      exodus_output.WriteStep(
          time_current,
          global_data,
          node_data_for_output,
          elem_data_labels_for_output,
          elem_data_for_output,
          derived_elem_data_labels,
          derived_elem_data);
    }

    macroscale_data.SwapStates();
  }

  return status;
}

int
QuasistaticTimeIntegrator(
    nimble::Parser&                   parser,
    nimble::GenesisMesh&              mesh,
    nimble::DataManager&              data_manager,
    nimble::BoundaryConditionManager& bc,
    nimble::ExodusOutput&             exodus_output,
    int                               num_mpi_ranks,
    int                               my_mpi_rank,
    nimble::comm_type const&          comm)
{
  int status = 0;

  int              dim             = mesh.GetDim();
  int              num_nodes       = mesh.GetNumNodes();
  int              num_blocks      = mesh.GetNumBlocks();
  const int* const global_node_ids = mesh.GetNodeGlobalIds();

  // Store various mappings from global to local node ids
  // Things get a bit tricky for periodic boundary conditions,
  // where dof are condensed out of the linear system
  std::vector<int>                linear_system_global_node_ids;
  std::map<int, std::vector<int>> map_from_linear_system;
  if (bc.IsPeriodicRVEProblem()) {
    int rve_corner_node_id(-1);
    mesh.CreatePeriodicRVELinearSystemMap(
        global_node_ids, linear_system_global_node_ids, map_from_linear_system, rve_corner_node_id);
    bc.CreateRVEFixedCornersBoundaryConditions(rve_corner_node_id);
  } else {
    for (int n = 0; n < num_nodes; ++n) {
      int global_node_id = global_node_ids[n];
      linear_system_global_node_ids.push_back(global_node_id);
      map_from_linear_system[global_node_id] = std::vector<int>();
      map_from_linear_system[global_node_id].push_back(n);
    }
  }
  int linear_system_num_nodes    = static_cast<int>(map_from_linear_system.size());
  int linear_system_num_unknowns = linear_system_num_nodes * dim;

  nimble::TpetraContainer tpetra_container;
  tpetra_container.Initialize(mesh, linear_system_global_node_ids, comm);

  int                 num_global_data = 0;
  std::vector<double> global_data(num_global_data);

  nimble::ModelData& macroscale_data = data_manager.GetMacroScaleData();

  int reference_coordinate_field_id     = macroscale_data.GetFieldId("reference_coordinate");
  int displacement_field_id             = macroscale_data.GetFieldId("displacement");
  int displacement_fluctuation_field_id = macroscale_data.GetFieldId("displacement_fluctuation");
  int velocity_field_id                 = macroscale_data.GetFieldId("velocity");
  int internal_force_field_id           = macroscale_data.GetFieldId("internal_force");
  int external_force_field_id           = macroscale_data.GetFieldId("external_force");

  // Set up the global vectors
  unsigned int num_unknowns             = num_nodes * mesh.GetDim();
  double*      reference_coordinate     = macroscale_data.GetNodeData(reference_coordinate_field_id);
  double*      physical_displacement    = macroscale_data.GetNodeData(displacement_field_id);
  double*      displacement_fluctuation = macroscale_data.GetNodeData(displacement_fluctuation_field_id);
  double*      velocity                 = macroscale_data.GetNodeData(velocity_field_id);
  double*      internal_force           = macroscale_data.GetNodeData(internal_force_field_id);
  double*      external_force           = macroscale_data.GetNodeData(external_force_field_id);

  std::vector<double> residual_vector(linear_system_num_unknowns, 0.0);
  std::vector<double> linear_solver_solution(linear_system_num_unknowns, 0.0);
  double*             displacement = physical_displacement;
  if (bc.IsPeriodicRVEProblem()) { displacement = displacement_fluctuation; }

  nimble::CRSMatrixContainer tangent_stiffness;
  nimble::CGScratchSpace     cg_scratch;
  std::vector<int>           i_index, j_index;
  nimble::DetermineTangentMatrixNonzeroStructure(mesh, linear_system_global_node_ids, i_index, j_index);
  tangent_stiffness.AllocateNonzeros(i_index, j_index);
  if (my_mpi_rank == 0) {
    std::cout << "Number of nonzeros in tangent stiffness matrix = " << tangent_stiffness.NumNonzeros() << "\n"
              << std::endl;
  }

  tpetra_container.AllocateTangentStiffnessMatrix(mesh);
  if (my_mpi_rank == 0) {
    std::cout << "Number of nonzeros in the crs tangent stiffness matrix = "
              << tpetra_container.TangentStiffnessMatrixNumNonzeros() << "\n"
              << std::endl;
  }

  std::map<int, std::vector<std::string>> const& elem_data_labels = macroscale_data.GetElementDataLabels();
  std::map<int, std::vector<std::string>> const& elem_data_labels_for_output =
      macroscale_data.GetElementDataLabelsForOutput();
  std::map<int, std::vector<std::string>> const& derived_elem_data_labels =
      macroscale_data.GetDerivedElementDataLabelsForOutput();

  double time_current(0.0), time_previous(0.0), delta_time(0.0);
  double final_time       = parser.FinalTime();
  int    num_load_steps   = parser.NumLoadSteps();
  int    output_frequency = parser.OutputFrequency();

  bc.ApplyKinematicBC(0.0, 0.0, Viewify(reference_coordinate, 3), Viewify(displacement, 3), Viewify(velocity, 3));

  std::vector<double> rve_macroscale_deformation_gradient(dim * dim, 0.0);
  std::vector<double> identity(dim * dim, 0.0);
  for (int i = 0; i < dim; i++) { identity[i] = 1.0; }
  std::vector<double> rve_center;
  if (bc.IsPeriodicRVEProblem()) { rve_center = mesh.BoundingBoxCenter(); }

  std::map<int, nimble::Block>&                   blocks = macroscale_data.GetBlocks();
  std::map<int, std::vector<std::vector<double>>> derived_elem_data;
  for (auto& block_it : blocks) {
    int                        block_id          = block_it.first;
    nimble::Block&             block             = block_it.second;
    int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    int const*                 elem_conn         = mesh.GetConnectivity(block_id);
    std::vector<double> const& elem_data_np1     = macroscale_data.GetElementDataNew(block_id);
    derived_elem_data[block_id]                  = std::vector<std::vector<double>>();
    block.ComputeDerivedElementData(
        reference_coordinate,
        displacement,
        num_elem_in_block,
        elem_conn,
        elem_data_labels.at(block_id).size(),
        elem_data_np1,
        derived_elem_data_labels.at(block_id).size(),
        derived_elem_data.at(block_id));
  }
  std::vector<std::vector<double>> node_data_for_output;
  macroscale_data.GetNodeDataForOutput(node_data_for_output);
  std::map<int, std::vector<std::vector<double>>> elem_data_for_output;
  macroscale_data.GetElementDataForOutput(elem_data_for_output);
  exodus_output.WriteStep(
      time_current,
      global_data,
      node_data_for_output,
      elem_data_labels_for_output,
      elem_data_for_output,
      derived_elem_data_labels,
      derived_elem_data);

  if (my_mpi_rank == 0) { std::cout << "Beginning quasistatic time integration:" << std::endl; }

  for (int step = 0; step < num_load_steps; step++) {
    time_previous = time_current;
    time_current += final_time / num_load_steps;
    delta_time = time_current - time_previous;

    bc.ApplyKinematicBC(
        time_current, time_previous, Viewify(reference_coordinate, 3), Viewify(displacement, 3), Viewify(velocity, 3));

    bc.GetRVEMacroscaleDeformationGradient(time_current, rve_macroscale_deformation_gradient.data());

    // Compute the residual, which is a norm of the (rearranged) nodal force
    // vector, with the dof associated with kinematic BC removed.
    for (int i = 0; i < num_unknowns; ++i) { internal_force[i] = 0.0; }
    for (auto& block_it : blocks) {
      int                        block_id          = block_it.first;
      int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      int const*                 elem_conn         = mesh.GetConnectivity(block_id);
      std::vector<int> const&    elem_global_ids   = mesh.GetElementGlobalIdsInBlock(block_id);
      nimble::Block&             block             = block_it.second;
      std::vector<double> const& elem_data_n       = macroscale_data.GetElementDataOld(block_id);
      std::vector<double>&       elem_data_np1     = macroscale_data.GetElementDataNew(block_id);
      block.ComputeInternalForce(
          reference_coordinate,
          displacement,
          velocity,
          rve_macroscale_deformation_gradient.data(),
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
          false);
    }
    for (int i = 0; i < linear_system_num_nodes * dim; i++) { residual_vector[i] = 0.0; }
    for (int n = 0; n < num_nodes; n++) {
      int ls_id = linear_system_global_node_ids[n];
      for (int dof = 0; dof < dim; dof++) {
        int local_index = n * dim + dof;
        int ls_index    = ls_id * dim + dof;
        if (ls_index < 0 || ls_index > residual_vector.size() - 1) {
          throw std::logic_error(
              "\nError:  Invalid index into residual vector in "
              "QuasistaticTimeIntegrator().\n");
        }
        residual_vector[ls_index] += -1.0 * internal_force[local_index];
      }
    }
    bc.ModifyRHSForKinematicBC(linear_system_global_node_ids.data(), residual_vector.data());
    double l2_norm(0.0), infinity_norm(0.0);
    for (int i = 0; i < linear_system_num_unknowns; i++) {
      double f = std::fabs(residual_vector[i]);
      l2_norm += f * f;
      if (f > infinity_norm) { infinity_norm = f; }
    }
    l2_norm                      = std::sqrt(l2_norm);
    double residual              = l2_norm + 20.0 * infinity_norm;
    double convergence_tolerance = 1.0e-6 * residual;
    int    iteration(0), max_nonlinear_iterations(50);

    if (my_mpi_rank == 0) {
      std::cout << "\nStep " << step + 1 << std::scientific << std::setprecision(3) << ", time " << time_current
                << ", delta_time " << delta_time << ", convergence tolerance " << convergence_tolerance << std::endl;
      std::cout << "  iteration " << iteration << ": residual = " << residual << std::endl;
    }

    bool is_output_step = false;
    if (step % output_frequency == 0 || step == num_load_steps - 1) { is_output_step = true; }

    while (residual > convergence_tolerance && iteration < max_nonlinear_iterations) {
      tangent_stiffness.SetAllValues(0.0);
      tpetra_container.TangentStiffnessMatrixSetScalar(0.0);
      for (auto& block_it : blocks) {
        int            block_id          = block_it.first;
        int            num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
        int const*     elem_conn         = mesh.GetConnectivity(block_id);
        nimble::Block& block             = block_it.second;
        block.ComputeTangentStiffnessMatrix(
            linear_system_num_unknowns,
            reference_coordinate,
            displacement,
            num_elem_in_block,
            elem_conn,
            linear_system_global_node_ids.data(),
            tangent_stiffness);
      }

      // For the dof with kinematic BC, zero out the rows and columns and put a
      // non-zero on the diagonal
      double diagonal_entry(0.0);
      for (int i = 0; i < linear_system_num_unknowns; ++i) { diagonal_entry += std::fabs(tangent_stiffness(i, i)); }
      diagonal_entry /= linear_system_num_unknowns;
      bc.ModifyTangentStiffnessMatrixForKinematicBC(
          linear_system_num_unknowns, linear_system_global_node_ids.data(), diagonal_entry, tangent_stiffness);

      bc.ModifyRHSForKinematicBC(linear_system_global_node_ids.data(), residual_vector.data());

      int num_cg_iterations(0);
      std::fill(linear_solver_solution.begin(), linear_solver_solution.end(), 0.0);
      bool success = nimble::CG_SolveSystem(
          tangent_stiffness, residual_vector.data(), cg_scratch, linear_solver_solution.data(), num_cg_iterations);
      if (!success) { throw std::logic_error("\nError:  CG linear solver failed to converge.\n"); }

      // todo: apply a line search

      // update the trial solution
      for (int n = 0; n < linear_system_num_nodes; n++) {
        std::vector<int> const& node_ids = map_from_linear_system[n];
        for (auto const& node_id : node_ids) {
          for (int dof = 0; dof < dim; dof++) {
            int local_index = node_id * dim + dof;
            int ls_index    = n * dim + dof;
            if (ls_index < 0 || ls_index > residual_vector.size() - 1) {
              throw std::logic_error(
                  "\nError:  Invalid index into residual vector in "
                  "QuasistaticTimeIntegrator().\n");
            }
            displacement[local_index] -= linear_solver_solution[ls_index];
          }
          // if we are solving a standard problem, then displacement is the
          // physical_displacement if we are solving a periodic RVE problem,
          // then displacement is the displacement_fluctuation, and we need to
          // set physical_displacement manually and include the macroscale
          // deformation gradient
          if (bc.IsPeriodicRVEProblem()) {
            std::vector<double> const& F   = rve_macroscale_deformation_gradient;
            int                        i_x = node_id * dim;
            int                        i_y = node_id * dim + 1;
            int                        i_z = node_id * dim + 2;
            // todo:  this is currently hard-coded to 3D
            physical_displacement[i_x] =
                displacement[i_x] + (F[K_F_XX] - identity[K_F_XX]) * (reference_coordinate[i_x] - rve_center.at(0)) +
                (F[K_F_XY] - identity[K_F_XY]) * (reference_coordinate[i_y] - rve_center.at(1)) +
                (F[K_F_XZ] - identity[K_F_XZ]) * (reference_coordinate[i_z] - rve_center.at(2));
            physical_displacement[i_y] =
                displacement[i_y] + (F[K_F_YX] - identity[K_F_YX]) * (reference_coordinate[i_x] - rve_center.at(0)) +
                (F[K_F_YY] - identity[K_F_YY]) * (reference_coordinate[i_y] - rve_center.at(1)) +
                (F[K_F_YZ] - identity[K_F_YZ]) * (reference_coordinate[i_z] - rve_center.at(2));
            physical_displacement[i_z] =
                displacement[i_z] + (F[K_F_ZX] - identity[K_F_ZX]) * (reference_coordinate[i_x] - rve_center.at(0)) +
                (F[K_F_ZY] - identity[K_F_ZY]) * (reference_coordinate[i_y] - rve_center.at(1)) +
                (F[K_F_ZZ] - identity[K_F_ZZ]) * (reference_coordinate[i_z] - rve_center.at(2));
          }
        }
      }

      // check convergence
      for (int i = 0; i < num_unknowns; ++i) { internal_force[i] = 0.0; }
      for (auto& block_it : blocks) {
        int                        block_id          = block_it.first;
        int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
        int const*                 elem_conn         = mesh.GetConnectivity(block_id);
        std::vector<int> const&    elem_global_ids   = mesh.GetElementGlobalIdsInBlock(block_id);
        nimble::Block&             block             = block_it.second;
        std::vector<double> const& elem_data_n       = macroscale_data.GetElementDataOld(block_id);
        std::vector<double>&       elem_data_np1     = macroscale_data.GetElementDataNew(block_id);
        block.ComputeInternalForce(
            reference_coordinate,
            displacement,
            velocity,
            rve_macroscale_deformation_gradient.data(),
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
            is_output_step);
      }
      for (int i = 0; i < linear_system_num_nodes * dim; i++) { residual_vector[i] = 0.0; }
      for (int n = 0; n < num_nodes; n++) {
        int ls_id = linear_system_global_node_ids[n];
        for (int dof = 0; dof < dim; dof++) {
          int local_index = n * dim + dof;
          int ls_index    = ls_id * dim + dof;
          if (ls_index < 0 || ls_index > residual_vector.size() - 1) {
            throw std::logic_error(
                "\nError:  Invalid index into residual vector in "
                "QuasistaticTimeIntegrator().\n");
          }
          residual_vector[ls_index] += -1.0 * internal_force[local_index];
        }
      }
      bc.ModifyRHSForKinematicBC(linear_system_global_node_ids.data(), residual_vector.data());
      double l2_norm(0.0), infinity_norm(0.0);
      for (int i = 0; i < num_nodes; i++) {
        double f = std::fabs(residual_vector[i]);
        l2_norm += f * f;
        if (f > infinity_norm) { infinity_norm = f; }
      }
      l2_norm  = std::sqrt(l2_norm);
      residual = l2_norm + 20.0 * infinity_norm;

      iteration += 1;

      if (my_mpi_rank == 0) {
        std::cout << "  iteration " << iteration << ": residual = " << residual
                  << ", linear cg iterations = " << num_cg_iterations << std::endl;
      }
    }  // while (residual > convergence_tolerance ... )

    if (iteration == max_nonlinear_iterations && my_mpi_rank == 0) {
      std::cout << "\n**** Nonlinear solver failed to converge in " << max_nonlinear_iterations << " iterations.\n"
                << std::endl;
      status = 1;
      return status;
    } else {
      if (is_output_step) {
        for (auto& block_it : blocks) {
          int                        block_id          = block_it.first;
          nimble::Block&             block             = block_it.second;
          int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
          int const*                 elem_conn         = mesh.GetConnectivity(block_id);
          std::vector<double> const& elem_data_np1     = macroscale_data.GetElementDataNew(block_id);
          block.ComputeDerivedElementData(
              reference_coordinate,
              displacement,
              num_elem_in_block,
              elem_conn,
              elem_data_labels.at(block_id).size(),
              elem_data_np1,
              derived_elem_data_labels.at(block_id).size(),
              derived_elem_data.at(block_id));
        }

        // Write output
        macroscale_data.GetNodeDataForOutput(node_data_for_output);
        macroscale_data.GetElementDataForOutput(elem_data_for_output);
        exodus_output.WriteStep(
            time_current,
            global_data,
            node_data_for_output,
            elem_data_labels_for_output,
            elem_data_for_output,
            derived_elem_data_labels,
            derived_elem_data);
      }

      macroscale_data.SwapStates();
    }
  }
  if (my_mpi_rank == 0) { std::cout << "\nComplete.\n" << std::endl; }

  return status;
}
