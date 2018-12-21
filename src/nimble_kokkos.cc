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
#include "nimble_parser.h"
#include "nimble_exodus_output.h"
#include "nimble_exodus_output_manager.h"
#include "nimble_boundary_condition_manager.h"
#include "nimble.mpi.utils.h"
#include "nimble_view.h"
#include "nimble_kokkos_defs.h"
#include "nimble_kokkos_data_manager.h"
#include "nimble_kokkos_block.h"
#include "nimble_contact.h"

#include <iostream>

void main_routine(int argc, char *argv[]) {

  int mpi_err;
  int num_mpi_ranks;
  int my_mpi_rank;

  mpi_err = MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_ranks);
  if (mpi_err != MPI_SUCCESS) {
    throw std::logic_error("\nError:  MPI_Comm_size() returned nonzero error code.\n");
  }
  mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);
  if (mpi_err != MPI_SUCCESS) {
    throw std::logic_error("\nError:  MPI_Comm_rank() returned nonzero error code.\n");
  }

  // Banner
  if (my_mpi_rank == 0) {

    std::string nimble_have_kokkos("false");
#ifdef NIMBLE_HAVE_KOKKOS
    nimble_have_kokkos = "true";
#endif
    std::string kokkos_enable_cuda("false");
#ifdef KOKKOS_ENABLE_CUDA
    kokkos_enable_cuda = "true";
#endif
    std::string kokkos_enable_cuda_uvm("false");
#ifdef KOKKOS_ENABLE_CUDA_UVM
    kokkos_enable_cuda_uvm = "true";
#endif

    std::cout << "\n-- NimbleSM_Kokkos" << std::endl;
    std::cout << "-- version " << nimble::NimbleVersion() << "\n" << std::endl;
    if (argc != 2) {
      std::cout << "Usage:  mpirun -np NP NimbleSM_Kokkos <input_deck.in>\n" << std::endl;
      Kokkos::finalize();
      MPI_Finalize();
      exit(1);
    }
    std::cout << "NimbleSM_Kokkos initialized on " << num_mpi_ranks << " mpi rank(s)." << std::endl;

    std::cout << "\nKokkos configuration:" << std::endl;
    std::cout << "  NIMBLE_HAVE_KOKKOS               " << nimble_have_kokkos << std::endl;
    std::cout << "  KOKKOS_ENABLE_CUDA               " << kokkos_enable_cuda << std::endl;
    std::cout << "  KOKKOS_ENABLE_CUDA_UVM           " << kokkos_enable_cuda_uvm << std::endl;
    std::cout << "  kokkos_host_execution_space      " << typeid(nimble_kokkos::kokkos_host_execution_space).name() << std::endl;
    std::cout << "  kokkos_host_mirror_memory_space  " << typeid(nimble_kokkos::kokkos_host_mirror_memory_space).name() << std::endl;
    std::cout << "  kokkos_host                      " << typeid(nimble_kokkos::kokkos_host).name() << std::endl;
    std::cout << "  kokkos_device_execution_space    " << typeid(nimble_kokkos::kokkos_device_execution_space).name() << std::endl;
    std::cout << "  kokkos_device_memory_space       " << typeid(nimble_kokkos::kokkos_device_memory_space).name() << std::endl;
    std::cout << "  kokkos_device                    " << typeid(nimble_kokkos::kokkos_device).name() << std::endl;
    std::cout << "  kokkos_layout                    " << typeid(nimble_kokkos::kokkos_layout).name() << std::endl;
    std::cout << std::endl;
  }

  std::string input_deck_name = argv[1];
  nimble::Parser parser;
  parser.Initialize(input_deck_name);

  // Read the mesh
  std::string genesis_file_name = nimble::IOFileName(parser.GenesisFileName(), "g", "", my_mpi_rank, num_mpi_ranks);
  std::string rve_genesis_file_name = nimble::IOFileName(parser.RVEGenesisFileName(), "g");
  nimble::GenesisMesh mesh;
  mesh.ReadFile(genesis_file_name);
  nimble::GenesisMesh rve_mesh;
  if (rve_genesis_file_name != "none") {
    rve_mesh.ReadFile(rve_genesis_file_name);
  }
  std::string tag = "kokkos";
  std::string output_exodus_name = nimble::IOFileName(parser.ExodusFileName(), "e", tag, my_mpi_rank, num_mpi_ranks);
  int dim = mesh.GetDim();
  int num_nodes = mesh.GetNumNodes();
  int num_blocks = mesh.GetNumBlocks();

  nimble_kokkos::DataManager data_manager;
  nimble_kokkos::ModelData & model_data = data_manager.GetMacroScaleData();

  int lumped_mass_field_id = model_data.AllocateNodeData(nimble::SCALAR, "lumped_mass", num_nodes);
  nimble_kokkos::HostScalarNodeView lumped_mass_h = model_data.GetHostScalarNodeData(lumped_mass_field_id);
  nimble_kokkos::DeviceScalarNodeView lumped_mass_d = model_data.GetDeviceScalarNodeData(lumped_mass_field_id);

  int reference_coordinate_field_id = model_data.AllocateNodeData(nimble::VECTOR, "reference_coordinate", num_nodes);
  nimble_kokkos::HostVectorNodeView reference_coordinate_h = model_data.GetHostVectorNodeData(reference_coordinate_field_id);
  nimble_kokkos::DeviceVectorNodeView reference_coordinate_d = model_data.GetDeviceVectorNodeData(reference_coordinate_field_id);

  int displacement_field_id = model_data.AllocateNodeData(nimble::VECTOR, "displacement", num_nodes);
  nimble_kokkos::HostVectorNodeView displacement_h = model_data.GetHostVectorNodeData(displacement_field_id);
  nimble_kokkos::DeviceVectorNodeView displacement_d = model_data.GetDeviceVectorNodeData(displacement_field_id);
  Kokkos::deep_copy(displacement_h, (double)(0.0));

  int velocity_field_id = model_data.AllocateNodeData(nimble::VECTOR, "velocity", num_nodes);
  nimble_kokkos::HostVectorNodeView velocity_h = model_data.GetHostVectorNodeData(velocity_field_id);
  nimble_kokkos::DeviceVectorNodeView velocity_d = model_data.GetDeviceVectorNodeData(velocity_field_id);
  Kokkos::deep_copy(velocity_h, (double)(0.0));

  int acceleration_field_id =  model_data.AllocateNodeData(nimble::VECTOR, "acceleration", num_nodes);
  nimble_kokkos::HostVectorNodeView acceleration_h = model_data.GetHostVectorNodeData(acceleration_field_id);
  Kokkos::deep_copy(acceleration_h, (double)(0.0));

  int internal_force_field_id =  model_data.AllocateNodeData(nimble::VECTOR, "internal_force", num_nodes);
  nimble_kokkos::HostVectorNodeView internal_force_h = model_data.GetHostVectorNodeData(internal_force_field_id);
  nimble_kokkos::DeviceVectorNodeView internal_force_d = model_data.GetDeviceVectorNodeData(internal_force_field_id);

  int contact_force_field_id =  model_data.AllocateNodeData(nimble::VECTOR, "contact_force", num_nodes);
  nimble_kokkos::HostVectorNodeView contact_force_h = model_data.GetHostVectorNodeData(contact_force_field_id);
  nimble_kokkos::DeviceVectorNodeView contact_force_d = model_data.GetDeviceVectorNodeData(contact_force_field_id);

  int deformation_gradient_field_id(-1);
  int stress_field_id(-1);

  // Blocks
  // std::map<int, nimble::Block>& blocks = model_data.GetBlocks();
  std::map<int, nimble_kokkos::Block> blocks;
  std::map<int, nimble_kokkos::Block>::iterator block_it;
  std::vector<int> block_ids = mesh.GetBlockIds();
  for (int i=0 ; i<num_blocks ; i++){
    int block_id = block_ids.at(i);
    std::string const & macro_material_parameters = parser.GetMacroscaleMaterialParameters(block_id);
    std::map<int, std::string> const & rve_material_parameters = parser.GetMicroscaleMaterialParameters();
    std::string rve_bc_strategy = parser.GetMicroscaleBoundaryConditionStrategy();
    int num_elements_in_block = mesh.GetNumElementsInBlock(block_id);
    blocks[block_id] = nimble_kokkos::Block();
    //blocks[block_id].Initialize(macro_material_parameters, rve_material_parameters, rve_mesh, rve_bc_strategy);
    blocks.at(block_id).Initialize(macro_material_parameters, num_elements_in_block);
    //int num_integration_points_per_element = blocks.at(block_id).GetHostElement()->NumIntegrationPointsPerElement();


    std::vector<double> initial_value(9, 0.0);
    initial_value[0] = initial_value[1] = initial_value[2] = 1.0;
    deformation_gradient_field_id = model_data.AllocateIntegrationPointData(block_id,
                                                                            nimble::FULL_TENSOR,
                                                                            "deformation_gradient",
                                                                            num_elements_in_block,
                                                                            initial_value);
    // volume-averaged quantities for I/O are stored as element data
    model_data.AllocateElementData(block_id,
                                   nimble::FULL_TENSOR,
                                   "deformation_gradient",
                                   num_elements_in_block);

    stress_field_id = model_data.AllocateIntegrationPointData(block_id,
                                                              nimble::SYMMETRIC_TENSOR,
                                                              "stress",
                                                              num_elements_in_block);

    // volume-averaged quantities for I/O are stored as element data
    model_data.AllocateElementData(block_id,
                                   nimble::SYMMETRIC_TENSOR,
                                   "stress",
                                   num_elements_in_block);

    if (parser.GetOutputFieldString().find("volume") != std::string::npos) {
      model_data.AllocateElementData(block_id,
                                     nimble::SCALAR,
                                     "volume",
                                     num_elements_in_block);
    }
  }

  // Initialize the initial- and boundary-condition manager
  std::map<int, std::string> const & node_set_names = mesh.GetNodeSetNames();
  std::map<int, std::vector<int> > const & node_sets = mesh.GetNodeSets();
  std::vector<std::string> const & bc_strings = parser.GetBoundaryConditionStrings();
  std::string const & time_integration_scheme = parser.TimeIntegrationScheme();

  nimble::BoundaryConditionManager boundary_condition_manager;
  boundary_condition_manager.Initialize(node_set_names, node_sets, bc_strings, dim, time_integration_scheme);

  // Initialize the output file
  nimble_kokkos::ExodusOutputManager exodus_output_manager;
  exodus_output_manager.SpecifyOutputFields(model_data, parser.GetOutputFieldString());
  std::vector<std::string> global_data_labels;
  std::vector<std::string> node_data_labels_for_output = exodus_output_manager.GetNodeDataLabelsForOutput();
  std::map<int, std::vector<std::string> > elem_data_labels_for_output = exodus_output_manager.GetElementDataLabelsForOutput();
  std::map<int, std::vector<std::string> > derived_elem_data_labels;
  for (int block_id : block_ids) {
    derived_elem_data_labels[block_id] = std::vector<std::string>(); // TODO elliminate this
  }
  // ****
  nimble::ExodusOutput exodus_output;
  exodus_output.Initialize(output_exodus_name, mesh);
  exodus_output.InitializeDatabase(mesh,
                                   global_data_labels,
                                   node_data_labels_for_output,
                                   elem_data_labels_for_output,
                                   derived_elem_data_labels);

  const double * const ref_coord_x = mesh.GetCoordinatesX();
  const double * const ref_coord_y = mesh.GetCoordinatesY();
  const double * const ref_coord_z = mesh.GetCoordinatesZ();
  for (int i=0 ; i<num_nodes ; i++) {
    reference_coordinate_h(i, 0) = ref_coord_x[i];
    reference_coordinate_h(i, 1) = ref_coord_y[i];
    reference_coordinate_h(i, 2) = ref_coord_z[i];
  }
  Kokkos::deep_copy(reference_coordinate_d, reference_coordinate_h);

  int block_index;

  // Containers for gathered data
  std::vector<nimble_kokkos::DeviceScalarNodeGatheredView> gathered_lumped_mass_d(num_blocks, nimble_kokkos::DeviceScalarNodeGatheredView("gathred_lumped_mass", 1));
  std::vector<nimble_kokkos::DeviceVectorNodeGatheredView> gathered_reference_coordinate_d(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_reference_coordinates", 1));
  std::vector<nimble_kokkos::DeviceVectorNodeGatheredView> gathered_displacement_d(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_displacement", 1));
  std::vector<nimble_kokkos::DeviceVectorNodeGatheredView> gathered_internal_force_d(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_internal_force", 1));
  std::vector<nimble_kokkos::DeviceVectorNodeGatheredView> gathered_contact_force_d(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_contact_force", 1));
  for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index+=1, block_it++) {
    int block_id = block_it->first;
    int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    Kokkos::resize(gathered_lumped_mass_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_reference_coordinate_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_displacement_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_internal_force_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_contact_force_d.at(block_index), num_elem_in_block);
  }

  Kokkos::deep_copy(lumped_mass_d, (double)(0.0));

  // Compute the lumped mass
  for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
    int block_id = block_it->first;
    nimble_kokkos::Block& block = block_it->second;
    nimble::Element* element_d = block.GetDeviceElement();
    double density = block.GetDensity();
    int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    int num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);
    int elem_conn_length = num_elem_in_block * num_nodes_per_elem;
    int const * elem_conn = mesh.GetConnectivity(block_id);

    nimble_kokkos::HostElementConnectivityView elem_conn_h("element_connectivity_h", elem_conn_length);
    for (int i=0 ; i<elem_conn_length ; i++) {
      elem_conn_h(i) = elem_conn[i];
    }
    nimble_kokkos::DeviceElementConnectivityView& elem_conn_d = block.GetDeviceElementConnectivityView();
    Kokkos::resize(elem_conn_d, elem_conn_length);
    Kokkos::deep_copy(elem_conn_d, elem_conn_h);

    nimble_kokkos::DeviceVectorNodeGatheredView gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
    nimble_kokkos::DeviceScalarNodeGatheredView gathered_lumped_mass_block_d = gathered_lumped_mass_d.at(block_index);

    model_data.GatherVectorNodeData(reference_coordinate_field_id,
                                    num_elem_in_block,
                                    num_nodes_per_elem,
                                    elem_conn_d,
                                    gathered_reference_coordinate_block_d);

    // COMPUTE LUMPED MASS
    Kokkos::parallel_for("Lumped Mass", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
        nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
        nimble_kokkos::DeviceScalarNodeGatheredSubView element_lumped_mass_d = Kokkos::subview(gathered_lumped_mass_block_d, i_elem, Kokkos::ALL);
        element_d->ComputeLumpedMass(density, element_reference_coordinate_d, element_lumped_mass_d);
      });

    // SCATTER TO NODE DATA
    model_data.ScatterScalarNodeData(lumped_mass_field_id,
                                     num_elem_in_block,
                                     num_nodes_per_elem,
                                     elem_conn_d,
                                     gathered_lumped_mass_block_d);
  }
  Kokkos::deep_copy(lumped_mass_h, lumped_mass_d);

  // Initialize the MPI container
  std::vector<int> global_node_ids(num_nodes);
  int const * const global_node_ids_ptr = mesh.GetNodeGlobalIds();
  for (int n=0 ; n<num_nodes ; ++n) {
    global_node_ids[n] = global_node_ids_ptr[n];
  }
  nimble::MPIContainer mpi_container;
  mpi_container.Initialize(global_node_ids);

  int mpi_scalar_dimension = 1;
  std::vector<double> mpi_scalar_buffer(mpi_scalar_dimension * num_nodes);
  int mpi_vector_dimension = 3;

  nimble::ContactManager contact_manager;
  bool contact_enabled = parser.HasContact();
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
                                          mpi_container,
                                          contact_master_block_ids,
                                          contact_slave_block_ids);
  }

  // MPI vector reduction on lumped mass
  for (unsigned int i=0 ; i<num_nodes ; i++) {
    mpi_scalar_buffer[i] = lumped_mass_h(i);
  }
  mpi_container.VectorReduction(mpi_scalar_dimension, mpi_scalar_buffer.data());
  for (int i=0 ; i<num_nodes ; i++) {
    lumped_mass_h(i) = mpi_scalar_buffer[i];
  }

  double time_current(0.0), time_previous(0.0);
  double final_time = parser.FinalTime();
  double delta_time, half_delta_time;
  int num_load_steps = parser.NumLoadSteps();
  int output_frequency = parser.OutputFrequency();

  // Apply the initial conditions and kinematic boundary conditions
  boundary_condition_manager.ApplyInitialConditions(reference_coordinate_h, velocity_h);
  boundary_condition_manager.ApplyKinematicBC(time_current, time_previous, reference_coordinate_h, displacement_h, velocity_h);
  Kokkos::deep_copy(displacement_d, displacement_h);
  Kokkos::deep_copy(velocity_d, velocity_h);

  // Output to Exodus file
  exodus_output_manager.ComputeElementData(mesh, model_data, blocks, gathered_reference_coordinate_d, gathered_displacement_d);
  std::vector<double> global_data;
  std::vector< std::vector<double> > const & node_data_for_output = exodus_output_manager.GetNodeDataForOutput(model_data);
  std::map<int, std::vector< std::vector<double> > > const & elem_data_for_output = exodus_output_manager.GetElementDataForOutput(model_data);
  std::map<int, std::vector< std::vector<double> > > derived_elem_data;
  exodus_output.WriteStep(time_current,
                          global_data,
                          node_data_for_output,
                          elem_data_labels_for_output,
                          elem_data_for_output,
                          derived_elem_data_labels,
                          derived_elem_data);

  for (int step=0 ; step<num_load_steps ; step++) {

    if (my_mpi_rank == 0) {
      if (10*(step+1) % num_load_steps == 0 && step != num_load_steps - 1) {
        std::cout << "   " << static_cast<int>( 100.0 * static_cast<double>(step+1)/num_load_steps ) << "% complete" << std::endl << std::flush;
      }
      else if (step == num_load_steps - 1) {
        std::cout << "  100% complete\n" << std::endl << std::flush;
      }
    }
    bool is_output_step = false;
    if (step%output_frequency == 0 || step == num_load_steps - 1) {
      is_output_step = true;
    }

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

    // Apply kinematic boundary conditions
    boundary_condition_manager.ApplyKinematicBC(time_current, time_previous, reference_coordinate_h, displacement_h, velocity_h);

    // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
    for (int i=0 ; i<num_nodes ; ++i) {
      displacement_h(i,0) += delta_time * velocity_h(i,0);
      displacement_h(i,1) += delta_time * velocity_h(i,1);
      displacement_h(i,2) += delta_time * velocity_h(i,2);
    }

    // Copy the current displacement and velocity value to device memory
    Kokkos::deep_copy(displacement_d, displacement_h);
    Kokkos::deep_copy(velocity_d, velocity_h);
    Kokkos::deep_copy(internal_force_d, (double)(0.0));
    Kokkos::deep_copy(contact_force_d, (double)(0.0));

    // Compute element-level kinematics
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

      model_data.GatherVectorNodeData(reference_coordinate_field_id, /* TODO SHOULD JUST PASS IN VIEW? */
                                      num_elem_in_block, /* TODO SHOULD BE ABLE TO GET THIS OFF VIEW "EXTENT" */
                                      num_nodes_per_elem,
                                      elem_conn_d,
                                      gathered_reference_coordinate_block_d);

      model_data.GatherVectorNodeData(displacement_field_id,
                                      num_elem_in_block,
                                      num_nodes_per_elem,
                                      elem_conn_d,
                                      gathered_displacement_block_d);

      nimble_kokkos::DeviceFullTensorIntPtView deformation_gradient_step_np1_d = model_data.GetDeviceFullTensorIntegrationPointData(block_id,
                                                                                                                                    deformation_gradient_field_id,
                                                                                                                                    nimble::STEP_NP1);
      // COMPUTE DEFORMATION GRADIENTS
      Kokkos::parallel_for("Deformation Gradient", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
          nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d = Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceFullTensorIntPtSubView element_deformation_gradient_step_np1_d = Kokkos::subview(deformation_gradient_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          element_d->ComputeDeformationGradients(element_reference_coordinate_d,
                                                 element_displacement_d,
                                                 element_deformation_gradient_step_np1_d);
        });
    }

    // Compute stress
    for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
      int block_id = block_it->first;
      nimble_kokkos::Block& block = block_it->second;
      nimble::Material* material_d = block.GetDeviceMaterialModel();
      int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      int num_integration_points_per_element = block.GetHostElement()->NumIntegrationPointsPerElement();

      nimble_kokkos::DeviceFullTensorIntPtView deformation_gradient_step_np1_d = model_data.GetDeviceFullTensorIntegrationPointData(block_id,
                                                                                                                                    deformation_gradient_field_id,
                                                                                                                                    nimble::STEP_NP1);

      nimble_kokkos::DeviceFullTensorIntPtView deformation_gradient_step_n_d = model_data.GetDeviceFullTensorIntegrationPointData(block_id,
                                                                                                                                  deformation_gradient_field_id,
                                                                                                                                  nimble::STEP_N);

      nimble_kokkos::DeviceSymTensorIntPtView stress_step_n_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id,
                                                                                                                  stress_field_id,
                                                                                                                  nimble::STEP_N);

      nimble_kokkos::DeviceSymTensorIntPtView stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id,
                                                                                                                    stress_field_id,
                                                                                                                    nimble::STEP_NP1);

      bool is_ngp_lame_model = block.GetHostMaterialModel()->IsNGPLAMEModel();

      if (!is_ngp_lame_model) {
        typedef typename Kokkos::MDRangePolicy< Kokkos::Rank<2> > MDPolicyType_2D;
        MDPolicyType_2D mdpolicy_2d( {{0,0}}, {{num_elem_in_block,num_integration_points_per_element}} );
        Kokkos::parallel_for("Stress", mdpolicy_2d, KOKKOS_LAMBDA (const int i_elem, const int i_ipt) {
            nimble_kokkos::DeviceFullTensorIntPtSingleEntryView element_deformation_gradient_step_n_d = Kokkos::subview(deformation_gradient_step_n_d, i_elem, i_ipt, Kokkos::ALL);
            nimble_kokkos::DeviceFullTensorIntPtSingleEntryView element_deformation_gradient_step_np1_d = Kokkos::subview(deformation_gradient_step_np1_d, i_elem, i_ipt, Kokkos::ALL);
            nimble_kokkos::DeviceSymTensorIntPtSingleEntryView element_stress_step_n_d = Kokkos::subview(stress_step_n_d, i_elem, i_ipt, Kokkos::ALL);
            nimble_kokkos::DeviceSymTensorIntPtSingleEntryView element_stress_step_np1_d = Kokkos::subview(stress_step_np1_d, i_elem, i_ipt, Kokkos::ALL);
            material_d->GetStress(time_previous,
                                  time_current,
                                  element_deformation_gradient_step_n_d,
                                  element_deformation_gradient_step_np1_d,
                                  element_stress_step_n_d,
                                  element_stress_step_np1_d);
          });
      }
      else {
#ifdef NIMBLE_HAVE_EXTRAS
        nimble::NGPLAMEMaterial::NGPLAMEData ngp_lame_data = *(block.GetNGPLAMEData());
        nimble_kokkos::DeviceSymTensorIntPtView stress_step_n_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id,
                                                                                                                    stress_field_id,
                                                                                                                    nimble::STEP_N);
        typedef typename Kokkos::MDRangePolicy< Kokkos::Rank<2> > MDPolicyType_2D;
        MDPolicyType_2D mdpolicy_2d( {{0,0}}, {{num_elem_in_block,num_integration_points_per_element}} );
        Kokkos::parallel_for("NGP LAME Kinematics", mdpolicy_2d, KOKKOS_LAMBDA (const int i_elem, const int i_ipt) {
            nimble_kokkos::DeviceFullTensorIntPtSingleEntryView element_deformation_gradient_step_n_d = Kokkos::subview(deformation_gradient_step_n_d, i_elem, i_ipt, Kokkos::ALL);
            nimble_kokkos::DeviceFullTensorIntPtSingleEntryView element_deformation_gradient_step_np1_d = Kokkos::subview(deformation_gradient_step_np1_d, i_elem, i_ipt, Kokkos::ALL);
            nimble_kokkos::DeviceSymTensorIntPtSingleEntryView element_stress_step_n_d = Kokkos::subview(stress_step_n_d, i_elem, i_ipt, Kokkos::ALL);
            double disp_grad_step_n[9], disp_grad_step_np1[9];
            for (int i=0 ; i<9 ; i++) {
              disp_grad_step_n[i] = element_deformation_gradient_step_n_d(i);
              disp_grad_step_np1[i] = element_deformation_gradient_step_np1_d(i);
              if (i<3) {
                disp_grad_step_n[i] -= 1.0;
                disp_grad_step_np1[i] -= 1.0;
              }
            }
            int i_mat_pt = i_elem*num_integration_points_per_element + i_ipt;
            for (int i=0 ; i<9 ; i++) {
              ngp_lame_data.disp_grad_(i_mat_pt, i) = disp_grad_step_np1[i];
              ngp_lame_data.velo_grad_(i_mat_pt, i) = (disp_grad_step_np1[i] - disp_grad_step_n[i]) / delta_time;
            }
            for (int i=0 ; i<6 ; i++) {
              ngp_lame_data.stress_old_(i_mat_pt, i) = element_stress_step_n_d(i);
            }
            // TODO handle state variables
          });
#endif
      }
    }

#ifdef NIMBLE_HAVE_EXTRAS
    if (blocks.begin()->second.GetHostMaterialModel()->IsNGPLAMEModel()) {

      // Create a view of NGP LAME material models for use in team loop
      Kokkos::View<mtk_ngp::DevicePtr<lame::ngp::Material>*> ngp_lame_materials_d("NGP LAME Material Models", num_blocks);
      auto ngp_lame_materials_h = Kokkos::create_mirror_view(ngp_lame_materials_d);
      for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
        nimble::NGPLAMEMaterial* nimble_ngp_lame_material_model = dynamic_cast<nimble::NGPLAMEMaterial*>(block_it->second.GetHostMaterialModel().get());
        ngp_lame_materials_h(block_index) = nimble_ngp_lame_material_model->GetNGPLAMEMaterialModel();
      }
      Kokkos::deep_copy(ngp_lame_materials_d, ngp_lame_materials_h);

      // Compute stress
      const Kokkos::TeamPolicy<> team_loop(num_blocks, Kokkos::AUTO);
      Kokkos::parallel_for("NGP LAME Stress", team_loop, KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
          int team_id = team.league_rank();
          ngp_lame_materials_d[team_id]->compute_stress(team);
        });

      // Copy stress back to the Nimble data structures
      for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
        int block_id = block_it->first;
        nimble_kokkos::Block& block = block_it->second;
        int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
        int num_integration_points_per_element = block.GetHostElement()->NumIntegrationPointsPerElement();
        nimble::NGPLAMEMaterial::NGPLAMEData ngp_lame_data = *(block.GetNGPLAMEData());
        nimble_kokkos::DeviceSymTensorIntPtView stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id,
                                                                                                                      stress_field_id,
                                                                                                                      nimble::STEP_NP1);
        typedef typename Kokkos::MDRangePolicy< Kokkos::Rank<2> > MDPolicyType_2D;
        MDPolicyType_2D mdpolicy_2d( {{0,0}}, {{num_elem_in_block,num_integration_points_per_element}} );
        Kokkos::parallel_for("NGP LAME Stress Copy", mdpolicy_2d, KOKKOS_LAMBDA (const int i_elem, const int i_ipt) {
            nimble_kokkos::DeviceSymTensorIntPtSingleEntryView element_stress_step_np1_d = Kokkos::subview(stress_step_np1_d, i_elem, i_ipt, Kokkos::ALL);
            int i_mat_pt = i_elem*num_integration_points_per_element + i_ipt;
            for (int i=0 ; i<6 ; i++) {
              element_stress_step_np1_d(i) = ngp_lame_data.stress_new_(i_mat_pt, i);
            }
            // TODO handle state variables
        });
      }
    }
#endif

    // Stress divergence
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
                                                                                                                    stress_field_id,
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

      model_data.ScatterVectorNodeData(internal_force_field_id,
                                       num_elem_in_block,
                                       num_nodes_per_elem,
                                       elem_conn_d,
                                       gathered_internal_force_block_d);
    } // loop over blocks

    Kokkos::deep_copy(internal_force_h, internal_force_d);

    // Perform a reduction to obtain correct values on MPI boundaries
    mpi_container.VectorReduction(mpi_vector_dimension, internal_force_h);


    // Evaluate the contact force
    if (contact_enabled) {
      contact_manager.ApplyDisplacements(displacement_d);

      double x_min, x_max, y_min, y_max, z_min, z_max;
      contact_manager.BoundingBox(x_min, x_max, y_min, y_max, z_min, z_max);
      // if (is_output_step) {
      //   contact_manager.WriteContactBoundingBoxToVTKFile("contact_bounding_box_", my_mpi_rank, step);
      // }

      contact_manager.ComputeContactForce(step+1, false);
      contact_manager.GetForces(contact_force_d);
    }
    Kokkos::deep_copy(contact_force_h, contact_force_d);

    // fill acceleration vector A^{n+1} = M^{-1} ( F^{n} + b^{n} )
    for (int i=0 ; i<num_nodes ; ++i) {
      acceleration_h(i,0) = (1.0/lumped_mass_h(i)) * (internal_force_h(i,0) + contact_force_h(i,0));
      acceleration_h(i,1) = (1.0/lumped_mass_h(i)) * (internal_force_h(i,1) + contact_force_h(i,1));
      acceleration_h(i,2) = (1.0/lumped_mass_h(i)) * (internal_force_h(i,2) + contact_force_h(i,2));
    }

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    for (int i=0 ; i<num_nodes ; ++i) {
      velocity_h(i,0) += half_delta_time * acceleration_h(i,0);
      velocity_h(i,1) += half_delta_time * acceleration_h(i,1);
      velocity_h(i,2) += half_delta_time * acceleration_h(i,2);
    }

    if (is_output_step) {

      boundary_condition_manager.ApplyKinematicBC(time_current, time_previous, reference_coordinate_h, displacement_h, velocity_h);
      Kokkos::deep_copy(displacement_d, displacement_h);
      Kokkos::deep_copy(velocity_d, velocity_h);
      exodus_output_manager.ComputeElementData(mesh, model_data, blocks, gathered_reference_coordinate_d, gathered_displacement_d);
      std::vector<double> global_data;
      std::vector< std::vector<double> > const & node_data_for_output = exodus_output_manager.GetNodeDataForOutput(model_data);
      std::map<int, std::vector< std::vector<double> > > const & elem_data_for_output = exodus_output_manager.GetElementDataForOutput(model_data);
      std::map<int, std::vector< std::vector<double> > > derived_elem_data;
      exodus_output.WriteStep(time_current,
                              global_data,
                              node_data_for_output,
                              elem_data_labels_for_output,
                              elem_data_for_output,
                              derived_elem_data_labels,
                              derived_elem_data);
    }
  }

  if (my_mpi_rank == 0) {
    std::cout << "\ncomplete.\n" << std::endl;
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  //Kokkos::print_configuration(std::cout);
  main_routine(argc, argv);
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
