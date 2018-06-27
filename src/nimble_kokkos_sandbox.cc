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
// #include "nimble_genesis_mesh.h"
// #include "nimble_exodus_output.h"
#include "nimble_exodus_output_manager.h"
#include "nimble_boundary_condition_manager.h"
#include "nimble_kokkos_defs.h"
#include "nimble_kokkos_data_manager.h"
// #include "nimble_linear_solver.h"
// #include "nimble_utils.h"
// #include "nimble_mesh_utils.h"
#include "nimble.mpi.utils.h"
// #include "nimble.quanta.stopwatch.h"
#include "nimble_view.h"

#include "nimble_kokkos_block.h"

#include <iostream>
// #include <iomanip>
// #include <cstdlib>
// #include <random>
// #include <limits>
// #include <fstream>
// #include <sstream>

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
    std::cout << "\n--KokkosSandbox" << std::endl;
    std::cout << "-- version " << nimble::NimbleVersion() << "\n" << std::endl;
    if (argc != 2) {
      std::cout << "Usage:  mpirun -np NP KokkosSandbox <input_deck.in>\n" << std::endl;
      Kokkos::finalize();
      MPI_Finalize();
      exit(1);
    }
    std::cout << "KokkosSandbox initialized on " << num_mpi_ranks << " mpi rank(s)." << std::endl;
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
  std::string tag = "mpi";
  std::string output_exodus_name = nimble::IOFileName(parser.ExodusFileName(), "e", tag, my_mpi_rank, num_mpi_ranks);
  int dim = mesh.GetDim();
  int num_nodes = mesh.GetNumNodes();
  int num_blocks = mesh.GetNumBlocks();

  nimble_kokkos::DataManager data_manager;
  nimble_kokkos::ModelData & model_data = data_manager.GetMacroScaleData();

  int lumped_mass_field_id = model_data.AllocateNodeData(nimble::SCALAR, "lumped_mass", num_nodes);
  nimble_kokkos::HostScalarView lumped_mass_h = model_data.GetHostScalarNodeData(lumped_mass_field_id);
  nimble_kokkos::DeviceScalarView lumped_mass_d = model_data.GetDeviceScalarNodeData(lumped_mass_field_id);

  int reference_coordinate_field_id = model_data.AllocateNodeData(nimble::VECTOR, "reference_coordinate", num_nodes);
  nimble_kokkos::HostVectorView reference_coordinate_h = model_data.GetHostVectorNodeData(reference_coordinate_field_id);
  nimble_kokkos::DeviceVectorView reference_coordinate_d = model_data.GetDeviceVectorNodeData(reference_coordinate_field_id);

  int displacement_field_id = model_data.AllocateNodeData(nimble::VECTOR, "displacement", num_nodes);
  nimble_kokkos::HostVectorView displacement_h = model_data.GetHostVectorNodeData(displacement_field_id);
  nimble_kokkos::DeviceVectorView displacement_d = model_data.GetDeviceVectorNodeData(displacement_field_id);
  Kokkos::deep_copy(displacement_h, (double)(0.0));

  int velocity_field_id = model_data.AllocateNodeData(nimble::VECTOR, "velocity", num_nodes);
  nimble_kokkos::HostVectorView velocity_h = model_data.GetHostVectorNodeData(velocity_field_id);
  nimble_kokkos::DeviceVectorView velocity_d = model_data.GetDeviceVectorNodeData(velocity_field_id);
  Kokkos::deep_copy(velocity_h, (double)(0.0));

  int acceleration_field_id =  model_data.AllocateNodeData(nimble::VECTOR, "acceleration", num_nodes);
  nimble_kokkos::HostVectorView acceleration_h = model_data.GetHostVectorNodeData(acceleration_field_id);
  Kokkos::deep_copy(acceleration_h, (double)(0.0));

  int internal_force_field_id =  model_data.AllocateNodeData(nimble::VECTOR, "internal_force", num_nodes);
  nimble_kokkos::HostVectorView internal_force_h = model_data.GetHostVectorNodeData(internal_force_field_id);
  nimble_kokkos::DeviceVectorView internal_force_d = model_data.GetDeviceVectorNodeData(internal_force_field_id);

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
    blocks[block_id] = nimble_kokkos::Block();
    //blocks[block_id].Initialize(macro_material_parameters, rve_material_parameters, rve_mesh, rve_bc_strategy);
    blocks.at(block_id).Initialize(macro_material_parameters);
    //int num_integration_points_per_element = blocks.at(block_id).GetHostElement()->NumIntegrationPointsPerElement();
    int num_elements_in_block = mesh.GetNumElementsInBlock(block_id);

    deformation_gradient_field_id =  model_data.AllocateIntegrationPointData(block_id,
                                                                             nimble::FULL_TENSOR,
                                                                             "deformation_gradient",
                                                                             num_elements_in_block);

    stress_field_id = model_data.AllocateIntegrationPointData(block_id,
                                                              nimble::SYMMETRIC_TENSOR,
                                                              "stress",
                                                              num_elements_in_block);

    // NEED TO MANAGE ALLOCATION OF DEFORMATION GRADIENT AND STRESS
    //   PROBABLY DON'T WANT THE MATERIALS TO DECLARE DEFORMATION GRADIENT AND STRESS AS STATE DATA, RATHER JUST SET THEM HERE
    //   NEED TO STORE FIELD IDS FOR DEFORMATION GRADIENT AND STRESS
    //   NEED TO THINK ABOUT STATE DATA

    // FOR STATE DATA (BUT CURRENTLY INCLUDE KINEMATICS & STRESS)
    // std::vector< std::pair<std::string, nimble::Length> > data_labels_and_lengths;
    // blocks.at(block_id).GetIntegrationPointDataLabelsAndLengths(data_labels_and_lengths);
    // for (auto & entry: data_labels_and_lengths) {
    //   std::cout << "DEBUGGING integration point data " << entry.first << std::endl;
    //   model_data.AllocateIntegrationPointData(block_id,
    //                                           entry.second,
    //                                           entry.first,
    //                                           num_integration_points);
    // }
  }
  // std::map<int, int> num_elem_in_each_block = mesh.GetNumElementsInBlock();
  // model_data.AllocateElementData(num_elem_in_each_block);
  // model_data.SpecifyOutputFields(parser.GetOutputFieldString());
  // std::map<int, std::vector<std::string> > const & elem_data_labels = model_data.GetElementDataLabels();
  // std::map<int, std::vector<std::string> > const & elem_data_labels_for_output = model_data.GetElementDataLabelsForOutput();
  // std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();

  // // Initialize the element data
  // std::vector<int> rve_output_elem_ids = parser.MicroscaleOutputElementIds();
  // for (block_it=blocks.begin(); block_it!=blocks.end() ; block_it++) {
  //   int block_id = block_it->first;
  //   int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
  //   std::vector<int> const & elem_global_ids = mesh.GetElementGlobalIdsInBlock(block_id);
  //   nimble::Block& block = block_it->second;
  //   std::vector<double> & elem_data_n = model_data.GetElementDataOld(block_id);
  //   std::vector<double> & elem_data_np1 = model_data.GetElementDataNew(block_id);
  //   block.InitializeElementData(num_elem_in_block,
  //                               elem_global_ids,
  //                               rve_output_elem_ids,
  //                               elem_data_labels.at(block_id),
  //                               derived_elem_data_labels.at(block_id),
  //                               elem_data_n,
  //                               elem_data_np1,
  //                               data_manager);
  // }

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
    derived_elem_data_labels[block_id] = std::vector<std::string>();
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
  std::vector<nimble_kokkos::DeviceScalarGatheredView> gathered_lumped_mass_d(num_blocks, nimble_kokkos::DeviceScalarGatheredView("gathred_lumped_mass", 1));
  std::vector<nimble_kokkos::DeviceVectorGatheredView> gathered_reference_coordinate_d(num_blocks, nimble_kokkos::DeviceVectorGatheredView("gathered_reference_coordinates", 1));
  std::vector<nimble_kokkos::DeviceVectorGatheredView> gathered_displacement_d(num_blocks, nimble_kokkos::DeviceVectorGatheredView("gathered_displacement", 1));
  std::vector<nimble_kokkos::DeviceVectorGatheredView> gathered_internal_force_d(num_blocks, nimble_kokkos::DeviceVectorGatheredView("gathered_internal_force", 1));
  for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index+=1, block_it++) {
    int block_id = block_it->first;
    int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    Kokkos::resize(gathered_lumped_mass_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_reference_coordinate_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_displacement_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_internal_force_d.at(block_index), num_elem_in_block);
  }

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

    nimble_kokkos::DeviceVectorGatheredView gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
    nimble_kokkos::DeviceScalarGatheredView gathered_lumped_mass_block_d = gathered_lumped_mass_d.at(block_index);

    model_data.GatherVectorNodeData(reference_coordinate_field_id,
                                    num_elem_in_block,
                                    num_nodes_per_elem,
                                    elem_conn_d,
                                    gathered_reference_coordinate_block_d);

    // COMPUTE LUMPED MASS
    Kokkos::parallel_for("Lumped Mass", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
        nimble_kokkos::DeviceVectorGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
        nimble_kokkos::DeviceScalarGatheredSubView element_lumped_mass_d = Kokkos::subview(gathered_lumped_mass_block_d, i_elem, Kokkos::ALL);
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
  mpi_container.Initialize(global_node_ids, parser.ReductionVersion());

  // MPI vector reduction on lumped mass
  int data_dimension = 1;
  std::vector<double> data(data_dimension * num_nodes);
  for (unsigned int i=0 ; i<data.size() ; i++) {
    data[i] = lumped_mass_h(i);
  }
  mpi_container.VectorReduction(data_dimension, data.data());
  for (int i=0 ; i<num_nodes ; i++) {
    lumped_mass_h(i) = data[i];
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

  // Compute derived element data
  //   std::map<int, Block> & blocks = model_data.GetBlocks();
  // std::map<int, std::vector<std::string> > const & elem_data_labels = model_data.GetElementDataLabels();
  // std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
  // int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
  // int displacement_field_id = model_data.GetFieldId("displacement");
  // const double * const reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
  // const double * const displacement = model_data.GetNodeData(displacement_field_id);

  // for (auto & entry : blocks) {
  //   int block_id = entry.first;
  //   Block & block = entry.second;
  //   if (derived_element_data_array[qthread_job_id].find(block_id) == derived_element_data_array[qthread_job_id].end()) {
  //     derived_element_data_array[qthread_job_id][block_id] = std::vector< std::vector<double> >();
  //   }
  //   int num_elem_in_block = mesh_array[qthread_job_id].GetNumElementsInBlock(block_id);
  //   const int * const elem_conn = mesh_array[qthread_job_id].GetConnectivity(block_id);
  //   std::vector<double> const & elem_data_np1 = model_data.GetElementDataNew(block_id);
  //   block.ComputeDerivedElementData(reference_coordinate,
  //                                   displacement,
  //                                   num_elem_in_block,
  //                                   elem_conn,
  //                                   elem_data_labels.at(block_id).size(),
  //                                   elem_data_np1,
  //                                   derived_elem_data_labels.at(block_id).size(),
  //                                   derived_element_data_array[qthread_job_id].at(block_id));
  // }

  // Output to Exodus file

  {
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

    // Compute the internal force
    for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
      int block_id = block_it->first;
      nimble_kokkos::Block& block = block_it->second;
      nimble::Element* element_d = block.GetDeviceElement();
      nimble::Material* material_d = block.GetDeviceMaterialModel();
      int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      int num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);
      int num_integration_points_per_element = block.GetHostElement()->NumIntegrationPointsPerElement();

      nimble_kokkos::DeviceElementConnectivityView elem_conn_d = block.GetDeviceElementConnectivityView();
      nimble_kokkos::DeviceVectorGatheredView gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
      nimble_kokkos::DeviceVectorGatheredView gathered_displacement_block_d = gathered_displacement_d.at(block_index);
      nimble_kokkos::DeviceVectorGatheredView gathered_internal_force_block_d = gathered_internal_force_d.at(block_index);

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

      nimble_kokkos::DeviceFullTensorView deformation_gradient_step_np1_d = model_data.GetDeviceFullTensorIntegrationPointData(block_id,
                                                                                                                               deformation_gradient_field_id,
                                                                                                                               nimble::STEP_NP1);

      // COMPUTE DEFORMATION GRADIENTS
      Kokkos::parallel_for("Deformation Gradient", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
          nimble_kokkos::DeviceVectorGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceVectorGatheredSubView element_displacement_d = Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceFullTensorSubView element_deformation_gradient_step_np1_d = Kokkos::subview(deformation_gradient_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          element_d->ComputeDeformationGradients(element_reference_coordinate_d,
                                                 element_displacement_d,
                                                 element_deformation_gradient_step_np1_d);
          // if(std::abs(element_deformation_gradient_step_np1_d(0,0)-1.0) > 1.0e-10)
          //   printf("\nDEF GRAD %e", element_deformation_gradient_step_np1_d(0,0));
        });

      nimble_kokkos::DeviceFullTensorView deformation_gradient_step_n_d = model_data.GetDeviceFullTensorIntegrationPointData(block_id,
                                                                                                                             deformation_gradient_field_id,
                                                                                                                             nimble::STEP_N);

      nimble_kokkos::DeviceSymTensorView stress_step_n_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id,
                                                                                                             stress_field_id,
                                                                                                             nimble::STEP_N);

      nimble_kokkos::DeviceSymTensorView stress_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id,
                                                                                                               stress_field_id,
                                                                                                               nimble::STEP_NP1);

      // COMPUTE STRESS
      Kokkos::parallel_for("Stress", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {

          // TODO THIS LOOP SHOULD PROBABLY SOMEHOW BE PART OF THE PARALLEL_FOR
          for (int i_int_pt=0 ; i_int_pt<num_integration_points_per_element; i_int_pt++) {
            nimble_kokkos::DeviceFullTensorSingleEntryView element_deformation_gradient_step_n_d = Kokkos::subview(deformation_gradient_step_n_d, i_elem, i_int_pt, Kokkos::ALL);
            nimble_kokkos::DeviceFullTensorSingleEntryView element_deformation_gradient_step_np1_d = Kokkos::subview(deformation_gradient_step_np1_d, i_elem, i_int_pt, Kokkos::ALL);
            nimble_kokkos::DeviceSymTensorSingleEntryView element_stress_step_n_d = Kokkos::subview(stress_step_n_d, i_elem, i_int_pt, Kokkos::ALL);
            nimble_kokkos::DeviceSymTensorSingleEntryView element_stress_step_np1_d = Kokkos::subview(stress_step_np1_d, i_elem, i_int_pt, Kokkos::ALL);

            // TODO RE-ADD CONST FOR FUNCTION ARGS?

            material_d->GetStress(time_previous,
                                  time_current,
                                  element_deformation_gradient_step_n_d,
                                  element_deformation_gradient_step_np1_d,
                                  element_stress_step_n_d,
                                  element_stress_step_np1_d);

          }
        });

      // COMPUTE NODAL FORCES
      Kokkos::parallel_for("Force", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
          nimble_kokkos::DeviceVectorGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceVectorGatheredSubView element_displacement_d = Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceSymTensorSubView element_stress_step_np1_d = Kokkos::subview(stress_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceVectorGatheredSubView element_internal_force_d = Kokkos::subview(gathered_internal_force_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
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
    int data_dimension = 3;
    std::vector<double> data(data_dimension * num_nodes);
    for (int i=0 ; i<num_nodes ; i++) {
      data[3*i  ] = internal_force_h(i,0);
      data[3*i+1] = internal_force_h(i,1);
      data[3*i+2] = internal_force_h(i,2);
    }
    mpi_container.VectorReduction(data_dimension, data.data());
    for (int i=0 ; i<num_nodes ; i++) {
      internal_force_h(i,0) = data[3*i  ];
      internal_force_h(i,1) = data[3*i+1];
      internal_force_h(i,2) = data[3*i+2];
    }

    // fill acceleration vector A^{n+1} = M^{-1} ( F^{n} + b^{n} )
    for (int i=0 ; i<num_nodes ; ++i) {
      acceleration_h(i,0) = (1.0/lumped_mass_h(i)) * (internal_force_h(i,0));
      acceleration_h(i,1) = (1.0/lumped_mass_h(i)) * (internal_force_h(i,1));
      acceleration_h(i,2) = (1.0/lumped_mass_h(i)) * (internal_force_h(i,2));
    }

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    for (int i=0 ; i<num_nodes ; ++i) {
      velocity_h(i,0) += half_delta_time * acceleration_h(i,0);
      velocity_h(i,1) += half_delta_time * acceleration_h(i,1);
      velocity_h(i,2) += half_delta_time * acceleration_h(i,2);
    }

    if (is_output_step) {

      boundary_condition_manager.ApplyKinematicBC(time_current, time_previous, reference_coordinate_h, displacement_h, velocity_h);

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
  main_routine(argc, argv);
  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}
