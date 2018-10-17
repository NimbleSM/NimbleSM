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
#include "nimble_genesis_mesh.h"
#include "nimble_exodus_output.h"
#include "nimble_boundary_condition_manager.h"
#include "nimble_data_manager.h"
#include "nimble_linear_solver.h"
#include "nimble_utils.h"
#include "nimble_mesh_utils.h"
#include "nimble.mpi.utils.h"
#include "nimble.quanta.stopwatch.h"
#include "nimble_view.h"

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <limits>
#include <fstream>
#include <sstream>
#include <qthread/qthread.hpp>

double time_previous_g, time_current_g;
bool is_output_step_g;
nimble::Parser parser_g;
nimble::GenesisMesh *mesh_array;
nimble::GenesisMesh *rve_mesh_array;
nimble::DataManager *data_manager_array;
nimble::BoundaryConditionManager *boundary_condition_manager_array;
nimble::ExodusOutput *exodus_output_array;
std::map<int, std::vector<std::vector<double> > > *derived_element_data_array;

typedef struct qthread_arg {
  int job_id;
  int mpi_rank;
  int mpi_num_ranks;
} qthread_arg_t;

static aligned_t QthreadsReadGenesisFiles(void *arg)
{
  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int mpi_rank = arg_struct->mpi_rank;
  int mpi_num_ranks = arg_struct->mpi_num_ranks;
  int qthread_job_id = arg_struct->job_id;
  int qthread_number_of_workers = qthread_num_workers();

  int name_mpi_rank = mpi_rank;
  int name_mpi_num_ranks = mpi_num_ranks;
  int name_qthread_job_id = qthread_job_id;
  int name_qthread_number_of_workers = qthread_number_of_workers;

  if (parser_g.UseTwoLevelMeshDecomposition() == false) {
    name_mpi_rank = mpi_rank * qthread_number_of_workers + qthread_job_id;
    name_mpi_num_ranks = mpi_num_ranks * qthread_number_of_workers;
    name_qthread_job_id = -1;
    name_qthread_number_of_workers = 0;
  }

  char serial_file_name[nimble::MAX_C_STR_SIZE];
  strcpy(serial_file_name, parser_g.GenesisFileName().c_str());

  char file_name_extension[nimble::MAX_C_STR_SIZE];
  file_name_extension[0] = 'g';
  file_name_extension[1] = '\0';

  char file_label[nimble::MAX_C_STR_SIZE];
  file_label[0] = '\0';

  char genesis_file_name[nimble::MAX_C_STR_SIZE];

  nimble::IOFileNameThreadSafe(serial_file_name,
                               file_name_extension,
                               file_label,
                               name_mpi_rank,
                               name_mpi_num_ranks,
                               name_qthread_job_id,
                               name_qthread_number_of_workers,
                               genesis_file_name);

  mesh_array[qthread_job_id].ReadFile(genesis_file_name);

  std::string parser_rve_file_name = parser_g.RVEGenesisFileName();
  if (parser_rve_file_name != "none") {
    strcpy(serial_file_name, parser_rve_file_name.c_str());
    char rve_genesis_file_name[nimble::MAX_C_STR_SIZE];
    nimble::IOFileNameThreadSafe(serial_file_name,
                                 file_name_extension,
                                 file_label,
                                 0,
                                 0,
                                 0,
                                 0,
                                 rve_genesis_file_name);

    rve_mesh_array[qthread_job_id].ReadFile(rve_genesis_file_name);
  }

  return 0;
}

static aligned_t QthreadsInitializeDataManager(void *arg)
{
  using nimble::ModelData;
  using nimble::Block;
  using nimble::Length;
  using nimble::SCALAR;
  using nimble::VECTOR;

  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int mpi_rank = arg_struct->mpi_rank;
  int mpi_num_ranks = arg_struct->mpi_num_ranks;
  int qthread_job_id = arg_struct->job_id;
  int qthread_number_of_workers = qthread_num_workers();

  int dim = mesh_array[qthread_job_id].GetDim();
  int num_nodes = mesh_array[qthread_job_id].GetNumNodes();
  int num_blocks = mesh_array[qthread_job_id].GetNumBlocks();
  std::vector<int> block_ids = mesh_array[qthread_job_id].GetBlockIds();

  bool rve_is_valid = rve_mesh_array[qthread_job_id].IsValid();

  // Global data
  int num_global_data = 0;
  std::vector<std::string> global_data_labels(num_global_data);
  std::vector<double> global_data(num_global_data);

  ModelData& model_data = data_manager_array[qthread_job_id].GetMacroScaleData();
  model_data.SetDimension(dim);

  int lumped_mass_field_id = model_data.AllocateNodeData(SCALAR, "lumped_mass", num_nodes);
  int reference_coordinate_field_id = model_data.AllocateNodeData(VECTOR, "reference_coordinate", num_nodes);
  int displacement_field_id = model_data.AllocateNodeData(VECTOR, "displacement", num_nodes);
  int velocity_field_id = model_data.AllocateNodeData(VECTOR, "velocity", num_nodes);
  int acceleration_field_id = model_data.AllocateNodeData(VECTOR, "acceleration", num_nodes);
  int internal_force_field_id = model_data.AllocateNodeData(VECTOR, "internal_force", num_nodes);
  int external_force_field_id = model_data.AllocateNodeData(VECTOR, "external_force", num_nodes);

  const double * const ref_coord_x = mesh_array[qthread_job_id].GetCoordinatesX();
  const double * const ref_coord_y = mesh_array[qthread_job_id].GetCoordinatesY();
  const double * const ref_coord_z = mesh_array[qthread_job_id].GetCoordinatesZ();
  double* reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
  for (int i=0 ; i<num_nodes ; i++) {
    reference_coordinate[3*i]   = ref_coord_x[i];
    reference_coordinate[3*i+1] = ref_coord_y[i];
    reference_coordinate[3*i+2] = ref_coord_z[i];
  }

  // Initialize the element blocks
  std::map<int, Block> & blocks = model_data.GetBlocks();
  for (int i=0 ; i<num_blocks ; i++){
    int block_id = block_ids[i];
    int num_elem_in_block = mesh_array[qthread_job_id].GetNumElementsInBlock(block_id);
    std::string const & macro_material_parameters = parser_g.GetMacroscaleMaterialParameters(block_id);
    std::map<int, std::string> const & rve_material_parameters = parser_g.GetMicroscaleMaterialParameters();
    std::string rve_bc_strategy = parser_g.GetMicroscaleBoundaryConditionStrategy();
    blocks[block_id] = Block();
    blocks[block_id].Initialize(macro_material_parameters, rve_material_parameters, rve_mesh_array[qthread_job_id], rve_bc_strategy);
    std::vector< std::pair<std::string, Length> > data_labels_and_lengths;
    blocks[block_id].GetDataLabelsAndLengths(data_labels_and_lengths);
    model_data.DeclareElementData(block_id, data_labels_and_lengths);
  }
  std::map<int, int> num_elem_in_each_block = mesh_array[qthread_job_id].GetNumElementsInBlock();
  model_data.AllocateElementData(num_elem_in_each_block);
  model_data.SpecifyOutputFields(parser_g.GetOutputFieldString());
  std::map<int, std::vector<std::string> > const & elem_data_labels = model_data.GetElementDataLabels();
  std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
  std::vector<int> rve_output_elem_ids = parser_g.MicroscaleOutputElementIds();

  for (auto & entry : blocks) {
    int block_id = entry.first;
    Block& block = entry.second;
    int num_elem_in_block = mesh_array[qthread_job_id].GetNumElementsInBlock(block_id);
    std::vector<int> const & elem_global_ids = mesh_array[qthread_job_id].GetElementGlobalIdsInBlock(block_id);
    std::vector<double> & elem_data_n = model_data.GetElementDataOld(block_id);
    std::vector<double> & elem_data_np1 = model_data.GetElementDataNew(block_id);
    block.InitializeElementData(num_elem_in_block,
                                elem_global_ids,
                                rve_output_elem_ids,
                                elem_data_labels.at(block_id),
                                derived_elem_data_labels.at(block_id),
                                elem_data_n,
                                elem_data_np1,
                                data_manager_array[qthread_job_id]);
  }
}

static aligned_t QthreadsInitializeBoundaryConditionManager(void *arg)
{
  using nimble::GenesisMesh;

  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int qthread_job_id = arg_struct->job_id;

  GenesisMesh const & macroscale_mesh = mesh_array[qthread_job_id];

  std::map<int, std::string> const & node_set_names = macroscale_mesh.GetNodeSetNames();
  std::map<int, std::vector<int> > const & node_sets = macroscale_mesh.GetNodeSets();
  std::vector<std::string> const & bc_strings = parser_g.GetBoundaryConditionStrings();
  int dim = macroscale_mesh.GetDim();
  std::string const & time_integration_scheme = parser_g.TimeIntegrationScheme();
  boundary_condition_manager_array[qthread_job_id].Initialize(node_set_names, node_sets, bc_strings, dim, time_integration_scheme);
}

static aligned_t QthreadsApplyInitialConditions(void *arg)
{
  using nimble::ModelData;

  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int qthread_job_id = arg_struct->job_id;

  ModelData& model_data = data_manager_array[qthread_job_id].GetMacroScaleData();

  int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
  int velocity_field_id = model_data.GetFieldId("velocity");
  /*const*/ double * const reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
  double* velocity = model_data.GetNodeData(velocity_field_id);
  boundary_condition_manager_array[qthread_job_id].ApplyInitialConditions(Viewify(reference_coordinate,3), Viewify(velocity,3));
}

static aligned_t QthreadsApplyKinematicBC(void *arg)
{
  using nimble::ModelData;

  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int qthread_job_id = arg_struct->job_id;

  ModelData& model_data = data_manager_array[qthread_job_id].GetMacroScaleData();

  int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
  int displacement_field_id = model_data.GetFieldId("displacement");
  int velocity_field_id = model_data.GetFieldId("velocity");
  /*const*/ double * const reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
  double* displacement = model_data.GetNodeData(displacement_field_id);
  double* velocity = model_data.GetNodeData(velocity_field_id);
  boundary_condition_manager_array[qthread_job_id].ApplyKinematicBC(time_current_g, time_previous_g, Viewify(reference_coordinate,3), Viewify(displacement,3), Viewify(velocity,3));
}

static aligned_t QthreadsInitializeExodusOutput(void *arg)
{
  using nimble::ModelData;

  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int mpi_rank = arg_struct->mpi_rank;
  int mpi_num_ranks = arg_struct->mpi_num_ranks;
  int qthread_job_id = arg_struct->job_id;
  int qthread_number_of_workers = qthread_num_workers();

  int name_mpi_rank = mpi_rank;
  int name_mpi_num_ranks = mpi_num_ranks;
  int name_qthread_job_id = qthread_job_id;
  int name_qthread_number_of_workers = qthread_number_of_workers;

  if (parser_g.UseTwoLevelMeshDecomposition() == false) {
    name_mpi_rank = mpi_rank * qthread_number_of_workers + qthread_job_id;
    name_mpi_num_ranks = mpi_num_ranks * qthread_number_of_workers;
    name_qthread_job_id = -1;
    name_qthread_number_of_workers = 0;
  }

  char serial_file_name[nimble::MAX_C_STR_SIZE];
  strcpy(serial_file_name, parser_g.ExodusFileName().c_str());

  char file_name_extension[nimble::MAX_C_STR_SIZE];
  file_name_extension[0] = 'e';
  file_name_extension[1] = '\0';

  char file_label[nimble::MAX_C_STR_SIZE];
  std::string file_label_str("qthreads");
  strcpy(file_label, file_label_str.c_str());

  char exodus_file_name[nimble::MAX_C_STR_SIZE];

  nimble::IOFileNameThreadSafe(serial_file_name,
                               file_name_extension,
                               file_label,
                               name_mpi_rank,
                               name_mpi_num_ranks,
                               name_qthread_job_id,
                               name_qthread_number_of_workers,
                               exodus_file_name);

  ModelData& model_data = data_manager_array[qthread_job_id].GetMacroScaleData();

  int num_global_data = 0;
  std::vector<std::string> global_data_labels(num_global_data);
  std::vector<std::string> const & node_data_labels_for_output = model_data.GetNodeDataLabelsForOutput();
  std::map<int, std::vector<std::string> > const & elem_data_labels_for_output = model_data.GetElementDataLabelsForOutput();
  std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
  exodus_output_array[qthread_job_id].Initialize(exodus_file_name, mesh_array[qthread_job_id]);
  exodus_output_array[qthread_job_id].InitializeDatabase(mesh_array[qthread_job_id],
                                                         global_data_labels,
                                                         node_data_labels_for_output,
                                                         elem_data_labels_for_output,
                                                         derived_elem_data_labels);
}

static aligned_t QthreadsComputeLumpedMass(void *arg)
{
  using nimble::GenesisMesh;
  using nimble::ModelData;
  using nimble::Block;

  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int qthread_job_id = arg_struct->job_id;

  GenesisMesh const & macroscale_mesh = mesh_array[qthread_job_id];
  ModelData& model_data = data_manager_array[qthread_job_id].GetMacroScaleData();

  int num_nodes = macroscale_mesh.GetNumNodes();
  int lumped_mass_field_id = model_data.GetFieldId("lumped_mass");
  int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
  double* lumped_mass = model_data.GetNodeData(lumped_mass_field_id);
  double* reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
  std::map<int, Block> & blocks = model_data.GetBlocks();

  for (auto const & entry : blocks) {
    int block_id = entry.first;
    Block const & block = entry.second;
    int num_elem_in_block = macroscale_mesh.GetNumElementsInBlock(block_id);
    int const * elem_conn = macroscale_mesh.GetConnectivity(block_id);
    block.ComputeLumpedMassMatrix(reference_coordinate,
                                  num_elem_in_block,
                                  elem_conn,
                                  lumped_mass);
  }
}

static aligned_t QthreadsComputeDerivedElementData(void *arg)
{
  using nimble::ModelData;
  using nimble::Block;

  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int qthread_job_id = arg_struct->job_id;

  ModelData& model_data = data_manager_array[qthread_job_id].GetMacroScaleData();

  std::map<int, Block> & blocks = model_data.GetBlocks();
  std::map<int, std::vector<std::string> > const & elem_data_labels = model_data.GetElementDataLabels();
  std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
  int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
  int displacement_field_id = model_data.GetFieldId("displacement");
  const double * const reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
  const double * const displacement = model_data.GetNodeData(displacement_field_id);

  for (auto & entry : blocks) {
    int block_id = entry.first;
    Block & block = entry.second;
    if (derived_element_data_array[qthread_job_id].find(block_id) == derived_element_data_array[qthread_job_id].end()) {
      derived_element_data_array[qthread_job_id][block_id] = std::vector< std::vector<double> >();
    }
    int num_elem_in_block = mesh_array[qthread_job_id].GetNumElementsInBlock(block_id);
    const int * const elem_conn = mesh_array[qthread_job_id].GetConnectivity(block_id);
    std::vector<double> const & elem_data_np1 = model_data.GetElementDataNew(block_id);
    block.ComputeDerivedElementData(reference_coordinate,
                                    displacement,
                                    num_elem_in_block,
                                    elem_conn,
                                    elem_data_labels.at(block_id).size(),
                                    elem_data_np1,
                                    derived_elem_data_labels.at(block_id).size(),
                                    derived_element_data_array[qthread_job_id].at(block_id));
  }
}

static aligned_t QthreadsExodusWriteStep(void *arg)
{
  using nimble::ModelData;

  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int qthread_job_id = arg_struct->job_id;

  ModelData& model_data = data_manager_array[qthread_job_id].GetMacroScaleData();

  int num_global_data = 0;
  std::vector<double> global_data(num_global_data);
  std::vector< std::vector<double> > node_data_for_output;
  model_data.GetNodeDataForOutput(node_data_for_output);
  std::map<int, std::vector<std::string> > const & elem_data_labels_for_output = model_data.GetElementDataLabelsForOutput();
  std::map<int, std::vector< std::vector<double> > > elem_data_for_output;
  model_data.GetElementDataForOutput(elem_data_for_output);
  std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
  exodus_output_array[qthread_job_id].WriteStep(time_current_g,
						    global_data,
						    node_data_for_output,
						    elem_data_labels_for_output,
						    elem_data_for_output,
						    derived_elem_data_labels,
                derived_element_data_array[qthread_job_id]);
}


static aligned_t QthreadsComputeInternalForce(void *arg)
{
  using nimble::GenesisMesh;
  using nimble::ModelData;
  using nimble::Block;

  qthread_arg_t* arg_struct = (qthread_arg_t*)arg;
  int qthread_job_id = arg_struct->job_id;

  GenesisMesh const & macroscale_mesh = mesh_array[qthread_job_id];
  ModelData& model_data = data_manager_array[qthread_job_id].GetMacroScaleData();

  int dim = macroscale_mesh.GetDim();
  int num_nodes = macroscale_mesh.GetNumNodes();
  int num_unknowns = num_nodes * dim;

  std::map<int, Block> & blocks = model_data.GetBlocks();
  std::map<int, std::vector<std::string> > const & elem_data_labels = model_data.GetElementDataLabels();
  std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
  int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
  int displacement_field_id = model_data.GetFieldId("displacement");
  int velocity_field_id = model_data.GetFieldId("velocity");
  int internal_force_field_id = model_data.GetFieldId("internal_force");
  const double * const reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
  const double * const displacement = model_data.GetNodeData(displacement_field_id);
  const double * const velocity = model_data.GetNodeData(velocity_field_id);
  double * const internal_force = model_data.GetNodeData(internal_force_field_id);

  // For explicit dynamics, the macroscale model is never treated as an RVE
  // so rve_macroscale_deformation_gradient will always be the identity matrix
  std::vector<double> rve_macroscale_deformation_gradient(dim*dim, 0.0);
  for (int i=0 ; i<dim ; i++) {
    rve_macroscale_deformation_gradient[i] = 1.0;
  }

  for (int i=0 ; i<num_unknowns ; ++i) {
    internal_force[i] = 0.0;
  }

  for (auto & entry : blocks) {
    int block_id = entry.first;
    int num_elem_in_block = macroscale_mesh.GetNumElementsInBlock(block_id);
    int const * elem_conn = macroscale_mesh.GetConnectivity(block_id);
    std::vector<int> const & elem_global_ids = macroscale_mesh.GetElementGlobalIdsInBlock(block_id);
    Block & block = entry.second;
    std::vector<double> const & elem_data_n = model_data.GetElementDataOld(block_id);
    std::vector<double> & elem_data_np1 = model_data.GetElementDataNew(block_id);
    block.ComputeInternalForce(reference_coordinate,
                               displacement,
                               velocity,
                               rve_macroscale_deformation_gradient.data(),
                               internal_force,
                               time_previous_g,
                               time_current_g,
                               num_elem_in_block,
                               elem_conn,
                               elem_global_ids.data(),
                               elem_data_labels.at(block_id),
                               elem_data_n,
                               elem_data_np1,
                               data_manager_array[qthread_job_id],
                               is_output_step_g);
  }
}

void SharedMemoryScatter(std::string const & variable_name,
                         std::map<int,int> const & global_node_id_to_shared_mem_node_id,
                         int num_threads,
                         std::vector<double>& shared_mem_vector) {

  using nimble::GenesisMesh;
  using nimble::DataManager;

  int dim = 3;

  for (unsigned int i=0 ; i<shared_mem_vector.size() ; ++i) {
    shared_mem_vector[i] = 0.0;
  }

  for(unsigned long i_thread = 0; i_thread < num_threads; i_thread++) {

    GenesisMesh const & macroscale_mesh = mesh_array[i_thread];
    int num_nodes = macroscale_mesh.GetNumNodes();
    int const * const global_node_ids_ptr = macroscale_mesh.GetNodeGlobalIds();
    nimble::ModelData& model_data = data_manager_array[i_thread].GetMacroScaleData();
    int field_id = model_data.GetFieldId(variable_name);
    int field_dimension = LengthToInt(model_data.GetField(field_id).length_, dim);

    double* data = model_data.GetNodeData(field_id);
    for (int i_node=0 ; i_node<num_nodes ; i_node++) {
      int global_node_id = global_node_ids_ptr[i_node];
      int shared_mem_node_id = global_node_id_to_shared_mem_node_id.at(global_node_id);
      for (int i=0 ; i<field_dimension ; i++) {
        shared_mem_vector[field_dimension*shared_mem_node_id+i] += data[field_dimension*i_node+i];
      }
    }
  }
}

void SharedMemoryGather(std::string const & variable_name,
                        std::map<int,int> const & global_node_id_to_shared_mem_node_id,
                        int num_threads,
                        std::vector<double>& shared_mem_vector) {

  using nimble::GenesisMesh;
  using nimble::DataManager;

  int dim = 3;

  for(unsigned long i_thread = 0; i_thread < num_threads; i_thread++) {

    GenesisMesh const & macroscale_mesh = mesh_array[i_thread];
    int num_nodes = macroscale_mesh.GetNumNodes();
    int const * const global_node_ids_ptr = macroscale_mesh.GetNodeGlobalIds();
    nimble::ModelData& model_data = data_manager_array[i_thread].GetMacroScaleData();
    int field_id = model_data.GetFieldId(variable_name);
    int field_dimension = LengthToInt(model_data.GetField(field_id).length_, dim);

    double* data = model_data.GetNodeData(field_id);
    for (int i_node=0 ; i_node<num_nodes ; i_node++) {
      int global_node_id = global_node_ids_ptr[i_node];
      int shared_mem_node_id = global_node_id_to_shared_mem_node_id.at(global_node_id);
      for (int i=0 ; i<field_dimension ; i++) {
        data[field_dimension*i_node+i] = shared_mem_vector[field_dimension*shared_mem_node_id+i];
      }
    }
  }
}

int ExplicitTimeIntegrator(nimble::Parser & parser,
			   nimble::GenesisMesh & mesh,
			   nimble::DataManager & data_manager,
			   nimble::BoundaryConditionManager & boundary_condition_manager,
         nimble::ExodusOutput & exodus_output,
         int num_mpi_ranks,
         int my_mpi_rank);

int QuasistaticTimeIntegrator(nimble::Parser & parser,
			      nimble::GenesisMesh & mesh,
			      nimble::DataManager & data_manager,
			      nimble::BoundaryConditionManager & boundary_condition_manager,
			      nimble::ExodusOutput & exodus_output,
            int num_mpi_ranks,
            int my_mpi_rank);

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  int ierr;
  int num_mpi_ranks;
  int my_mpi_rank;

  // todo: check error codes
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_ranks);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);

  if (argc != 6 && my_mpi_rank == 0) {
    std::cout << "\nUsage:  mpirun -np NP NimbleSM_Qthreads -num_shepherds <num_shepherds> -num_workers_per_shepherd <num_workers_per_shepherd> <input_deck.in>\n" << std::endl;
    exit(1);
  }

  // Parse command line for qthreads settings
  int qt_num_shepherds = 1;
  int qt_num_workers_per_shepherd = 1;
  for (int i_arg=0 ; i_arg<argc-1 ; i_arg++) {
    if (std::string(argv[i_arg]) == "-num_shepherds") {
      qt_num_shepherds = std::atoi(argv[i_arg+1]);
    }
    if(std::string(argv[i_arg]) == "-num_workers_per_shepherd") {
      qt_num_workers_per_shepherd = std::atoi(argv[i_arg+1]);
    }
  }

  int setenv_overwrite = 1;

  std::stringstream qt_num_shepherds_ss;
  qt_num_shepherds_ss << qt_num_shepherds;
  setenv("QT_NUM_SHEPHERDS", qt_num_shepherds_ss.str().c_str(), setenv_overwrite);

  std::stringstream qt_num_workers_per_shepherd_ss;
  qt_num_workers_per_shepherd_ss << qt_num_workers_per_shepherd;
  setenv("QT_NUM_WORKERS_PER_SHEPHERD", qt_num_workers_per_shepherd_ss.str().c_str(), setenv_overwrite);

  int qt_stack_size_power_of_two = 18;
  int qt_stack_size = 2;
  for (int i=0 ; i<qt_stack_size_power_of_two ; i++) {
    qt_stack_size *= 2;
  }
  std::stringstream qt_stack_size_ss;
  qt_stack_size_ss << qt_stack_size;
  setenv("QT_STACK_SIZE", qt_stack_size_ss.str().c_str(), setenv_overwrite);

  ierr = qthread_initialize();
  assert(ierr == QTHREAD_SUCCESS);

  // Banner
  if (my_mpi_rank == 0) {
    std::cout << "\n-- NimbleSM" << std::endl;
    std::cout << "-- version " << nimble::NimbleVersion() << std::endl;
    std::cout << "     " << num_mpi_ranks << " mpi ranks" << std::endl;
    std::cout << "     " << qthread_num_shepherds() << " qthread shepherds per mpi rank" << std::endl;
    std::cout << "     " << qthread_num_workers() << " qthread workers per mpi rank" << std::endl;
    std::cout << "     " << num_mpi_ranks*qthread_num_workers() << " total qthread workers over all mpi ranks" << std::endl;
  }

  int status = 0;

  std::string input_deck_name = argv[argc-1];

  int qt_num_workers = qthread_num_workers();
  aligned_t qt_return_value[qt_num_workers];
  qthread_arg_t qthread_arg_struct[qt_num_workers];

  time_previous_g = time_current_g = 0.0;
  mesh_array = new nimble::GenesisMesh[qt_num_workers];
  rve_mesh_array = new nimble::GenesisMesh[qt_num_workers];
  data_manager_array = new nimble::DataManager[qt_num_workers];
  boundary_condition_manager_array = new nimble::BoundaryConditionManager[qt_num_workers];
  exodus_output_array = new nimble::ExodusOutput[qt_num_workers];
  derived_element_data_array = new std::map<int, std::vector<std::vector<double> > >[qt_num_workers];

  parser_g.Initialize(input_deck_name);
  double final_time = parser_g.FinalTime();
  double delta_time, half_delta_time;
  int num_load_steps = parser_g.NumLoadSteps();
  int output_frequency = parser_g.OutputFrequency();

  for(unsigned long i = 0; i < qt_num_workers; i++) {
    qthread_arg_struct[i].mpi_rank = my_mpi_rank;
    qthread_arg_struct[i].mpi_num_ranks = num_mpi_ranks;
    qthread_arg_struct[i].job_id = i;
  }

  // Read the genesis mesh files
  for(unsigned long i = 0; i < qt_num_workers; i++) {
  	ierr = qthread_fork(QthreadsReadGenesisFiles, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }
  for (int i = 0; i < qt_num_workers; i++) {
    ierr = qthread_readFF(NULL, &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }

  // TODO read rve mesh files

  // Initialize the data manager
  for(unsigned long i = 0; i < qt_num_workers; i++) {
  	ierr = qthread_fork(QthreadsInitializeDataManager, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }
  for (int i = 0; i < qt_num_workers; i++) {
    ierr = qthread_readFF(NULL, &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }

  // Create shared memory data containers
  std::set<int> global_node_ids_shared_mem_set;
  for(unsigned long i = 0; i < qt_num_workers; i++) {
    nimble::GenesisMesh const & macroscale_mesh = mesh_array[i];
    int num_nodes = macroscale_mesh.GetNumNodes();
    int const * const global_node_ids_ptr = macroscale_mesh.GetNodeGlobalIds();
    for (int i_node=0 ; i_node<num_nodes ; i_node++) {
      global_node_ids_shared_mem_set.insert(global_node_ids_ptr[i_node]);
    }
  }
  std::vector<int> global_node_ids_shared_mem(global_node_ids_shared_mem_set.begin(), global_node_ids_shared_mem_set.end());
  std::map<int, int> global_node_id_to_shared_mem_node_id;
  for (unsigned int i=0 ; i<global_node_ids_shared_mem.size() ; i++) {
    global_node_id_to_shared_mem_node_id[global_node_ids_shared_mem[i]] = i;
  }
  std::vector<double> lumped_mass_shared_mem(global_node_ids_shared_mem.size(), 0.0);
  std::vector<double> internal_force_shared_mem(3*global_node_ids_shared_mem.size(), 0.0);

  // Initialize the MPI container
  nimble::MPIContainer mpi_container;
  int mpi_reduction_version = parser_g.ReductionVersion();
  mpi_container.Initialize(global_node_ids_shared_mem, mpi_reduction_version);

  // Initialize the boundary condition manager
  for(unsigned long i = 0; i < qt_num_workers; i++) {
    ierr = qthread_fork(QthreadsInitializeBoundaryConditionManager, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }
  for (int i = 0; i < qt_num_workers; i++) {
    ierr = qthread_readFF(NULL, &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }

  // Initialize the exodus output files
  for(unsigned long i = 0; i < qt_num_workers; i++) {
  	ierr = qthread_fork(QthreadsInitializeExodusOutput, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }
  for (int i = 0; i < qt_num_workers; i++) {
    ierr = qthread_readFF(NULL, &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }

  // Compute lumped mass
  for(unsigned long i = 0; i < qt_num_workers; i++) {
  	ierr = qthread_fork(QthreadsComputeLumpedMass, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }
  for (int i = 0; i < qt_num_workers; i++) {
    ierr = qthread_readFF(NULL, &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }

  // Move lumped mass from the per-thread data structures to a per-MPI-rank data structure
  SharedMemoryScatter("lumped_mass",
                      global_node_id_to_shared_mem_node_id,
                      qt_num_workers,
                      lumped_mass_shared_mem);

  // Perform a reduction to obtain correct values on MPI boundaries
  int scalar_dimension = 1;
  mpi_container.VectorReduction(scalar_dimension, lumped_mass_shared_mem.data());

  // Move the properly-summed data from the per-MPI-rank data structure to the per-thread data structures
  SharedMemoryGather("lumped_mass",
                     global_node_id_to_shared_mem_node_id,
                     qt_num_workers,
                     lumped_mass_shared_mem);

  // TODO Compute critical time step

  // Apply initial conditions
  for(unsigned long i = 0; i < qt_num_workers; i++) {
  	ierr = qthread_fork(QthreadsApplyInitialConditions, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }
  for (int i = 0; i < qt_num_workers; i++) {
    ierr = qthread_readFF(NULL, &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }

  // Apply kinematic boundary conditions
  for(unsigned long i = 0; i < qt_num_workers; i++) {
  	ierr = qthread_fork(QthreadsApplyKinematicBC, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }
  for (int i = 0; i < qt_num_workers; i++) {
    ierr = qthread_readFF(NULL, &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }

  // Compute derived element data
  for(unsigned long i = 0; i < qt_num_workers; i++) {
  	ierr = qthread_fork(QthreadsComputeDerivedElementData, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }
  for (int i = 0; i < qt_num_workers; i++) {
    ierr = qthread_readFF(NULL, &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }

  // Write initial state to the exodus files
  for(unsigned long i = 0; i < qt_num_workers; i++) {
  	ierr = qthread_fork(QthreadsExodusWriteStep, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }
  for (int i = 0; i < qt_num_workers; i++) {
    ierr = qthread_readFF(NULL, &qt_return_value[i]);
    assert(ierr == QTHREAD_SUCCESS);
  }

  double user_specified_time_step = final_time/num_load_steps;
  if (my_mpi_rank == 0) {
    std::cout << "\nUser specified time step:              " << user_specified_time_step << std::endl;
  //std::cout << "Approximate maximum stable time step:  " << critical_time_step << "\n" << std::endl;
  // if (user_specified_time_step > critical_time_step) {
  //   std::cout << "**** WARNING:  The user specified time step exceeds the computed maximum stable time step.\n" << std::endl;
  // }
  std::cout << "\nExplicit time integration:\n    0% complete" << std::endl;
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
    is_output_step_g = false;
    if (step%output_frequency == 0 || step == num_load_steps - 1) {
      is_output_step_g = true;
    }

    time_previous_g = time_current_g;
    time_current_g += final_time/num_load_steps;
    double delta_time = time_current_g - time_previous_g;
    double half_delta_time = 0.5*delta_time;

    // V^{n+1/2} = V^{n} + (dt/2) * A^{n}
    for(unsigned long i_thread = 0; i_thread < qt_num_workers; i_thread++) {
      int num_unknowns = mesh_array[i_thread].GetNumNodes() * mesh_array[i_thread].GetDim();
      nimble::ModelData& model_data = data_manager_array[i_thread].GetMacroScaleData();
      int velocity_field_id = model_data.GetFieldId("velocity");
      int acceleration_field_id = model_data.GetFieldId("acceleration");
      double * const velocity = model_data.GetNodeData(velocity_field_id);
      const double * const acceleration = model_data.GetNodeData(acceleration_field_id);
      for (int i=0 ; i<num_unknowns ; ++i) {
        velocity[i] += half_delta_time * acceleration[i];
      }
    }

    // Apply kinematic boundary conditions
    for(unsigned long i = 0; i < qt_num_workers; i++) {
      ierr = qthread_fork(QthreadsApplyKinematicBC, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
      assert(ierr == QTHREAD_SUCCESS);
    }
    for (int i = 0; i < qt_num_workers; i++) {
      ierr = qthread_readFF(NULL, &qt_return_value[i]);
      assert(ierr == QTHREAD_SUCCESS);
    }

    // // Evaluate external body forces
    // for (int i=0 ; i<num_unknowns ; ++i) {
    //   external_force[i] = 0.0;
    // }

    // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
    for(unsigned long i_thread = 0; i_thread < qt_num_workers; i_thread++) {
      int num_unknowns = mesh_array[i_thread].GetNumNodes() * mesh_array[i_thread].GetDim();
      nimble::ModelData& model_data = data_manager_array[i_thread].GetMacroScaleData();
      int displacement_field_id = model_data.GetFieldId("displacement");
      int velocity_field_id = model_data.GetFieldId("velocity");
      double * const displacement = model_data.GetNodeData(displacement_field_id);
      const double * const velocity = model_data.GetNodeData(velocity_field_id);
      for (int i=0 ; i<num_unknowns ; ++i) {
        displacement[i] += delta_time * velocity[i];
      }
    }


    // // Evaluate the internal force
    for(unsigned long i = 0; i < qt_num_workers; i++) {
      ierr = qthread_fork(QthreadsComputeInternalForce, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
      assert(ierr == QTHREAD_SUCCESS);
    }
    for (int i = 0; i < qt_num_workers; i++) {
      ierr = qthread_readFF(NULL, &qt_return_value[i]);
      assert(ierr == QTHREAD_SUCCESS);
    }


    // Move internal force from the per-thread data structures to a per-MPI-rank data structure
    SharedMemoryScatter("internal_force",
                        global_node_id_to_shared_mem_node_id,
                        qt_num_workers,
                        internal_force_shared_mem);

    // Perform a reduction to obtain correct values on MPI boundaries
    int vector_dimension = 3;
    mpi_container.VectorReduction(vector_dimension, internal_force_shared_mem.data());

    // Move the properly-summed data from the per-MPI-rank data structure to the per-thread data structures
    SharedMemoryGather("internal_force",
                       global_node_id_to_shared_mem_node_id,
                       qt_num_workers,
                       internal_force_shared_mem);

    // fill acceleration vector A^{n+1} = M^{-1} ( F^{n} + b^{n} )
    for(unsigned long i_thread = 0; i_thread < qt_num_workers; i_thread++) {
      int num_unknowns = mesh_array[i_thread].GetNumNodes() * mesh_array[i_thread].GetDim();
      nimble::ModelData& model_data = data_manager_array[i_thread].GetMacroScaleData();
      int lumped_mass_field_id = model_data.GetFieldId("lumped_mass");
      int internal_force_field_id = model_data.GetFieldId("internal_force");
      int external_force_field_id = model_data.GetFieldId("external_force");
      int acceleration_field_id = model_data.GetFieldId("acceleration");
      const double * const lumped_mass = model_data.GetNodeData(lumped_mass_field_id);
      const double * const internal_force = model_data.GetNodeData(internal_force_field_id);
      const double * const external_force = model_data.GetNodeData(external_force_field_id);
      double * const acceleration = model_data.GetNodeData(acceleration_field_id);
      for (int i=0 ; i<num_unknowns ; ++i) {
        acceleration[i] = (1.0/lumped_mass[i/3]) * (internal_force[i] + external_force[i]);
      }
    }

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    for(unsigned long i_thread = 0; i_thread < qt_num_workers; i_thread++) {
      int num_unknowns = mesh_array[i_thread].GetNumNodes() * mesh_array[i_thread].GetDim();
      nimble::ModelData& model_data = data_manager_array[i_thread].GetMacroScaleData();
      int acceleration_field_id = model_data.GetFieldId("acceleration");
      int velocity_field_id = model_data.GetFieldId("velocity");
      const double * const acceleration = model_data.GetNodeData(acceleration_field_id);
      double * const velocity = model_data.GetNodeData(velocity_field_id);
      for (int i=0 ; i<num_unknowns ; ++i) {
        velocity[i] += half_delta_time * acceleration[i];
      }
    }

    if (is_output_step_g) {

      // Compute derived element data
      for(unsigned long i = 0; i < qt_num_workers; i++) {
        ierr = qthread_fork(QthreadsComputeDerivedElementData, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
        assert(ierr == QTHREAD_SUCCESS);
      }
      for (int i = 0; i < qt_num_workers; i++) {
        ierr = qthread_readFF(NULL, &qt_return_value[i]);
        assert(ierr == QTHREAD_SUCCESS);
      }

      // Apply kinematic boundary conditions
      for(unsigned long i = 0; i < qt_num_workers; i++) {
        ierr = qthread_fork(QthreadsApplyKinematicBC, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
        assert(ierr == QTHREAD_SUCCESS);
      }
      for (int i = 0; i < qt_num_workers; i++) {
        ierr = qthread_readFF(NULL, &qt_return_value[i]);
        assert(ierr == QTHREAD_SUCCESS);
      }

      // Write initial state to the exodus files
      for(unsigned long i = 0; i < qt_num_workers; i++) {
        ierr = qthread_fork(QthreadsExodusWriteStep, (void*)&(qthread_arg_struct[i]), &qt_return_value[i]);
        assert(ierr == QTHREAD_SUCCESS);
      }
      for (int i = 0; i < qt_num_workers; i++) {
        ierr = qthread_readFF(NULL, &qt_return_value[i]);
        assert(ierr == QTHREAD_SUCCESS);
      }
    }
  }

  qthread_finalize();
  MPI_Finalize();

  delete[] mesh_array;
  delete[] rve_mesh_array;
  delete[] data_manager_array;
  delete[] exodus_output_array;
  delete[] derived_element_data_array;

  exit(0);
}
