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

#include "nimble_darma_utils.h"
#include "nimble_utils.h"
#include "nimble_view.h"
#include "nimble_material_factory.h"
#include <iomanip>

#ifdef NIMBLE_HAVE_EXTRAS
#include "nimble_extras_material_factory.h"
#endif

using namespace darma_runtime;

namespace nimble {

  void PrintProgress(ProgressBarFlag progress_bar_flag,
                     int step,
                     int num_steps) {
    if (progress_bar_flag == FIRST_STEP) {
      std::cout << "Explicit time integration:\n    0% complete" << std::endl << std::flush;
    }
    else if(progress_bar_flag == PRINT_PROGRESS) {
      std::cout << "   " << static_cast<int>( 100.0 * static_cast<double>(step+1)/num_steps ) << "% complete" << std::endl << std::flush;
    }
    else if(progress_bar_flag == LAST_STEP) {
      std::cout << "  100% complete\n" << std::endl << std::flush;
    }
  }

  void ReadGenesisFiles::operator()(
    Index1D<int> index,
    AccessHandleCollection<GenesisMesh, Range1D<int>> macroscale_mesh_collection,
    AccessHandleCollection<GenesisMesh, Range1D<int>> rve_mesh_collection,
    Parser parser
  ) const {

    int num_ranks = index.max_value + 1;
    int my_rank = index.value;
    int min_rank = index.min_value;

    if (my_rank == min_rank) {
      std::cout << "\nReading genesis file" << std::endl << std::flush;
    }

    auto macroscale_mesh_handle = macroscale_mesh_collection[index].local_access();
    std::string genesis_file_name = IOFileName(parser.GenesisFileName(), "g", "", my_rank, num_ranks);
    macroscale_mesh_handle.get_reference().ReadFile(genesis_file_name);

    auto rve_mesh_handle = rve_mesh_collection[index].local_access();
    std::string rve_genesis_file_name = IOFileName(parser.RVEGenesisFileName(), "g");
    rve_mesh_handle.get_reference().ReadFile(rve_genesis_file_name);

    if (my_rank == min_rank) {
      std::cout << "  complete" << std::endl << std::flush;
    }
  }

  void InitializeDataManager::operator()(
    Index1D<int> index,
    AccessHandleCollection<GenesisMesh, Range1D<int>> macroscale_mesh_collection,
    AccessHandleCollection<GenesisMesh, Range1D<int>> rve_mesh_collection,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
    Parser parser
  ) const {

#ifdef NIMBLE_HAVE_EXTRAS
  using MaterialFactoryType = nimble::ExtrasMaterialFactory;
#else
  using MaterialFactoryType = nimble::MaterialFactory;
#endif

    int num_ranks = index.max_value + 1;
    int my_rank = index.value;
    int min_rank = index.min_value;

    if (my_rank == min_rank) {
      std::cout << "\nInitializing DataManager (num nodes = " << macroscale_mesh_collection[index].local_access().get_value().GetNumNodes() << ")" << std::endl << std::flush;
    }

    auto macroscale_mesh_handle = macroscale_mesh_collection[index].local_access();
    auto rve_mesh_handle = rve_mesh_collection[index].local_access();
    auto data_manager_handle = data_manager_collection[index].local_access();

    int dim = macroscale_mesh_handle.get_value().GetDim();
    int num_nodes = macroscale_mesh_handle.get_value().GetNumNodes();
    int num_blocks = macroscale_mesh_handle.get_value().GetNumBlocks();
    std::vector<int> block_ids = macroscale_mesh_handle.get_value().GetBlockIds();

    bool rve_is_valid = rve_mesh_handle.get_value().IsValid();

    // Global data
    int num_global_data = 0;
    std::vector<std::string> global_data_labels(num_global_data);
    std::vector<double> global_data(num_global_data);

    ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();
    model_data.SetDimension(dim);

    int lumped_mass_field_id = model_data.AllocateNodeData(SCALAR, "lumped_mass", num_nodes);
    int reference_coordinate_field_id = model_data.AllocateNodeData(VECTOR, "reference_coordinate", num_nodes);
    int displacement_field_id = model_data.AllocateNodeData(VECTOR, "displacement", num_nodes);
    int velocity_field_id =  model_data.AllocateNodeData(VECTOR, "velocity", num_nodes);
    int acceleration_field_id =  model_data.AllocateNodeData(VECTOR, "acceleration", num_nodes);
    int internal_force_field_id =  model_data.AllocateNodeData(VECTOR, "internal_force", num_nodes);
    int external_force_field_id =  model_data.AllocateNodeData(VECTOR, "external_force", num_nodes);

    const double * const ref_coord_x = macroscale_mesh_handle.get_value().GetCoordinatesX();
    const double * const ref_coord_y = macroscale_mesh_handle.get_value().GetCoordinatesY();
    const double * const ref_coord_z = macroscale_mesh_handle.get_value().GetCoordinatesZ();
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
      int num_elem_in_block = macroscale_mesh_handle.get_value().GetNumElementsInBlock(block_id);
      std::string const & macro_material_parameters = parser.GetMacroscaleMaterialParameters(block_id);
      std::map<int, std::string> const & rve_material_parameters = parser.GetMicroscaleMaterialParameters();
      std::string rve_bc_strategy = parser.GetMicroscaleBoundaryConditionStrategy();
      blocks[block_id] = Block();
      MaterialFactoryType factory;
      blocks[block_id].Initialize(macro_material_parameters, rve_material_parameters, rve_mesh_handle.get_value(), rve_bc_strategy, factory);
      std::vector< std::pair<std::string, Length> > data_labels_and_lengths;
      blocks[block_id].GetDataLabelsAndLengths(data_labels_and_lengths);
      model_data.DeclareElementData(block_id, data_labels_and_lengths);
    }
    std::map<int, int> num_elem_in_each_block = macroscale_mesh_handle.get_value().GetNumElementsInBlock();
    model_data.AllocateElementData(num_elem_in_each_block);
    model_data.SpecifyOutputFields(parser.GetOutputFieldString());
    std::map<int, std::vector<std::string> > const & elem_data_labels = model_data.GetElementDataLabels();
    std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
    std::vector<int> rve_output_elem_ids = parser.MicroscaleOutputElementIds();

    for (auto & entry : blocks) {
      int block_id = entry.first;
      Block& block = entry.second;
      int num_elem_in_block = macroscale_mesh_handle.get_value().GetNumElementsInBlock(block_id);
      std::vector<int> const & elem_global_ids = macroscale_mesh_handle.get_value().GetElementGlobalIdsInBlock(block_id);
      std::vector<double> & elem_data_n = model_data.GetElementDataOld(block_id);
      std::vector<double> & elem_data_np1 = model_data.GetElementDataNew(block_id);
      MaterialFactoryType factory;
      block.InitializeElementData(num_elem_in_block,
                                  elem_global_ids,
                                  rve_output_elem_ids,
                                  elem_data_labels.at(block_id),
                                  derived_elem_data_labels.at(block_id),
                                  elem_data_n,
                                  elem_data_np1,
                                  factory,
                                  data_manager_handle.get_reference());
    }

    if (my_rank == min_rank) {
      std::cout << "  complete" << std::endl << std::flush;
    }
  }

  void IdentifyGloballySharedNodes::operator()(
    Index1D<int> index,
    AccessHandleCollection<GenesisMesh, Range1D<int>> mesh_collection,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> boundary_node_global_ids_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> boundary_ranks_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> my_global_node_ids_collection
  ) const {

    using keyword_arguments_for_publication::version;
    using keyword_arguments_for_publication::n_readers;

    int num_ranks = index.max_value + 1;
    int my_rank = index.value;
    int min_rank = index.min_value;

    if (num_ranks == 1) {
      return;
    }

     if (my_rank == min_rank) {
       std::cout << "\nIdentifying globally-shared nodes (num nodes = " << mesh_collection[index].local_access().get_value().GetNumNodes() << ")" << std::endl << std::flush;
     }

    auto data_manager_handle = data_manager_collection[index].local_access();
    auto boundary_node_global_ids_handle = boundary_node_global_ids_collection[index].local_access();
    auto boundary_ranks_handle = boundary_ranks_collection[index].local_access();
    auto my_global_node_ids_handle = my_global_node_ids_collection[index].local_access();
    auto mesh_handle = mesh_collection[index].local_access();

    // Create a list of all the global node ids on this rank
    int num_nodes = mesh_handle.get_value().GetNumNodes();
    const int * const global_node_ids = mesh_handle.get_value().GetNodeGlobalIds();
    my_global_node_ids_handle.get_reference().resize(num_nodes);
    for (int i=0 ; i<num_nodes ; ++i) {
      my_global_node_ids_handle.get_reference().at(i) = global_node_ids[i];
    }

    my_global_node_ids_handle.publish(version=0, n_readers=num_ranks-1);

    // todo Be smarter about this n^2 search, perhaps some bounding box stuff or skinning

    // Read the lists of global node ids from all other processors
    // and record nodes that are present on both this rank and another rank
    for (int rank=0 ; rank<num_ranks ; rank++) {
      if (rank != my_rank) {
        auto off_processor_global_node_ids_handle = my_global_node_ids_collection[rank].read_access(version=0);
        create_work(reads(my_global_node_ids_handle, off_processor_global_node_ids_handle), [
           my_global_node_ids_handle, off_processor_global_node_ids_handle, boundary_node_global_ids_handle,
           boundary_ranks_handle, rank
        ]{
            std::vector<int> const & my_global_node_ids = my_global_node_ids_handle.get_value();
            std::vector<int> const & off_processor_global_node_ids = off_processor_global_node_ids_handle.get_value();
            std::vector<int> & boundary_node_global_ids = boundary_node_global_ids_handle.get_reference();
            std::vector<int> & boundary_ranks = boundary_ranks_handle.get_reference();
            for (auto const & node_id : my_global_node_ids) {
              if (std::find(off_processor_global_node_ids.begin(), off_processor_global_node_ids.end(), node_id) != off_processor_global_node_ids.end() ) {
                if (std::find(boundary_node_global_ids.begin(), boundary_node_global_ids.end(), node_id) == boundary_node_global_ids.end()) {
                  boundary_node_global_ids.push_back(node_id);
                }
                if (std::find(boundary_ranks.begin(), boundary_ranks.end(), rank) == boundary_ranks.end()) {
                  boundary_ranks.push_back(rank);
                }
              }
            }
          });
      }
    }

    // Store a map from global node id to local node it
    create_work(reads(mesh_handle), [
        mesh_handle, data_manager_handle, num_nodes
    ]{
        const int * const global_node_ids_cw = mesh_handle.get_value().GetNodeGlobalIds();
        std::map<int, int>& global_node_id_to_local_node_id = data_manager_handle.get_reference().GetMacroScaleData().GetGlobalNodeIdToLocalNodeIdMap();
        for (int i=0 ; i<num_nodes ; i++) {
          global_node_id_to_local_node_id[global_node_ids_cw[i]] = i;
        }
      });

    if (my_rank == min_rank) {
      std::cout << "  complete" << std::endl << std::flush;
    }
  }

  void InitializeBoundaryConditionManager::operator()(
     Index1D<int> index,
     AccessHandleCollection<GenesisMesh, Range1D<int>> mesh_collection,
     AccessHandleCollection<BoundaryConditionManager, Range1D<int>> boundary_condition_manager_collection,
     Parser parser
   ) const {

    auto mesh_handle = mesh_collection[index].local_access();
    auto boundary_condition_manager_handle = boundary_condition_manager_collection[index].local_access();

    std::map<int, std::string> const & node_set_names = mesh_handle.get_value().GetNodeSetNames();
    std::map<int, std::vector<int> > const & node_sets = mesh_handle.get_value().GetNodeSets();
    std::vector<std::string> const & bc_strings = parser.GetBoundaryConditionStrings();
    int dim = mesh_handle.get_value().GetDim();
    std::string const & time_integration_scheme = parser.TimeIntegrationScheme();
    boundary_condition_manager_handle.get_reference().Initialize(node_set_names, node_sets, bc_strings, dim, time_integration_scheme);
  }

   void InitializeExodusOutput::operator()(
     Index1D<int> index,
     AccessHandleCollection<GenesisMesh, Range1D<int>> macroscale_mesh_collection,
     AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
     AccessHandleCollection<ExodusOutput, Range1D<int>> exodus_output_collection,
     Parser parser
   ) const {

     int num_ranks = index.max_value + 1;
     int my_rank = index.value;
     int min_rank = index.min_value;

     if (my_rank == min_rank) {
       std::cout << "\nInitializing exodus output (num nodes = " << macroscale_mesh_collection[index].local_access().get_value().GetNumNodes() << ")" << std::endl << std::flush;
     }

     auto macroscale_mesh_handle = macroscale_mesh_collection[index].local_access();
     auto data_manager_handle = data_manager_collection[index].local_access();
     auto exodus_output_handle = exodus_output_collection[index].local_access();

     ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();

     std::string exodus_file_name = IOFileName(parser.GenesisFileName(), "e", "darma", my_rank, num_ranks);
     int num_global_data = 0;
     std::vector<std::string> global_data_labels(num_global_data);
     std::vector<std::string> const & node_data_labels_for_output = model_data.GetNodeDataLabelsForOutput();
     std::map<int, std::vector<std::string> > const & elem_data_labels_for_output = model_data.GetElementDataLabelsForOutput();
     std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
     exodus_output_handle.get_reference().Initialize(exodus_file_name, macroscale_mesh_handle.get_value());
     exodus_output_handle.get_reference().InitializeDatabase(macroscale_mesh_handle.get_value(),
		  				     global_data_labels,
		  				     node_data_labels_for_output,
		  				     elem_data_labels_for_output,
		  				     derived_elem_data_labels);

     if (my_rank == min_rank) {
       std::cout << "  complete" << std::endl << std::flush;
     }
  }

  void ComputeLumpedMass::operator()(
    Index1D<int> index,
    AccessHandleCollection<GenesisMesh, Range1D<int>> macroscale_mesh_collection,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection
  ) const {

    int my_rank = index.value;
    int min_rank = index.min_value;

    if (my_rank == min_rank) {
      std::cout << "\nComputing lumped mass" << std::endl << std::flush;
    }

    auto macroscale_mesh_handle = macroscale_mesh_collection[index].local_access();
    auto data_manager_handle = data_manager_collection[index].local_access();

    ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();

    int num_nodes = macroscale_mesh_handle.get_value().GetNumNodes();
    int lumped_mass_field_id = model_data.GetFieldId("lumped_mass");
    int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
    double* lumped_mass = model_data.GetNodeData(lumped_mass_field_id);
    double* reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
    std::map<int, Block> & blocks = model_data.GetBlocks();

    for (auto const & entry : blocks) {
      int block_id = entry.first;
      Block const & block = entry.second;
      int num_elem_in_block = macroscale_mesh_handle.get_value().GetNumElementsInBlock(block_id);
      int const * elem_conn = macroscale_mesh_handle.get_value().GetConnectivity(block_id);
      block.ComputeLumpedMassMatrix(reference_coordinate,
                                    num_elem_in_block,
                                    elem_conn,
                                    lumped_mass);
    }

    if (my_rank == min_rank) {
      std::cout << "  complete" << std::endl << std::flush;
    }
  }

  void DistributedVectorReduction::operator()(
    Index1D<int> index,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> boundary_node_global_ids_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> boundary_ranks_collection,
    AccessHandleCollection<std::vector<double>, Range1D<int>> global_reduction_buffer_collection,
    DistributedVectorReductionQuantity distributed_vector_reduction_quantity,
    int step
  ) const {

    using keyword_arguments_for_publication::version;
    using keyword_arguments_for_publication::n_readers;

    int num_ranks = index.max_value + 1;
    int my_rank = index.value;
    int min_rank = index.min_value;

    if (num_ranks == 1) {
      return;
    }

    // if (my_rank == min_rank) {
    //   std::cout << "\nVector reduction on min rank" << std::endl << std::flush;
    // }
    //std::cout << "\nVector reduction on rank " << my_rank << std::endl << std::flush;

    auto data_manager_handle = data_manager_collection[index].local_access();
    auto boundary_node_global_ids_handle = boundary_node_global_ids_collection[index].local_access();
    auto boundary_ranks_handle = boundary_ranks_collection[index].local_access();
    auto global_reduction_buffer_handle = global_reduction_buffer_collection[index].local_access();

    unsigned int num_neighboring_ranks = boundary_ranks_handle.get_value().size();

    // todo: using enum here because compiler won't let me pass a std::string as a functor argument
    std::string quantity_label = "none";
    Length quantity_length = UNDEFINED_LENGTH;
    if (distributed_vector_reduction_quantity == LUMPED_MASS) {
      quantity_label = "lumped_mass";
      quantity_length = SCALAR;
    }
    else if(distributed_vector_reduction_quantity == INTERNAL_FORCE) {
      quantity_label = "internal_force";
      quantity_length = VECTOR;
    }

    ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();

    // Load the local values for the shared nodes into a buffer that all ranks can access
    int dim = model_data.GetDimension();
    int len = LengthToInt(quantity_length, dim);
    std::map<int, int> const & global_node_id_to_local_node_id = model_data.GetGlobalNodeIdToLocalNodeIdMap();
    int num_shared_nodes = boundary_node_global_ids_handle.get_value().size();
    global_reduction_buffer_handle.get_reference().resize(num_shared_nodes * len, 0.0);
    int field_id = model_data.GetFieldId(quantity_label);
    const double * const data = model_data.GetNodeData(field_id);
    for (auto i = 0 ; i < num_shared_nodes ; i++) {
      int global_node_id = boundary_node_global_ids_handle.get_value().at(i);
      int local_node_id = global_node_id_to_local_node_id.at(global_node_id);
      for (int j = 0 ; j < len ; j++) {
        global_reduction_buffer_handle.get_reference().at(len*i + j) = data[len*local_node_id + j];
      }
    }

    boundary_node_global_ids_handle.publish(version=step, n_readers=num_neighboring_ranks); // todo: does not need to be communicated every time, should do in initialization phase
    global_reduction_buffer_handle.publish(version=step, n_readers=num_neighboring_ranks);

    // Sum the values from the other ranks into the local data
    for (auto const & boundary_rank : boundary_ranks_handle.get_value()) {
      auto off_processor_global_ids_handle = boundary_node_global_ids_collection[boundary_rank].read_access(version=step);
      auto off_processor_global_reduction_buffer_handle = global_reduction_buffer_collection[boundary_rank].read_access(version=step);
      create_work(reads(off_processor_global_ids_handle, off_processor_global_reduction_buffer_handle), [
          off_processor_global_reduction_buffer_handle, off_processor_global_ids_handle, data_manager_handle,
          quantity_length, field_id
      ]{
          std::vector<double> const & off_processor_data = off_processor_global_reduction_buffer_handle.get_value();
          std::vector<int> const & off_processor_global_ids = off_processor_global_ids_handle.get_value();
          ModelData& on_processor_model_data = data_manager_handle.get_reference().GetMacroScaleData();
          int dim = on_processor_model_data.GetDimension();
          int len = LengthToInt(quantity_length, dim);
          std::map<int, int> const & global_node_id_to_local_node_id = on_processor_model_data.GetGlobalNodeIdToLocalNodeIdMap();
          double* data = on_processor_model_data.GetNodeData(field_id);
          for (auto i = 0 ; i < off_processor_global_ids.size() ; i++) {
            int global_node_id = off_processor_global_ids.at(i);
            if (global_node_id_to_local_node_id.find(global_node_id) != global_node_id_to_local_node_id.end()) {
              int local_node_id = global_node_id_to_local_node_id.at(global_node_id);
              for (int j = 0 ; j < len ; j++) {
                data[len*local_node_id + j] += off_processor_data.at(len*i + j);
              }
            }
          }
        });
    }

    //if (my_rank == min_rank) {
    //  std::cout << "  complete" << std::endl << std::flush;
    // }
    //std::cout << "  complete rank " << my_rank << std::endl << std::flush;
  }

  void ComputeCriticalTimeStep::operator()(
    Index1D<int> index,
    AccessHandleCollection<GenesisMesh, Range1D<int>> macroscale_mesh_collection,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
    Parser parser
  ) const {

    using keyword_arguments_for_collectives::tag;
    using keyword_arguments_for_collectives::in_out;
    using keyword_arguments_for_collectives::piece;
    using keyword_arguments_for_collectives::n_pieces;

    int num_ranks = index.max_value + 1;
    int my_rank = index.value;
    int min_rank = index.min_value;

    if (my_rank == min_rank) {
      std::cout << "\nComputing critical time step" << std::endl << std::flush;
    }

    auto macroscale_mesh_handle = macroscale_mesh_collection[index].local_access();
    auto data_manager_handle = data_manager_collection[index].local_access();

    ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();

    int displacement_field_id = model_data.GetFieldId("displacement");
    const double * const displacement = model_data.GetNodeData(displacement_field_id);
    int coordinate_field_id = model_data.GetFieldId("reference_coordinate");
    const double * const reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
    std::map<int, Block> & blocks = model_data.GetBlocks();

    double user_time_step = parser.FinalTime()/parser.NumLoadSteps();
    double max_stable_time_step = std::numeric_limits<double>::max();

    for (auto const & entry : blocks) {
      int block_id = entry.first;
      Block const & block = entry.second;
      int num_elem_in_block = macroscale_mesh_handle.get_value().GetNumElementsInBlock(block_id);
      const int * const elem_conn = macroscale_mesh_handle.get_value().GetConnectivity(block_id);
      double block_critical_time_step = block.ComputeCriticalTimeStep(reference_coordinate,
                                                                      displacement,
                                                                      num_elem_in_block,
                                                                      elem_conn);
      if (block_critical_time_step < max_stable_time_step) {
        max_stable_time_step = block_critical_time_step;
      }
    }

    auto global_max_stable_time_step = initial_access<double>("global max stable time step", my_rank);
    create_work([
        global_max_stable_time_step, max_stable_time_step
    ]{
      global_max_stable_time_step.set_value(max_stable_time_step);
    });

    if (num_ranks > 1) {
      std::string tag_str("all reduce for critical time step");
      allreduce<ReductionMin>(tag=tag_str, in_out=global_max_stable_time_step, piece=my_rank, n_pieces=num_ranks);
    }

    create_work(reads(global_max_stable_time_step), [
        global_max_stable_time_step, data_manager_handle, min_rank, my_rank, user_time_step
    ]{
        double time_step = global_max_stable_time_step.get_value();
        data_manager_handle.get_reference().GetMacroScaleData().SetCriticalTimeStep(time_step);
        if (my_rank == min_rank) {
          std::cout << "\nUser specified time step:              " << user_time_step << std::endl;
          std::cout << "Approximate maximum stable time step:  " << time_step << "\n" << std::endl;
          if (user_time_step > time_step) {
            std::cout << "**** WARNING:  The user specified time step exceeds the computed maximum stable time step.\n" << std::endl;
          }
        }
      });

    if (my_rank == min_rank) {
      std::cout << "  complete" << std::endl << std::flush;
    }
  }

  void ApplyInitialConditions::operator()(
    Index1D<int> index,
    AccessHandleCollection<BoundaryConditionManager, Range1D<int>> boundary_condition_manager_collection,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection
  ) const {

    auto boundary_condition_manager_handle = boundary_condition_manager_collection[index].local_access();
    auto data_manager_handle = data_manager_collection[index].local_access();

    ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();

    int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
    int velocity_field_id = model_data.GetFieldId("velocity");
    const double * const reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
    double* velocity = model_data.GetNodeData(velocity_field_id);
    boundary_condition_manager_handle.get_reference().ApplyInitialConditions(Viewify(reference_coordinate,3), Viewify(velocity,3));
  }

  void ApplyKinematicBC::operator()(
    Index1D<int> index,
    AccessHandleCollection<BoundaryConditionManager, Range1D<int>> boundary_condition_manager_collection,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
    double time_current,
    double time_previous
  ) const {

    int my_rank = index.value;
    int min_rank = index.min_value;

    //    if (my_rank == min_rank) {
    //      std::cout << "\nApplying kinematic boundary conditions" << std::endl << std::flush;
    //    }

    auto boundary_condition_manager_handle = boundary_condition_manager_collection[index].local_access();
    auto data_manager_handle = data_manager_collection[index].local_access();

    ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();

    int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
    int displacement_field_id = model_data.GetFieldId("displacement");
    int velocity_field_id = model_data.GetFieldId("velocity");
    const double * const reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
    double* displacement = model_data.GetNodeData(displacement_field_id);
    double* velocity = model_data.GetNodeData(velocity_field_id);
    boundary_condition_manager_handle.get_reference().ApplyKinematicBC(time_current, time_previous, Viewify(reference_coordinate,3), Viewify(displacement,3), Viewify(velocity,3));

    //    if (my_rank == min_rank) {
    //      std::cout << "  complete" << std::endl << std::flush;
    //    }
  }

  void ComputeDerivedElementData::operator()(
    Index1D<int> index,
    AccessHandleCollection<GenesisMesh, Range1D<int>> macroscale_mesh_collection,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
    AccessHandleCollection<DerivedElementDataType, Range1D<int>> derived_element_data_collection
   ) const {

    int my_rank = index.value;
    int min_rank = index.min_value;

    //if (my_rank == min_rank) {
    //  std::cout << "\nComputing derived element data" << std::endl << std::flush;
    //}

    auto macroscale_mesh_handle = macroscale_mesh_collection[index].local_access();
    auto data_manager_handle = data_manager_collection[index].local_access();
    auto derived_element_data_handle = derived_element_data_collection[index].local_access();

    ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();

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
      if (derived_element_data_handle.get_value().find(block_id) == derived_element_data_handle.get_value().end()) {
        derived_element_data_handle.get_reference()[block_id] = std::vector< std::vector<double> >();
      }
      int num_elem_in_block = macroscale_mesh_handle.get_value().GetNumElementsInBlock(block_id);
      const int * const elem_conn = macroscale_mesh_handle.get_value().GetConnectivity(block_id);
      std::vector<double> const & elem_data_np1 = model_data.GetElementDataNew(block_id);
      block.ComputeDerivedElementData(reference_coordinate,
                                      displacement,
                                      num_elem_in_block,
                                      elem_conn,
                                      elem_data_labels.at(block_id).size(),
                                      elem_data_np1,
                                      derived_elem_data_labels.at(block_id).size(),
                                      derived_element_data_handle.get_reference().at(block_id));
    }

    //if (my_rank == min_rank) {
    //  std::cout << "  complete" << std::endl << std::flush;
    //}
  }

   void ExodusWriteStep::operator()(
     Index1D<int> index,
     AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
     AccessHandleCollection<ExodusOutput, Range1D<int>> exodus_output_collection,
     AccessHandleCollection<DerivedElementDataType, Range1D<int>> derived_element_data_collection,
     double time_current
   ) const {

     int num_ranks = index.max_value + 1;
     int my_rank = index.value;
     int min_rank = index.min_value;

     //if (my_rank == min_rank) {
     //  std::cout << "\nWriting to exodus output" << std::endl << std::flush;
     //}

     auto data_manager_handle = data_manager_collection[index].local_access();
     auto exodus_output_handle = exodus_output_collection[index].local_access();
     auto derived_element_data_handle = derived_element_data_collection[index].local_access();

     ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();

     int num_global_data = 0;
     std::vector<double> global_data(num_global_data);
     std::vector< std::vector<double> > node_data_for_output;
     model_data.GetNodeDataForOutput(node_data_for_output);
     std::map<int, std::vector<std::string> > const & elem_data_labels_for_output = model_data.GetElementDataLabelsForOutput();
     std::map<int, std::vector< std::vector<double> > > elem_data_for_output;
     model_data.GetElementDataForOutput(elem_data_for_output);
     std::map<int, std::vector<std::string> > const & derived_elem_data_labels = model_data.GetDerivedElementDataLabelsForOutput();
     exodus_output_handle.get_reference().WriteStep(time_current,
						    global_data,
						    node_data_for_output,
						    elem_data_labels_for_output,
						    elem_data_for_output,
						    derived_elem_data_labels,
						    derived_element_data_handle.get_value());

     //if (my_rank == min_rank) {
     //  std::cout << "  complete" << std::endl << std::flush;
     //}
   }

  void Initialize::operator()(
    Index1D<int> index,
    AccessHandleCollection<GenesisMesh, Range1D<int>> macroscale_mesh_collection,
    AccessHandleCollection<GenesisMesh, Range1D<int>> rve_mesh_collection,
    AccessHandleCollection<BoundaryConditionManager, Range1D<int>> boundary_condition_manager_collection,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> boundary_node_global_ids_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> boundary_ranks_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> my_global_node_ids_collection,
    AccessHandleCollection<ExodusOutput, Range1D<int>> exodus_output_collection,
    Parser parser
  ) const {

    int num_ranks = index.max_value + 1;
    int my_rank = index.value;
    int min_rank = index.min_value;

    ReadGenesisFiles()(index, macroscale_mesh_collection, rve_mesh_collection, parser);

    InitializeDataManager()(index, macroscale_mesh_collection, rve_mesh_collection, data_manager_collection, parser);

    create_work([=]{
	IdentifyGloballySharedNodes()(index,
				      macroscale_mesh_collection,
				      data_manager_collection,
				      boundary_node_global_ids_collection,
				      boundary_ranks_collection,
				      my_global_node_ids_collection);
      });

    create_work([=]{
	InitializeBoundaryConditionManager()(index,
					     macroscale_mesh_collection,
					     boundary_condition_manager_collection,
					     parser);
      });

    create_work([=]{
	InitializeExodusOutput()(index,
				 macroscale_mesh_collection,
				 data_manager_collection,
				 exodus_output_collection,
				 parser);
      });
  }

  void ExplicitTimeStep::operator()(
    Index1D<int> index,
    AccessHandleCollection<GenesisMesh, Range1D<int>> macroscale_mesh_collection,
    AccessHandleCollection<BoundaryConditionManager, Range1D<int>> boundary_condition_manager_collection,
    AccessHandleCollection<DataManager, Range1D<int>> data_manager_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> boundary_node_global_ids_collection,
    AccessHandleCollection<std::vector<int>, Range1D<int>> boundary_ranks_collection,
    AccessHandleCollection<std::vector<double>, Range1D<int>> global_reduction_buffer_collection,
    int step,
    int num_steps,
    double time_current,
    double time_previous,
    ProgressBarFlag progress_bar_flag,
    bool is_output_step
  ) const {

    int num_ranks = index.max_value + 1;
    int my_rank = index.value;
    int min_rank = index.min_value;

    if (my_rank == min_rank) {
      PrintProgress(progress_bar_flag, step, num_steps);
    }

    auto mesh_handle = macroscale_mesh_collection[index].local_access();
    auto boundary_condition_manager_handle = boundary_condition_manager_collection[index].local_access();
    auto data_manager_handle = data_manager_collection[index].local_access();

    ModelData& model_data = data_manager_handle.get_reference().GetMacroScaleData();

    double delta_time = time_current - time_previous;
    double half_delta_time = 0.5 * delta_time;

    int dim = mesh_handle.get_value().GetDim();
    int num_nodes = mesh_handle.get_value().GetNumNodes();
    unsigned int num_unknowns = num_nodes * dim;

    int lumped_mass_field_id = model_data.GetFieldId("lumped_mass");
    int reference_coordinate_field_id = model_data.GetFieldId("reference_coordinate");
    int displacement_field_id = model_data.GetFieldId("displacement");
    int velocity_field_id = model_data.GetFieldId("velocity");
    int acceleration_field_id = model_data.GetFieldId("acceleration");
    int internal_force_field_id = model_data.GetFieldId("internal_force");
    int external_force_field_id = model_data.GetFieldId("external_force");
    double* lumped_mass = model_data.GetNodeData(lumped_mass_field_id);
    double* reference_coordinate = model_data.GetNodeData(reference_coordinate_field_id);
    double* displacement = model_data.GetNodeData(displacement_field_id);
    double* velocity = model_data.GetNodeData(velocity_field_id);
    double* acceleration = model_data.GetNodeData(acceleration_field_id);
    double* internal_force = model_data.GetNodeData(internal_force_field_id);
    double* external_force = model_data.GetNodeData(external_force_field_id);

    std::map<int, Block> & blocks = model_data.GetBlocks();
    std::map<int, std::vector<std::string> > const & elem_data_labels = model_data.GetElementDataLabels();

    // For explicit dynamics, the macroscale model is never treated as an RVE
    // so rve_macroscale_deformation_gradient will always be the identity matrix
    std::vector<double> rve_macroscale_deformation_gradient(dim*dim, 0.0);
    for (int i=0 ; i<dim ; i++) {
      rve_macroscale_deformation_gradient[i] = 1.0;
    }

    // V^{n+1/2} = V^{n} + (dt/2) * A^{n}
    for (int i=0 ; i<num_unknowns ; ++i) {
      velocity[i] += half_delta_time * acceleration[i];
    }

    ApplyKinematicBC()(index, boundary_condition_manager_collection, data_manager_collection, time_current, time_previous);

    // Evaluate external body forces
    for (int i=0 ; i<num_unknowns ; ++i) {
      external_force[i] = 0.0;
    }

    // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
    for (int i=0 ; i<num_unknowns ; ++i) {
      displacement[i] += delta_time * velocity[i];
    }

    // Evaluate the internal force
    for (int i=0 ; i<num_unknowns ; ++i) {
      internal_force[i] = 0.0;
    }

    for (auto const & entry : blocks) {
      int block_id = entry.first;
      Block const & block = entry.second;
      int num_elem_in_block = mesh_handle.get_value().GetNumElementsInBlock(block_id);
      std::vector<int> const & elem_global_ids = mesh_handle.get_value().GetElementGlobalIdsInBlock(block_id);
      int const * elem_conn = mesh_handle.get_value().GetConnectivity(block_id);
      std::vector<double> const & elem_data_n = model_data.GetElementDataOld(block_id);
      std::vector<double> & elem_data_np1 = model_data.GetElementDataNew(block_id);

      block.ComputeInternalForce(reference_coordinate,
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
                                 data_manager_handle.get_reference(),
                                 is_output_step);
    }

    // Check for NaNs and Infs
    CheckVectorSanity(num_unknowns, internal_force, "internal force");

    DistributedVectorReduction()(index,
				 data_manager_collection,
				 boundary_node_global_ids_collection,
				 boundary_ranks_collection,
				 global_reduction_buffer_collection,
				 INTERNAL_FORCE,
				 step+1);

    create_work([data_manager_handle, num_unknowns, half_delta_time]{

        // the data_manager_handle is valid within this create_work, but the pointers from
        // above are not

        ModelData& model_data_cw = data_manager_handle.get_reference().GetMacroScaleData();
        int lumped_mass_field_id_cw = model_data_cw.GetFieldId("lumped_mass");
        int velocity_field_id_cw = model_data_cw.GetFieldId("velocity");
        int acceleration_field_id_cw = model_data_cw.GetFieldId("acceleration");
        int internal_force_field_id_cw = model_data_cw.GetFieldId("internal_force");
        int external_force_field_id_cw = model_data_cw.GetFieldId("external_force");
        double* lumped_mass_cw = model_data_cw.GetNodeData(lumped_mass_field_id_cw);
        double* velocity_cw = model_data_cw.GetNodeData(velocity_field_id_cw);
        double* acceleration_cw = model_data_cw.GetNodeData(acceleration_field_id_cw);
        double* internal_force_cw = model_data_cw.GetNodeData(internal_force_field_id_cw);
        double* external_force_cw = model_data_cw.GetNodeData(external_force_field_id_cw);

        // A^{n+1} = M^{-1} ( F^{n} + b^{n} )
        for (int i=0 ; i<num_unknowns ; ++i) {
          acceleration_cw[i] = (1.0/lumped_mass_cw[i/3]) * (internal_force_cw[i] + external_force_cw[i]);
        }

        // V^{n+1} = V^{n+1/2} + (dt/2)*A^{n+1}
        for (int i=0 ; i<num_unknowns ; ++i) {
          velocity_cw[i] += half_delta_time * acceleration_cw[i];
        }

      });
  }

} // namespace nimble
