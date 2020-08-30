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

#include "nimble_exodus_output.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <ctime>
#include <cstring>
#include <set>

#ifdef NIMBLE_HAVE_EXODUS
  #include "exodusII.h"
#endif

namespace nimble {

  void ExodusOutput::Initialize(std::string const & filename, GenesisMesh const & genesis_mesh) {
    filename_ = filename;
    dim_ = genesis_mesh.GetDim();
    num_nodes_ = static_cast<int>( genesis_mesh.GetNumNodes() );
    num_elements_ = static_cast<int>( genesis_mesh.GetNumElements() );
    num_blocks_ = static_cast<int>( genesis_mesh.GetNumBlocks() );
    num_global_blocks_ = static_cast<int>( genesis_mesh.GetNumGlobalBlocks() );
    all_block_ids_ = genesis_mesh.GetAllBlockIds();
    block_ids_ = genesis_mesh.GetBlockIds();
    num_node_sets_ = genesis_mesh.GetNumNodeSets();
  }

  void ExodusOutput::InitializeDatabase(GenesisMesh const & genesis_mesh,
                                        std::vector<std::string> const & global_data_names,
                                        std::vector<std::string> const & node_data_names,
                                        std::map< int, std::vector<std::string> > const & elem_data_names,
                                        std::map< int, std::vector<std::string> > const & derived_elem_data_names) {

#ifndef NIMBLE_HAVE_EXODUS
    InitializeDatabaseTextFile(genesis_mesh, global_data_names, node_data_names, elem_data_names, derived_elem_data_names);
#else

    int num_global_data = static_cast<int>(global_data_names.size());
    int num_node_data = static_cast<int>(node_data_names.size());

    // Initialize exodus database; Overwrite any existing file with this name
    int exodus_file_id = ex_create(filename_.c_str(), EX_CLOBBER, &CPU_word_size_, &IO_word_size_);
    if (exodus_file_id < 0) ReportExodusError(exodus_file_id, "InitializeDatabase", "ex_create");

    // Write the Quality Assurance (QA) record
    int retval = ex_put_init(exodus_file_id, "NimbleSM", dim_, num_nodes_, num_elements_, num_global_blocks_, num_node_sets_, num_side_sets_);
    if (retval != 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_init");
    WriteQARecord(exodus_file_id);

    // Write nodal coordinate names to database
    const char *coord_names[3] = {"x", "y", "z"};
    retval = ex_put_coord_names(exodus_file_id,const_cast<char**>(coord_names));
    if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_coord_names");

    // Write the coordinates
    retval = ex_put_coord(exodus_file_id, genesis_mesh.GetCoordinatesX(), genesis_mesh.GetCoordinatesY(), genesis_mesh.GetCoordinatesZ());
    if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_coord");

    // Get the list of global IDs for element blocks
    std::set<int> gid_list;
    for(int i=0 ; i < num_blocks_ ; ++i)
      gid_list.insert(block_ids_[i]);

    // Write element block parameters
    //
    // Note: EPU expects to see all the global blocks in each file.
    // It is allowed to have 0 elements in a block in a processor.
    //
    for(int i=0 ; i < num_global_blocks_ ; ++i) {
      const auto &gblock_list = genesis_mesh.GetAllBlockIds();
      int gid = gblock_list[i];
      int num_elements_in_block = 0;
      int num_nodes_per_elem = 0;
      std::string elem_type;
      if (gid_list.count(gid) > 0) {
        num_elements_in_block = genesis_mesh.GetNumElementsInBlock(gid);
        num_nodes_per_elem = genesis_mesh.GetNumNodesPerElement(gid);
        elem_type = genesis_mesh.GetElementType(gid);
      }
      retval = ex_put_block(exodus_file_id, EX_ELEM_BLOCK, gid, elem_type.c_str(), num_elements_in_block, num_nodes_per_elem, 0, 0, 0);
      if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_block");
    }

    // Write the block names
    char **block_names = new char*[num_global_blocks_];
    for(int i=0 ; i<num_global_blocks_ ; ++i){
      block_names[i] = new char[MAX_STR_LENGTH+1];
      const int id = all_block_ids_[i];
      const std::string block_name = genesis_mesh.GetBlockName(id);
      strcpy(block_names[i], block_name.c_str());
    }
    retval = ex_put_names(exodus_file_id, EX_ELEM_BLOCK, block_names);
    if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_names EX_ELEM_BLOCK");

    // Write element connectivity
    for(int i=0 ; i<num_blocks_ ; ++i){
      int id = block_ids_[i];
      int num_elements_in_block = genesis_mesh.GetNumElementsInBlock(id);
      if (num_elements_in_block > 0) {
        const int* conn = genesis_mesh.GetConnectivity(id);
        // Switch from 0-based indexing to 1-based indexing
        int num_node_in_elem = genesis_mesh.GetNumNodesPerElement(id);
        std::vector<int> exodus_conn(num_elements_in_block*num_node_in_elem);
        for (int j=0 ; j<num_elements_in_block*num_node_in_elem ; j++) {
          exodus_conn[j] = conn[j] + 1;
        }
        retval = ex_put_conn(exodus_file_id, EX_ELEM_BLOCK, id, &exodus_conn[0], nullptr, nullptr);
        if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_conn");
      }
    }

    // Write global node number map (global node IDs)
    const int* node_global_ids = genesis_mesh.GetNodeGlobalIds();
    // Switch to 1-based indexing
    std::vector<int> temp_node_global_ids(num_nodes_);
    for (int i = 0 ; i < num_nodes_ ; ++i) {
      temp_node_global_ids[i] = node_global_ids[i] + 1;
    }
    retval = ex_put_id_map(exodus_file_id, EX_NODE_MAP, &temp_node_global_ids[0]);
    if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_id_map");

    // Write global element number map (global element IDs)
    const int* elem_global_ids = genesis_mesh.GetElementGlobalIds();
    // Switch to 1-based indexing
    std::vector<int> temp_elem_global_ids(num_elements_);
    for (int i = 0 ; i < num_elements_ ; ++i) {
      temp_elem_global_ids[i] = elem_global_ids[i] + 1;
    }
    retval = ex_put_id_map(exodus_file_id, EX_ELEM_MAP, &temp_elem_global_ids[0]);
    if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_id_map");

    // Write node sets
    if(num_node_sets_ > 0) {
      std::vector<int> node_set_ids = genesis_mesh.GetNodeSetIds();
      std::map<int, std::string> node_set_names = genesis_mesh.GetNodeSetNames();
      std::map<int, std::vector<int> > node_sets = genesis_mesh.GetNodeSets();
      std::map<int, std::vector<double> > distribution_factors = genesis_mesh.GetDistributionFactors();
      std::vector<int> num_nodes_per_set(num_node_sets_);
      std::vector<int> num_distribution_factors_per_set(num_node_sets_);
      int total_num_node(0);
      int total_num_dist(0);
      for (int i=0 ; i<num_node_sets_; i++) {
        int id = node_set_ids[i];
        int nn = static_cast<int>( node_sets[id].size() );
        int ndf = static_cast<int>( distribution_factors[id].size() );
        num_nodes_per_set[i] = nn;
        num_distribution_factors_per_set[i] = ndf;
        total_num_node += nn;
        total_num_dist += ndf;
      }
      std::vector<int> node_sets_node_index(num_node_sets_);
      std::vector<int> node_sets_dist_index(num_node_sets_);
      std::vector<int> node_sets_node_list(total_num_node);
      std::vector<double> node_sets_dist_fact(total_num_dist);
      int node_index(0);
      int dist_index(0);
      for (int i=0 ; i<num_node_sets_; i++) {
        int id = node_set_ids[i];
        node_sets_node_index[i] = node_index;
        node_sets_dist_index[i] = dist_index;
        for (const auto &n_set : node_sets[id]) {
          // Switch from 0-based indexing to 1-based indexing
          node_sets_node_list[node_index++] = n_set + 1;
        }
        for (const auto &j_factor : distribution_factors[id]) {
          node_sets_dist_fact[dist_index++] = j_factor;
        }
      }

      ex_set_specs set_specs;
      set_specs.sets_ids = node_set_ids.data();
      set_specs.num_entries_per_set = num_nodes_per_set.data();
      set_specs.num_dist_per_set = num_distribution_factors_per_set.data();
      set_specs.sets_entry_index = node_sets_node_index.data();
      set_specs.sets_dist_index = node_sets_dist_index.data();
      set_specs.sets_entry_list = node_sets_node_list.data();
      set_specs.sets_extra_list = nullptr;
      set_specs.sets_dist_fact = node_sets_dist_fact.data();
      retval = ex_put_concat_sets(exodus_file_id, EX_NODE_SET, &set_specs);
      if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_concat_sets");
    }

    // Write global data info
    char **global_var_names = nullptr;
    if (num_global_data > 0) {
      global_var_names = new char*[num_global_data];
      for (int i=0 ; i<num_global_data ; i++) {
        global_var_names[i] = new char[MAX_STR_LENGTH+1];
        strcpy(global_var_names[i], global_data_names[i].c_str());
      }
      retval = ex_put_variable_param(exodus_file_id, EX_GLOBAL, num_global_data);
      if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_variable_param");
      retval = ex_put_variable_names(exodus_file_id, EX_GLOBAL, num_global_data, global_var_names);
      if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_variable_names");
    }

    // Write node data info
    char **node_var_names = nullptr;
    if (num_node_data > 0 && num_nodes_ > 0) {
      node_var_names = new char*[num_node_data];
      for (int i=0 ; i<num_node_data ; i++) {
        node_var_names[i] = new char[MAX_STR_LENGTH+1];
        strcpy(node_var_names[i], node_data_names[i].c_str());
      }
      retval = ex_put_variable_param(exodus_file_id, EX_NODAL, num_node_data);
      if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_variable_param");
      retval = ex_put_variable_names(exodus_file_id, EX_NODAL, num_node_data, node_var_names);
      if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_variable_names");
    }

    // Write element data info
    std::set<std::string> unique_elem_var_names;
    for (int i=0 ; i<num_blocks_ ; ++i){
      int id = block_ids_[i];
      for (const auto & jname : elem_data_names.at(id)) {
        unique_elem_var_names.insert(jname);
      }
      for (const auto & jname : derived_elem_data_names.at(id)) {
        unique_elem_var_names.insert(jname);
      }
    }

    // Create map from data name to exodus element data index
    std::vector<std::string> elem_var_names;
    for (std::set<std::string>::const_iterator it = unique_elem_var_names.begin() ; it != unique_elem_var_names.end() ; it++) {
      elem_var_names.push_back(*it);
      elem_data_index_[*it] = static_cast<int>(elem_var_names.size());
    }

    int num_element_vars = elem_var_names.size();
    char **element_var_names = nullptr;
    if (num_element_vars > 0) {
      element_var_names = new char*[num_element_vars];
      for (int i=0 ; i<num_element_vars ; i++) {
        element_var_names[i] = new char[MAX_STR_LENGTH+1];
        strcpy(element_var_names[i], elem_var_names[i].c_str());
      }
      retval = ex_put_variable_param(exodus_file_id, EX_ELEM_BLOCK, num_element_vars);
      if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_variable_param");
      retval = ex_put_variable_names(exodus_file_id, EX_ELEM_BLOCK, num_element_vars, element_var_names);
      if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_put_variable_names");
    }

    // Close exodus file
    retval = ex_update(exodus_file_id);
    if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_update");
    retval = ex_close(exodus_file_id);
    if (retval!= 0) ReportExodusError(retval, "InitializeDatabase", "ex_close");

    // Clean up
    if (block_names != nullptr) {
      for (int i = num_blocks_; i>0 ; i--) {
        delete[] block_names[i-1];
      }
      delete[] block_names;
    }
    if (global_var_names != nullptr) {
      for (int i = num_global_data; i>0 ; i--) {
        delete[] global_var_names[i-1];
      }
      delete[] global_var_names;
    }
    if (node_var_names != nullptr) {
      for (int i = num_node_data; i>0 ; i--) {
        delete[] node_var_names[i-1];
      }
      delete[] node_var_names;
    }
    if (element_var_names != nullptr) {
      for (int i = num_element_vars; i>0 ; i--) {
        delete[] element_var_names[i-1];
      }
      delete[] element_var_names;
    }
#endif
  }

  void
  ExodusOutput::WriteStep(double time,
			  std::vector<double> const & global_data,
			  std::vector< std::vector<double> > const & node_data,
        std::map< int, std::vector< std::string > > const & elem_data_names,
        std::map< int, std::vector< std::vector<double> > > const & elem_data,
        std::map< int, std::vector< std::string > > const & derived_elem_data_names,
        std::map< int, std::vector< std::vector<double> > > const & derived_elem_data) {

    exodus_write_count_ += 1;

#ifndef NIMBLE_HAVE_EXODUS
    WriteStepTextFile(time,
                      global_data,
                      node_data,
                      elem_data_names,
                      elem_data,
                      derived_elem_data_names,
                      derived_elem_data);
#else

    float exodus_version;
    int exodus_file_id = ex_open(filename_.c_str(), EX_WRITE, &CPU_word_size_, &IO_word_size_, &exodus_version);
    if (exodus_file_id < 0) ReportExodusError(exodus_file_id, "WriteStep", "ex_open");

    // Write time value
    int retval = ex_put_time(exodus_file_id, exodus_write_count_, &time);
    if (retval!= 0) ReportExodusError(retval, "WriteStep", "ex_put_time");

    // Write global data
    int num_global_vars = static_cast<int>(global_data.size());
    if (num_global_vars > 0) {
      retval = ex_put_var(exodus_file_id, exodus_write_count_, EX_GLOBAL, 1, 0, num_global_vars, &global_data[0]);
      if (retval!= 0) ReportExodusError(retval, "WriteStep", "ex_put_var");
    }

    // Write node data
    if (num_nodes_ > 0) {
      for(unsigned int i=0 ; i<node_data.size() ; ++i) {
        int variable_index = static_cast<int>(i+1);
        retval = ex_put_var(exodus_file_id, exodus_write_count_, EX_NODAL, variable_index, 1, num_nodes_, &node_data[i][0]);
        if (retval!= 0) ReportExodusError(retval, "WriteStep", "ex_put_var");
      }
    }

    // Write element data
    for(int i=0 ; i<num_blocks_ ; ++i){
      int block_id = block_ids_[i];
      int num_elem_data = elem_data_names.at(block_id).size();
      for(int j=0 ; j<num_elem_data ; ++j) {
        std::string data_name = elem_data_names.at(block_id).at(j);
        int variable_index = elem_data_index_.at(data_name);
        int block_num_elem = static_cast<int>(elem_data.at(block_id).at(j).size());
        retval = ex_put_var(exodus_file_id, exodus_write_count_, EX_ELEM_BLOCK, variable_index, block_id, block_num_elem, &elem_data.at(block_id).at(j)[0]);
        if (retval!= 0) ReportExodusError(retval, "WriteStep", "ex_put_var");
      }
      int num_derived_elem_data = derived_elem_data_names.at(block_id).size();
      for(int j=0 ; j<num_derived_elem_data ; ++j) {
        std::string data_name = derived_elem_data_names.at(block_id).at(j);
        int variable_index = elem_data_index_.at(data_name);
        int block_num_elem = static_cast<int>(derived_elem_data.at(block_id).at(j).size());
        retval = ex_put_var(exodus_file_id, exodus_write_count_, EX_ELEM_BLOCK, variable_index, block_id, block_num_elem, &derived_elem_data.at(block_id).at(j)[0]);
        if (retval!= 0) ReportExodusError(retval, "WriteStep", "ex_put_var");
      }
    }

    retval = ex_update(exodus_file_id);
    if (retval!= 0) ReportExodusError(retval, "WriteStep", "ex_update");
    retval = ex_close(exodus_file_id);
    if (retval!= 0) ReportExodusError(retval, "WriteStep", "ex_close");
#endif
  }

  void ExodusOutput::InitializeDatabaseTextFile(GenesisMesh const & genesis_mesh,
                                                std::vector<std::string> const & global_data_names,
                                                std::vector<std::string> const & node_data_names,
                                                std::map< int, std::vector<std::string> > const & elem_data_names,
                                                std::map< int, std::vector<std::string> > const & derived_elem_data_names) {

    std::string text_filename = filename_ + ".txt";
    std::ofstream output_file;
    output_file.open(text_filename);
    output_file << "number_of_nodes " << num_nodes_ << std::endl;
    output_file << "number_of_elements " << num_elements_ << std::endl;
    output_file << "number_of_blocks " << num_blocks_ << std::endl;
    output_file << "number_of_node_sets " << num_node_sets_ << std::endl;
    output_file << "node_coordinates" << std::endl;

    // Write the coordinates
    const int* node_global_ids = genesis_mesh.GetNodeGlobalIds();
    const double * const x_coord = genesis_mesh.GetCoordinatesX();
    const double * const y_coord = genesis_mesh.GetCoordinatesY();
    const double * const z_coord = genesis_mesh.GetCoordinatesZ();
    for (int i=0 ; i<num_nodes_ ; ++i) {
      // switch to 1-based indexing
      int global_id = node_global_ids[i] + 1;
      output_file << std::setprecision(16) << global_id << " " <<
        x_coord[i] << " " << y_coord[i] << " " << z_coord[i] << std::endl;
    }

    // Write the element blocks
    for (int i_block=0 ; i_block<num_blocks_; ++i_block) {
      int block_id = block_ids_[i_block];
      std::string block_name = genesis_mesh.GetBlockName(block_id);
      int num_elem_in_block = genesis_mesh.GetNumElementsInBlock(block_id);
      int num_nodes_per_elem = genesis_mesh.GetNumNodesPerElement(block_id);
      output_file << "element_block\n" << block_id << " " << block_name << " "
                  << num_elem_in_block << " " << num_nodes_per_elem << std::endl;
      std::vector<int> const & global_elem_ids = genesis_mesh.GetElementGlobalIdsInBlock(block_id);
      const int * const elem_conn = genesis_mesh.GetConnectivity(block_id);
      for (int i_elem=0 ; i_elem<num_elem_in_block ; i_elem++) {
        // switch to 1-based indexing
        int global_id = global_elem_ids[i_elem] + 1;
        output_file << global_id;
        for (int i_node=0 ; i_node<num_nodes_per_elem ; ++i_node) {
          // switch to 1-based indexing
          int node_id = elem_conn[i_elem*num_nodes_per_elem + i_node] + 1;
          output_file << " " << node_id;
        }
        output_file << std::endl;
      }
    }

    // Write node sets
    std::vector<int> node_set_ids = genesis_mesh.GetNodeSetIds();
    std::map<int, std::string> node_set_names = genesis_mesh.GetNodeSetNames();
    std::map<int, std::vector<int> > node_sets = genesis_mesh.GetNodeSets();
    for (auto & node_set_id : node_set_ids) {
      std::vector<int> const & node_set = node_sets[node_set_id];
      output_file << "nodeset" << std::endl;
      output_file << node_set_id << " " << node_set_names[node_set_id] << " " << node_set.size() << std::endl;
      for (auto & node_id : node_set) {
        // switch to 1-based indexing
        output_file << node_id + 1 << std::endl;
      }
    }

    // Write information on global data
    output_file << "global_data" << std::endl;
    output_file << global_data_names.size() << std::endl;
    for (auto & name : global_data_names) {
      output_file << name << std::endl;
    }

    // Write information on node data
    output_file << "node_data" << std::endl;
    output_file << node_data_names.size() << std::endl;
    for (auto & name : node_data_names) {
      output_file << name << std::endl;
    }

    // Write information on element data
    std::set<std::string> unique_elem_var_names;
    for (int i=0 ; i<num_blocks_ ; ++i){
      int id = block_ids_[i];
      for (const auto &my_name : elem_data_names.at(id))
        unique_elem_var_names.insert(my_name);
      for (const auto &my_name : derived_elem_data_names.at(id))
          unique_elem_var_names.insert(my_name);
    }
    // Create map from data name to exodus element data index
    std::vector<std::string> elem_var_names;
    for (std::set<std::string>::const_iterator it = unique_elem_var_names.begin() ; it != unique_elem_var_names.end() ; it++) {
      elem_var_names.push_back(*it);
      elem_data_index_[*it] = static_cast<int>(elem_var_names.size());
    }
    int num_element_vars = elem_var_names.size();
    output_file << "element_data" << std::endl;
    output_file << num_element_vars << std::endl;
    for (auto const & name : elem_var_names) {
      output_file << name << std::endl;
    }

    output_file.close();
  }

  void ExodusOutput::WriteStepTextFile(double time,
                                       std::vector<double> const & global_data,
                                       std::vector< std::vector<double> > const & node_data,
                                       std::map< int, std::vector< std::string > > const & elem_data_names,
                                       std::map< int, std::vector< std::vector<double> > > const & elem_data,
                                       std::map< int, std::vector< std::string > > const & derived_elem_data_names,
                                       std::map< int, std::vector< std::vector<double> > > const & derived_elem_data) {

    std::string text_filename = filename_ + ".txt";
    std::ofstream output_file;
    output_file.open(text_filename, std::ios_base::app);

    // Write the time
    output_file << "time" << std::endl;
    output_file << std::setprecision(16) << time << std::endl;

    // Global data
    output_file << "global_data_values" << std::endl;
    for (auto const & datum : global_data) {
      output_file << datum << std::endl;
    }

    // Node data
    output_file << "node_data_values" << std::endl;
    for (auto const & data : node_data) {
      for (auto const & datum : data) {
        output_file << std::setprecision(16) << datum << std::endl;
      }
    }

    // Element data
    output_file << "element_data_values" << std::endl;
    for (int i=0 ; i<num_blocks_ ; ++i){
      int block_id = block_ids_[i];
      int num_elem_data = elem_data_names.at(block_id).size();
      for (int j=0 ; j<num_elem_data ; ++j) {
        std::string data_name = elem_data_names.at(block_id).at(j);
        int variable_index = elem_data_index_.at(data_name);
        int block_num_elem = static_cast<int>(elem_data.at(block_id).at(j).size());
        output_file << data_name << " " << variable_index << " " << block_id << " " << block_num_elem << std::endl;
        std::vector<double> const & data = elem_data.at(block_id).at(j);
        for (auto const & datum : data) {
          output_file << std::setprecision(16) << datum << std::endl;
        }
      }
      int num_derived_elem_data = derived_elem_data_names.at(block_id).size();
      for(int j=0 ; j<num_derived_elem_data ; ++j) {
        std::string data_name = derived_elem_data_names.at(block_id).at(j);
        int variable_index = elem_data_index_.at(data_name);
        int block_num_elem = static_cast<int>(derived_elem_data.at(block_id).at(j).size());
        output_file << data_name << " " << variable_index << " " << block_id << " " << block_num_elem << std::endl;
        std::vector<double> const & data = derived_elem_data.at(block_id).at(j);
        for (auto const & datum : data) {
          output_file << std::setprecision(16) << datum << std::endl;
        }
      }
    }

    output_file.close();
  }

  void
  ExodusOutput::WriteQARecord(int exodus_file_id) {

#ifdef NIMBLE_HAVE_EXODUS
    // Get the current system date and time
    // DJL:  NOT THREAD SAFE, SO JUST PRINT NOTHING FOR THE DATE AND TIME
    // std::chrono::time_point<std::chrono::system_clock> time_now = std::chrono::system_clock::now();
    // std::time_t time_now_t = std::chrono::system_clock::to_time_t(time_now);
    std::stringstream datestream, timestream;
    // datestream << std::put_time(std::localtime(&time_now_t), "%F");
    // timestream << std::put_time(std::localtime(&time_now_t), "%T");

    // Quality assurance (QA) data
    std::string qa_name_string, qa_descriptor_string, qa_date_string, qa_time_string;
    qa_name_string = "NimbleSM";
    qa_descriptor_string = "unknown";
    qa_date_string = datestream.str();
    qa_time_string = timestream.str();

    // Copy to required c-style arrays
    int num_qa_records = 1;
    char* qa_records[1][4];
    char qa_name[MAX_STR_LENGTH+1], qa_descriptor[MAX_STR_LENGTH+1], qa_date[MAX_STR_LENGTH+1], qa_time[MAX_STR_LENGTH+1];
    qa_records[0][0] = qa_name;
    qa_records[0][1] = qa_descriptor;
    qa_records[0][2] = qa_date;
    qa_records[0][3] = qa_time;

    strcpy(qa_name, qa_name_string.c_str());
    strcpy(qa_descriptor, qa_descriptor_string.c_str());
    strcpy(qa_date, qa_date_string.c_str());
    strcpy(qa_time, qa_time_string.c_str());

    int retval = ex_put_qa(exodus_file_id, num_qa_records, qa_records);
    if (retval!= 0) ReportExodusError(retval, "WriteQARecord", "ex_put_qa");
#endif
  }

  void
  ExodusOutput::ReportExodusError(int error_code,
				  const char *method_name,
				  const char *exodus_method_name) {
    std::stringstream ss;
    if (error_code < 0) {
      ss << "**Error ExodusOutput::" << method_name << "(), Error code: " << error_code << " (" << exodus_method_name << ")\n";
      throw std::logic_error(ss.str());
    }
    else {
      ss << "**Warning ExodusOutput::" << method_name << "(), Warning code: " << error_code << " (" << exodus_method_name << ")\n";
      std::cout << ss.str() << std::endl;
    }
  }

}
