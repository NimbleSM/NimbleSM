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

#include "nimble_genesis_mesh.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "nimble_macros.h"

#ifdef NIMBLE_HAVE_EXODUS
#include "exodusII.h"
#endif

namespace nimble {

void
GenesisMesh::ReadFile(std::string file_name)
{
#ifndef NIMBLE_HAVE_EXODUS
  ReadTextFile(file_name);
#else

  file_name_ = file_name;

  if (!IsValid()) { return; }

  // Open the genesis file
  int   word_size    = sizeof(double);
  int   io_word_size = 0;
  float exodus_version;
  int   exodus_file_id = ex_open(file_name_.c_str(), EX_READ, &word_size, &io_word_size, &exodus_version);
  if (exodus_file_id < 0) {
    std::cout << "\n****Error: unable to open file: " << file_name_.c_str() << "\n" << std::endl;
    ReportExodusError(exodus_file_id, "GenesisMesh::ReadFile", "ex_open");
  }

  // Read the initialization parameters
  int  num_nodes, num_elem, num_blocks, num_node_sets, num_side_sets;
  char title[MAX_LINE_LENGTH + 1];
  int  retval =
      ex_get_init(exodus_file_id, title, &dim_, &num_nodes, &num_elem, &num_blocks, &num_node_sets, &num_side_sets);
  if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_init");

  // Node coordinates
  node_x_.resize(num_nodes);
  node_y_.resize(num_nodes);
  if (dim_ == 3) { node_z_.resize(num_nodes); }
  retval = ex_get_coord(exodus_file_id, &node_x_[0], &node_y_[0], &node_z_[0]);
  if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_coord");

  // Global node numbering (default map name)
  node_global_id_.resize(num_nodes);
  retval = ex_get_id_map(exodus_file_id, EX_NODE_MAP, &node_global_id_[0]);
  if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_id_map");
  for (int i = 0; i < num_nodes; ++i) {
    node_global_id_[i] -= 1;  // Switch from 1-based indexing to 0-based indexing
  }

  // Global element numbering (default map name)
  elem_global_id_.resize(num_elem);
  retval = ex_get_id_map(exodus_file_id, EX_ELEM_MAP, &elem_global_id_[0]);
  if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_id_map");
  for (int i = 0; i < num_elem; ++i) {
    elem_global_id_[i] -= 1;  // Switch from 1-based indexing to 0-based indexing
  }

  // Check for auxiliary node maps and element maps
  int num_node_maps, num_elem_maps;
  retval = ex_get_map_param(exodus_file_id, &num_node_maps, &num_elem_maps);
  if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_map_param");
  if (num_node_maps > 1 || num_elem_maps > 1) {
    throw std::logic_error(
        "GenesisMesh::ReadFile(), multiple auxiliary node/element maps not "
        "supported!");
  }

  // If maps are named "original_global_id_map", use them instead of the maps
  // obtained above The name "original_global_id_map" is used by seacas
  // utilities like decomp
  if (num_node_maps > 0) {
    char map_name[MAX_STR_LENGTH + 1];
    retval = ex_get_name(exodus_file_id, EX_NODE_MAP, 1, map_name);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_name");
    if (std::string(map_name) != std::string("original_global_id_map")) {
      throw std::logic_error("GenesisMesh::ReadFile(), unsupported auxiliary node map!");
    }
    retval = ex_get_num_map(exodus_file_id, EX_NODE_MAP, 1, &node_global_id_[0]);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_num_map");
    for (int i = 0; i < num_nodes; ++i) {
      node_global_id_[i] -= 1;  // Switch from 1-based indexing to 0-based indexing
    }
  }
  if (num_elem_maps > 0) {
    char map_name[MAX_STR_LENGTH + 1];
    retval = ex_get_name(exodus_file_id, EX_ELEM_MAP, 1, map_name);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_name");
    if (std::string(map_name) != std::string("original_global_id_map")) {
      throw std::logic_error("GenesisMesh::ReadFile(), unsupported auxiliary element map!");
    }
    retval = ex_get_num_map(exodus_file_id, EX_ELEM_MAP, 1, &elem_global_id_[0]);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_num_map");
    for (int i = 0; i < num_elem; ++i) {
      elem_global_id_[i] -= 1;  // Switch from 1-based indexing to 0-based indexing
    }
  }

  // Read the node sets
  if (num_node_sets > 0) {
    node_set_ids_.resize(num_node_sets);
    retval = ex_get_ids(exodus_file_id, EX_NODE_SET, &node_set_ids_[0]);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_ids");
  }
  for (int i = 0; i < num_node_sets; ++i) {
    int  id = node_set_ids_[i];
    char name[MAX_STR_LENGTH + 1];
    retval = ex_get_name(exodus_file_id, EX_NODE_SET, id, name);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_name");
    // If the node set name came back blank, create one that looks like
    // "nodelist_1", "nodelist_2", etc.
    std::string node_set_name(name);
    if (node_set_name.size() == 0) {
      std::stringstream ss;
      ss << "nodelist_" << id;
      node_set_name = ss.str();
    }
    node_set_names_[id] = node_set_name;

    int num_nodes_in_ns;
    int num_dist_factors_in_ns;
    retval = ex_get_set_param(exodus_file_id, EX_NODE_SET, id, &num_nodes_in_ns, &num_dist_factors_in_ns);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_set_param");
    node_sets_[id] = std::vector<int>();
    if (num_nodes_in_ns > 0 && num_dist_factors_in_ns == num_nodes_in_ns) {
      node_sets_[id] = std::vector<int>(num_nodes_in_ns);
      retval         = ex_get_set(exodus_file_id, EX_NODE_SET, id, &node_sets_[id][0], NULL);
      if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_set");
      // convert from 1-based indexing to 0-based indexing
      for (unsigned int j = 0; j < node_sets_[id].size(); j++) { node_sets_[id][j] -= 1; }
    }
    ns_distribution_factors_[id] = std::vector<double>();
    if (num_nodes_in_ns > 0) {
      ns_distribution_factors_[id] = std::vector<double>(num_dist_factors_in_ns);
      retval = ex_get_set_dist_fact(exodus_file_id, EX_NODE_SET, id, &ns_distribution_factors_[id][0]);
      if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_set_dist_fact");
    }
  }

  // Read the side sets
  if (num_side_sets > 0) {
    side_set_ids_.resize(num_side_sets);
    retval = ex_get_ids(exodus_file_id, EX_SIDE_SET, &side_set_ids_[0]);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_ids");
  }
  for (int i = 0; i < num_side_sets; ++i) {
    int  id = side_set_ids_[i];
    char name[MAX_STR_LENGTH + 1];
    retval = ex_get_name(exodus_file_id, EX_SIDE_SET, id, name);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_name");
    // If the side set name came back blank, create one that looks like
    // "sideset_1", "sideset_2", etc.
    std::string side_set_name(name);
    if (side_set_name.size() == 0) {
      std::stringstream ss;
      ss << "sideset_" << id;
      side_set_name = ss.str();
    }
    side_set_names_[id] = side_set_name;

    int num_nodes_in_ss;
    int num_dist_factors_in_ss;
    retval = ex_get_set_param(exodus_file_id, EX_SIDE_SET, id, &num_nodes_in_ss, &num_dist_factors_in_ss);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_set_param");
    side_sets_[id] = std::vector<int>();
    if (num_nodes_in_ss > 0 && num_dist_factors_in_ss == num_nodes_in_ss) {
      side_sets_[id] = std::vector<int>(num_nodes_in_ss);
      retval         = ex_get_set(exodus_file_id, EX_SIDE_SET, id, &side_sets_[id][0], NULL);
      if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_set");
      // convert from 1-based indexing to 0-based indexing
      for (unsigned int j = 0; j < side_sets_[id].size(); j++) { side_sets_[id][j] -= 1; }
    }
    ss_distribution_factors_[id] = std::vector<double>();
    if (num_nodes_in_ss > 0) {
      ss_distribution_factors_[id] = std::vector<double>(num_dist_factors_in_ss);
      retval = ex_get_set_dist_fact(exodus_file_id, EX_SIDE_SET, id, &ss_distribution_factors_[id][0]);
      if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_set_dist_fact");
    }
  }

  // Process the element blocks
  std::vector<int> all_block_ids(num_blocks);
  retval = ex_get_ids(exodus_file_id, EX_ELEM_BLOCK, &all_block_ids[0]);
  if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_ids");

  // Load information only for blocks with elements on this partition
  for (auto block_id : all_block_ids) {
    int  num_elem_this_block(0), num_nodes_per_elem(0), num_edges_per_elem(0), num_faces_per_elem(0), num_attributes(0);
    char elem_type[MAX_STR_LENGTH + 1];
    retval = ex_get_block(
        exodus_file_id,
        EX_ELEM_BLOCK,
        block_id,
        elem_type,
        &num_elem_this_block,
        &num_nodes_per_elem,
        &num_edges_per_elem,
        &num_faces_per_elem,
        &num_attributes);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_block");
    if (num_elem_this_block > 0) { block_ids_.push_back(block_id); }
    all_block_ids_.push_back(block_id);

    // Get the block name, if there is one
    char exodus_block_name[MAX_STR_LENGTH + 1];
    retval = ex_get_name(exodus_file_id, EX_ELEM_BLOCK, block_id, exodus_block_name);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_name");
    // If the block name came back blank, create one that looks like "block_1",
    // "block_2", etc.
    std::string block_name(exodus_block_name);
    if (block_name.size() == 0) {
      std::stringstream ss;
      ss << "block_" << block_id;
      block_name = ss.str();
    }
    if (num_elem_this_block > 0) { block_names_[block_id] = block_name; }
    all_block_names_[block_id] = block_name;
  }

  num_blocks = static_cast<int>(block_ids_.size());

  int elem_local_index = 0;

  for (int i_block = 0; i_block < num_blocks; i_block++) {
    int block_id = block_ids_.at(i_block);

    // Get the block parameters and the element connectivity
    int  num_elem_this_block(0), num_nodes_per_elem(0), num_edges_per_elem(0), num_faces_per_elem(0), num_attributes(0);
    char elem_type[MAX_STR_LENGTH + 1];
    retval = ex_get_block(
        exodus_file_id,
        EX_ELEM_BLOCK,
        block_id,
        elem_type,
        &num_elem_this_block,
        &num_nodes_per_elem,
        &num_edges_per_elem,
        &num_faces_per_elem,
        &num_attributes);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_block");
    block_num_nodes_per_elem_[block_id] = num_nodes_per_elem;

    // global element ids for this block
    block_elem_global_ids_[block_id] = std::vector<int>(num_elem_this_block);
    for (int i = 0; i < num_elem_this_block; i++) {
      block_elem_global_ids_.at(block_id).at(i) = elem_global_id_[elem_local_index++];
    }

    // element connectivity for this block
    block_elem_connectivity_[block_id] = std::vector<int>(num_elem_this_block * num_nodes_per_elem);
    retval = ex_get_conn(exodus_file_id, EX_ELEM_BLOCK, block_id, &block_elem_connectivity_.at(block_id)[0], 0, 0);
    if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_get_conn");
    // Switch from 1-based indexing to 0-based indexing
    for (unsigned int i = 0; i < block_elem_connectivity_.at(block_id).size(); i++) {
      block_elem_connectivity_.at(block_id).at(i) -= 1;
    }
  }

  retval = ex_close(exodus_file_id);
  if (retval != 0) ReportExodusError(retval, "GenesisMesh::ReadFile()", "ex_close");
#endif
}

void
GenesisMesh::ReadTextFile(std::string file_name)
{
  file_name_ = file_name;
  dim_       = 3;

  if (!IsValid()) { return; }

  std::stringstream error_msg_ss;

  std::string   text_file_name = file_name + ".txt";
  std::ifstream mesh_file(text_file_name.c_str());

  if (!mesh_file.is_open()) {
    std::stringstream error_msg_ss;
    error_msg_ss << "\n** Error, failed to open mesh file " << text_file_name << "\n";
    throw std::logic_error(error_msg_ss.str());
  }

  int num_nodes, num_elem, num_blocks, num_node_sets;

  int global_elem_index = 0;

  while (mesh_file.good()) {
    std::string str;
    getline(mesh_file, str);

    std::string       key;
    std::stringstream key_ss(str);
    key_ss >> key;

    std::stringstream ss(str);

    if (key == "number_of_nodes") {
      ss >> key;
      ss >> num_nodes;
    } else if (key == "number_of_elements") {
      ss >> key;
      ss >> num_elem;
      elem_global_id_.resize(num_elem);
    } else if (key == "number_of_blocks") {
      ss >> key;
      ss >> num_blocks;
    } else if (key == "number_of_node_sets") {
      ss >> key;
      ss >> num_node_sets;
    } else if (key == "node_coordinates") {
      node_global_id_.resize(num_nodes);
      node_x_.resize(num_nodes);
      node_y_.resize(num_nodes);
      node_z_.resize(num_nodes);
      for (int i = 0; i < num_nodes; i++) {
        getline(mesh_file, str);
        std::stringstream ss(str);
        int               global_node_id;
        double            x, y, z;
        ss >> global_node_id >> x >> y >> z;
        // switch to 0-based indexing
        node_global_id_[i] = global_node_id - 1;
        node_x_[i]         = x;
        node_y_[i]         = y;
        node_z_[i]         = z;
      }
    } else if (key == "element_block") {
      getline(mesh_file, str);
      std::stringstream ss(str);
      int               block_id, num_elem_in_block, num_nodes_per_elem;
      std::string       block_name;
      ss >> block_id >> block_name >> num_elem_in_block >> num_nodes_per_elem;
      block_ids_.push_back(block_id);
      block_names_[block_id]              = block_name;
      block_num_nodes_per_elem_[block_id] = num_nodes_per_elem;
      block_elem_connectivity_[block_id]  = std::vector<int>(num_elem_in_block * num_nodes_per_elem);
      block_elem_global_ids_[block_id]    = std::vector<int>(num_elem_in_block);
      std::vector<int>& conn              = block_elem_connectivity_[block_id];
      for (int i = 0; i < num_elem_in_block; i++) {
        getline(mesh_file, str);
        std::stringstream ss(str);
        int               elem_global_id;
        ss >> elem_global_id;
        // switch to 0-based indexing
        elem_global_id -= 1;
        elem_global_id_[global_elem_index++] = elem_global_id;
        block_elem_global_ids_[block_id][i]  = elem_global_id;
        for (int j = 0; j < num_nodes_per_elem; j++) {
          ss >> conn[i * num_nodes_per_elem + j];
          // switch to 0-based indexing
          conn[i * num_nodes_per_elem + j] -= 1;
        }
      }
    } else if (key == "nodeset") {
      getline(mesh_file, str);
      std::stringstream ss(str);
      int               node_set_id, node_set_number_of_nodes;
      std::string       node_set_name;
      ss >> node_set_id >> node_set_name >> node_set_number_of_nodes;
      node_set_ids_.push_back(node_set_id);
      node_set_names_[node_set_id] = node_set_name;
      node_sets_[node_set_id]      = std::vector<int>(node_set_number_of_nodes);
      std::vector<int>& node_set   = node_sets_[node_set_id];
      for (int i = 0; i < node_set_number_of_nodes; i++) {
        getline(mesh_file, str);
        std::stringstream ss(str);
        ss >> node_set[i];
        // switch to 0-based indexing
        node_set[i] -= 1;
      }
    }
  }
}

void
GenesisMesh::Initialize(
    std::string const&                     file_name,
    std::vector<int> const&                node_global_id,
    std::vector<double> const&             node_x,
    std::vector<double> const&             node_y,
    std::vector<double> const&             node_z,
    std::vector<int> const&                elem_global_id,
    std::vector<int> const&                block_ids,
    std::map<int, std::string> const&      block_names,
    std::map<int, std::vector<int>> const& block_elem_global_ids,
    std::map<int, int> const&              block_num_nodes_per_elem,
    std::map<int, std::vector<int>> const& block_elem_connectivity)
{
  file_name_                = file_name;
  dim_                      = 3;
  node_global_id_           = node_global_id;
  node_x_                   = node_x;
  node_y_                   = node_y;
  node_z_                   = node_z;
  elem_global_id_           = elem_global_id;
  block_ids_                = block_ids;
  all_block_ids_            = block_ids;
  block_names_              = block_names;
  all_block_names_          = block_names;
  block_elem_global_ids_    = block_elem_global_ids;
  block_num_nodes_per_elem_ = block_num_nodes_per_elem;
  block_elem_connectivity_  = block_elem_connectivity;
}

int
GenesisMesh::GetNumElementsInBlock(int block_id) const
{
  unsigned int conn_size          = block_elem_connectivity_.at(block_id).size();
  int          num_nodes_per_elem = block_num_nodes_per_elem_.at(block_id);
  int          num_elem_in_block  = 0;
  if (conn_size > 0) { num_elem_in_block = static_cast<int>(conn_size) / num_nodes_per_elem; }
  return num_elem_in_block;
}

std::map<int, int>
GenesisMesh::GetNumElementsInBlock() const
{
  std::map<int, int> num_elem_in_each_block;
  for (std::vector<int>::const_iterator it = block_ids_.begin(); it != block_ids_.end(); it++) {
    num_elem_in_each_block[*it] = GetNumElementsInBlock(*it);
  }
  return num_elem_in_each_block;
}

std::string
GenesisMesh::GetElementType(int block_id) const
{
  int         dim               = GetDim();
  int         num_node_per_elem = GetNumNodesPerElement(block_id);
  std::string elem_type         = "UNKNOWN";

  if (dim == 2) {
    switch (num_node_per_elem) {
      case 3: elem_type = "TRIANGLE:"; break;
      case 4: elem_type = "QUAD"; break;
      default: elem_type = "UNKNOWN"; break;
    }
  } else if (dim == 3) {
    switch (num_node_per_elem) {
      case 1: elem_type = "SPHERE"; break;
      case 3: elem_type = "TRIANGLE"; break;
      case 4: elem_type = "TETRA"; break;
      case 8: elem_type = "HEX"; break;
      default: elem_type = "UNKNOWN"; break;
    }
  }

  if (elem_type == "UNKNOWN") { throw std::logic_error("Error processing input mesh, unknown element type."); }

  return elem_type;
}

bool
GenesisMesh::HasBlock(std::string const& block_name) const
{
  for (auto& entry : block_names_) {
    if (entry.second == block_name) { return true; }
  }
  return false;
}

int
GenesisMesh::GetBlockId(std::string const& block_name) const
{
  for (auto& entry : block_names_) {
    if (entry.second == block_name) { return entry.first; }
  }
  throw std::logic_error("\n**** Error in GenesisMesh::GetBlockId(), block name not found.\n");
}

void
GenesisMesh::BlockNamesToOnProcessorBlockIds(std::vector<std::string> const& block_names, std::vector<int>& block_ids)
{
  block_ids.clear();
  for (auto& name : block_names) {
    // if the block name is not on this processor, do not add it to the list of
    // block ids
    if (HasBlock(name)) { block_ids.push_back(GetBlockId(name)); }
  }
}

void
GenesisMesh::BoundingBox(double& x_min, double& x_max, double& y_min, double& y_max, double& z_min, double& z_max) const
{
  // THIS WORKS ONLY IN SERIAL

  unsigned int num_nodes = GetNumNodes();
  double       big       = std::numeric_limits<double>::max();
  x_min                  = big;
  x_max                  = -1.0 * big;
  y_min                  = big;
  y_max                  = -1.0 * big;
  z_min                  = big;
  z_max                  = -1.0 * big;
  for (unsigned int n = 0; n < num_nodes; n++) {
    double x = node_x_[n];
    double y = node_y_[n];
    double z = node_z_[n];
    if (x < x_min) x_min = x;
    if (x > x_max) x_max = x;
    if (y < y_min) y_min = y;
    if (y > y_max) y_max = y;
    if (z < z_min) z_min = z;
    if (z > z_max) z_max = z;
  }
}

std::vector<double>
GenesisMesh::BoundingBoxCenter() const
{
  double x_min, x_max, y_min, y_max, z_min, z_max;
  BoundingBox(x_min, x_max, y_min, y_max, z_min, z_max);
  std::vector<double> center(3);
  center[0] = (x_max + x_min) / 2.0;
  center[1] = (y_max + y_min) / 2.0;
  center[2] = (z_max + z_min) / 2.0;
  return center;
}

void
GenesisMesh::AppendPeriodicPair(
    int                 local_primary_node_id,
    int                 local_secondary_node_id,
    const int* const    global_node_ids,
    std::map<int, int>& global_node_id_secondary_to_primary) const
{
  int global_primary_node_id                                    = global_node_ids[local_primary_node_id];
  int global_secondary_node_id                                  = global_node_ids[local_secondary_node_id];
  global_node_id_secondary_to_primary[global_secondary_node_id] = global_primary_node_id;
}

void
GenesisMesh::Print(bool verbose, int my_rank) const
{
  int num_elem = GetNumElements();

  std::cout << "\nExodus mesh on rank " << my_rank << ":" << std::endl;
  std::cout << "  file name:  " << file_name_ << std::endl;
  std::cout << "  dimension:  " << dim_ << std::endl;
  std::cout << "  num nodes:  " << node_x_.size() << std::endl;
  std::cout << "  num elems:  " << num_elem << std::endl;
  std::cout << "  num blocks: " << block_ids_.size() << std::endl;

  if (verbose) {
    std::cout << "  Block information:" << std::endl;
    for (int i = 0; i < block_ids_.size(); i++) {
      int block_id = block_ids_[i];
      std::cout << "    block id:   " << block_id << std::endl;
      std::cout << "    block name: " << block_names_.at(block_id) << std::endl;
      std::cout << "    elem type:  " << GetElementType(block_id) << std::endl;
      std::cout << "    num elems:  "
                << block_elem_connectivity_.at(block_id).size() / block_num_nodes_per_elem_.at(block_id) << std::endl;
    }
  }
  std::cout << std::endl;
  std::cout << std::flush;
}

void
GenesisMesh::ReportExodusError(int error_code, const char* method_name, const char* exodus_method_name) const
{
  std::stringstream ss;
  if (error_code < 0) {
    ss << "\n**Error GenesisMesh::" << method_name << ", Error code: " << error_code << " (" << exodus_method_name
       << ")\n";
    throw std::logic_error(ss.str());
  } else {
    ss << "\n**Warning GenesisMesh::" << method_name << ", Warning code: " << error_code << " (" << exodus_method_name
       << ")\n";
    std::cout << ss.str() << std::endl;
  }
}

}  // namespace nimble
