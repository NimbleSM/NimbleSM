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

#include "nimble_contact_manager.h"
#include "nimble_contact_interface.h"
#include "nimble_utils.h"
#include "mpi_buckets/src/CollisionManager.hpp"

#ifdef NIMBLE_HAVE_MPI
  #include "mpi.h"
#endif

#ifdef NIMBLE_HAVE_BVH
  #include <bvh/broadphase.hpp>
  #include <bvh/narrowphase.hpp>
  #include <bvh/tree.hpp>
  #include <bvh/vis/vis_bvh.hpp>
  #include <bvh/patch.hpp>
#ifdef NIMBLE_HAVE_MPI
  #include <bvh/vt/mpi_interop.hpp>
  #include <bvh/vt/broadphase.hpp>
  #include <bvh/vt/tree_build.hpp>
#endif
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>

namespace nimble {

  void ParseContactCommand(std::string const & command,
                           std::vector<std::string> & master_block_names,
                           std::vector<std::string> & slave_block_names,
                           double & penalty_parameter) {

    std::string contact_master_key;
    std::string contact_slave_key;

    std::stringstream ss(command);

    ss >> contact_master_key;
    if (contact_master_key != "master_blocks") {
      std::stringstream error_ss;
      error_ss << "\n**** Error processing contact command, unknown key: " << contact_master_key << std::endl;
      throw std::logic_error(error_ss.str());
    }

    bool slave_key_found = false;
    while (ss.good() && !slave_key_found) {
      std::string temp;
      ss >> temp;
      if (temp == "slave_blocks") {
        slave_key_found = true;
      }
      else {
        master_block_names.push_back(temp);
      }
    }

    if (!slave_key_found) {
      throw std::logic_error("\n**** Error processing contact command, expected \"slave_blocks\".\n");
    }

    bool penalty_parameter_key_found = false;
    while (ss.good() && !penalty_parameter_key_found) {
      std::string temp;
      ss >> temp;
      if (temp == "penalty_parameter") {
        penalty_parameter_key_found = true;
      }
      else {
        slave_block_names.push_back(temp);
      }
    }

    if (!penalty_parameter_key_found) {
      throw std::logic_error("\n**** Error processing contact command, expected \"penalty_parameter\".\n");
    }

    ss >> penalty_parameter;
  }

  void
  ContactManager::SkinBlocks(GenesisMesh const & mesh,
                             std::vector<int> const & block_ids,
                             int entity_id_offset,
                             std::vector< std::vector<int> > & skin_faces,
                             std::vector<int> & entity_ids) {

    std::map< std::vector<int>, std::vector<int> > faces;
    std::map< std::vector<int>, std::vector<int> >::iterator face_it;

    for (auto & block_id : block_ids) {

      int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      int num_node_per_elem = mesh.GetNumNodesPerElement(block_id);
      const int * const conn = mesh.GetConnectivity(block_id);
      std::vector<int> const & elem_global_ids = mesh.GetElementGlobalIdsInBlock(block_id);
      int conn_index = 0;

      // key is sorted node list for a face
      int num_node_per_face = 4;
      std::vector<int> key(num_node_per_face);
      // value is count, unsorted node list, exodus element id, and face ordinal
      std::vector<int> value(num_node_per_face + 3);

      for (int i_elem = 0 ; i_elem < num_elem_in_block ; i_elem++) {

        // switch from 0-based indexing to 1-based indexing
        // so that the ids will be valid exodus ids in the contact visualization output
        int elem_global_id = elem_global_ids[i_elem] + 1;

        // Examine each face, following the Exodus node-ordering convention

        // face 0: 0 1 5 4
        value[0] = 1;
        key[0] = value[1] = conn[conn_index + 0];
        key[1] = value[2] = conn[conn_index + 1];
        key[2] = value[3] = conn[conn_index + 5];
        key[3] = value[4] = conn[conn_index + 4];
        value[5] = elem_global_id;
        value[6] = 0;
        std::sort(key.begin(), key.end());
        face_it = faces.find(key);
        if (face_it == faces.end())
          faces[key] = value;
        else
          face_it->second[0] += 1;

        // face 1: 1 2 6 5
        value[0] = 1;
        key[0] = value[1] = conn[conn_index + 1];
        key[1] = value[2] = conn[conn_index + 2];
        key[2] = value[3] = conn[conn_index + 6];
        key[3] = value[4] = conn[conn_index + 5];
        value[5] = elem_global_id;
        value[6] = 1;
        std::sort(key.begin(), key.end());
        face_it = faces.find(key);
        if (face_it == faces.end())
          faces[key] = value;
        else
          face_it->second[0] += 1;

        // face 2: 2 3 7 6
        value[0] = 1;
        key[0] = value[1] = conn[conn_index + 2];
        key[1] = value[2] = conn[conn_index + 3];
        key[2] = value[3] = conn[conn_index + 7];
        key[3] = value[4] = conn[conn_index + 6];
        value[5] = elem_global_id;
        value[6] = 2;
        std::sort(key.begin(), key.end());
        face_it = faces.find(key);
        if (face_it == faces.end())
          faces[key] = value;
        else
          face_it->second[0] += 1;

        // face 3: 0 4 7 3
        value[0] = 1;
        key[0] = value[1] = conn[conn_index + 0];
        key[1] = value[2] = conn[conn_index + 4];
        key[2] = value[3] = conn[conn_index + 7];
        key[3] = value[4] = conn[conn_index + 3];
        value[5] = elem_global_id;
        value[6] = 3;
        std::sort(key.begin(), key.end());
        face_it = faces.find(key);
        if (face_it == faces.end())
          faces[key] = value;
        else
          face_it->second[0] += 1;

        // face 4: 0 3 2 1
        value[0] = 1;
        key[0] = value[1] = conn[conn_index + 0];
        key[1] = value[2] = conn[conn_index + 3];
        key[2] = value[3] = conn[conn_index + 2];
        key[3] = value[4] = conn[conn_index + 1];
        value[5] = elem_global_id;
        value[6] = 4;
        std::sort(key.begin(), key.end());
        face_it = faces.find(key);
        if (face_it == faces.end())
          faces[key] = value;
        else
          face_it->second[0] += 1;

        // face 5: 4 5 6 7
        value[0] = 1;
        key[0] = value[1] = conn[conn_index + 4];
        key[1] = value[2] = conn[conn_index + 5];
        key[2] = value[3] = conn[conn_index + 6];
        key[3] = value[4] = conn[conn_index + 7];
        value[5] = elem_global_id;
        value[6] = 5;
        std::sort(key.begin(), key.end());
        face_it = faces.find(key);
        if (face_it == faces.end())
          faces[key] = value;
        else
          face_it->second[0] += 1;

        conn_index += num_node_per_elem;
      }
    }

    skin_faces.clear();
    entity_ids.clear();
    for (auto face : faces) {
      if (face.second[0] == 1) {
        std::vector<int> skin_face;
        for (int i=0 ; i<face.second.size()-3 ; i++) {
          int id = face.second.at(i+1);
          skin_face.push_back(id);
        }
        skin_faces.push_back(skin_face);
        int entity_id = (face.second[5] + entity_id_offset) << 5;  // 59 bits for the genesis element id plus an offset value
        entity_id |= face.second[6] << 2;                          // 3 bits for the face ordinal
        entity_id |= 0;                                            // 2 bits for triangle ordinal (unknown until face is subdivided downstream)
        entity_ids.push_back(entity_id);
      }
      else if (face.second[0] != 2) {
        throw std::logic_error("Error in mesh skinning routine, face found more than two times!\n");
      }
    }

    return;
  }

  void
  ContactManager::RemoveInternalSkinFaces(GenesisMesh const & mesh,
                                          std::vector< std::vector<int> >& faces,
                                          std::vector<int>& entity_ids) {
#ifdef NIMBLE_HAVE_MPI

    using face_iterator_t = std::vector< std::vector<int> >::iterator;
    using face_id_iterator_t = std::vector<int>::iterator;

    int mpi_rank, num_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

    int num_nodes_in_face = 4;

    const int* genesis_node_global_ids = mesh.GetNodeGlobalIds();
    std::vector<int> face_global_ids;
    for (auto& face : faces) {
      for (int i=0 ; i<face.size() ; i++) {
        face_global_ids.push_back(genesis_node_global_ids[face[i]]);
      }
    }

    int max_buff_size = face_global_ids.size();
    int global_max_buff_size = 0;
    MPI_Allreduce(&max_buff_size, &global_max_buff_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    std::vector<int> mpi_buffer(global_max_buff_size);
    int mpi_buffer_num_faces = global_max_buff_size/num_nodes_in_face;

    for (int i_rank=0 ; i_rank<num_ranks ; i_rank++) {
      std::fill(mpi_buffer.begin(), mpi_buffer.end(), -1);
      if (mpi_rank == i_rank) {
        for (int i=0 ; i<face_global_ids.size() ; i++) {
          mpi_buffer[i] = face_global_ids.at(i);
        }
      }
      // DJL MPI COMMUNICATION SHOULD BE DONE ONLY WITH NEIGHBORHING PARTITIONS
      MPI_Bcast(mpi_buffer.data(), mpi_buffer.size(), MPI_INT, i_rank, MPI_COMM_WORLD);
      if (mpi_rank != i_rank) {
        for (int i_mpi_buff_face = 0 ; i_mpi_buff_face < mpi_buffer_num_faces ; i_mpi_buff_face++) {

          for (int i_face = 0 ; i_face < faces.size() ; i_face++) {
            std::vector<int> face_global_ids = faces[i_face];
            // DJL INEFFICIENT HANDLING OF GLOBAL IDS
            for (int i=0; i<face_global_ids.size() ; i++) {
              face_global_ids[i] = genesis_node_global_ids[face_global_ids[i]];
            }
            // DJL INCORRECTLY SEARCHING EMPTY ENTRIES AT THE END OF mpi_buffer
            if (std::find(face_global_ids.begin(), face_global_ids.end(), mpi_buffer.at(i_mpi_buff_face*num_nodes_in_face))   != face_global_ids.end() &&
                std::find(face_global_ids.begin(), face_global_ids.end(), mpi_buffer.at(i_mpi_buff_face*num_nodes_in_face+1)) != face_global_ids.end() &&
                std::find(face_global_ids.begin(), face_global_ids.end(), mpi_buffer.at(i_mpi_buff_face*num_nodes_in_face+2)) != face_global_ids.end() &&
                std::find(face_global_ids.begin(), face_global_ids.end(), mpi_buffer.at(i_mpi_buff_face*num_nodes_in_face+3)) != face_global_ids.end()) {
              face_iterator_t face_it = faces.begin() + i_face;
              face_id_iterator_t face_id_it = entity_ids.begin() + i_face;
              faces.erase(face_it);
              entity_ids.erase(face_id_it);
              break;
            }
          }
        }
      }
    }

#endif
  }

  ContactManager::ContactManager(std::shared_ptr<ContactInterface> interface,
                                 size_t dicing_factor)
    : penalty_parameter_(0.0),
      dicing_factor_(dicing_factor)
#if defined(NIMBLE_HAVE_MPI) && defined(NIMBLE_HAVE_BVH)
      , collision_world_(bvh::vt::context::current()->num_ranks() * dicing_factor, bvh::vt::context::current()->num_ranks() * dicing_factor),
      face_patch_collection_(bvh::vt::index_1d(bvh::vt::context::current()->num_ranks() * static_cast<int>(dicing_factor))),
      node_patch_collection_(bvh::vt::index_1d(bvh::vt::context::current()->num_ranks() * static_cast<int>(dicing_factor)))
#endif
  , contact_interface(interface)
  {
  }

  void
  ContactManager::CreateContactEntities(GenesisMesh const & mesh,
                                        nimble::MPIContainer & mpi_container,
                                        std::vector<int> const & master_block_ids,
                                        std::vector<int> const & slave_block_ids) {

    int mpi_rank = 0;
    int num_ranks = 1;
#ifdef NIMBLE_HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#endif

    contact_enabled_ = true;

    const double* coord_x = mesh.GetCoordinatesX();
    const double* coord_y = mesh.GetCoordinatesY();
    const double* coord_z = mesh.GetCoordinatesZ();

    int max_node_global_id = mesh.GetMaxNodeGlobalId();
#ifdef NIMBLE_HAVE_MPI
    int global_max_node_global_id = max_node_global_id;
    MPI_Allreduce(&max_node_global_id, &global_max_node_global_id, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    max_node_global_id = global_max_node_global_id;
#endif

    // find all the element faces on the master and slave contact blocks
    // the entity ids created here will be used downstream for the contact faces
    // the contact nodes will not use these entity ids, instead they will just use the exodus id of the node as the entity id
    std::vector< std::vector<int> > master_skin_faces, slave_skin_faces;
    std::vector<int> master_skin_entity_ids, slave_skin_entity_ids;
    int contact_entity_id_offset = max_node_global_id; // ensure that there are no duplicate entity ids between nodes and faces
    SkinBlocks(mesh, master_block_ids, contact_entity_id_offset, master_skin_faces, master_skin_entity_ids);
    SkinBlocks(mesh, slave_block_ids, contact_entity_id_offset, slave_skin_faces, slave_skin_entity_ids);

    // remove faces that are along partition boundaries
    RemoveInternalSkinFaces(mesh, master_skin_faces, master_skin_entity_ids);
    RemoveInternalSkinFaces(mesh, slave_skin_faces, slave_skin_entity_ids);

    // create a list of ghosted nodes (e.g., nodes that are owned by a different processor)
    std::vector<int> partition_boundary_node_local_ids;
    std::vector<int> min_rank_containing_partition_boundary_nodes;
    mpi_container.GetPartitionBoundaryNodeLocalIds(partition_boundary_node_local_ids,
                                                   min_rank_containing_partition_boundary_nodes);
    std::vector<int> ghosted_node_local_ids;
    for (int i=0 ; i<partition_boundary_node_local_ids.size() ; i++) {
      if (min_rank_containing_partition_boundary_nodes[i] != mpi_rank) {
        ghosted_node_local_ids.push_back(partition_boundary_node_local_ids[i]);
      }
    }

    // construct containers for the subset of the model that is involved with contact
    // this constitutes a submodel that is stored in the ContactManager
    std::set<int> node_ids_set;
    for (auto const & face : master_skin_faces) {
      for (auto const & id : face) {
        node_ids_set.insert(id);
      }
    }
    for (auto const & face : slave_skin_faces) {
      for (auto const & id : face) {
        node_ids_set.insert(id);
      }
    }
    node_ids_ = std::vector<int>(node_ids_set.begin(), node_ids_set.end());

    std::map<int, int> genesis_mesh_node_id_to_contact_submodel_id;
    for (unsigned int i_node=0 ; i_node<node_ids_.size() ; ++i_node) {
      genesis_mesh_node_id_to_contact_submodel_id[node_ids_[i_node]] = i_node;
    }

    // replace the node ids that correspond to the genesis mesh
    // with node ids that correspond to the contact submodel
    for (unsigned int i_face=0 ; i_face<master_skin_faces.size() ; ++i_face) {
      for (unsigned int i_node=0 ; i_node<master_skin_faces[i_face].size() ; ++i_node) {
        int genesis_mesh_node_id = master_skin_faces[i_face][i_node];
        master_skin_faces[i_face][i_node] = genesis_mesh_node_id_to_contact_submodel_id.at(genesis_mesh_node_id);
      }
    }
    for (unsigned int i_face=0 ; i_face<slave_skin_faces.size() ; ++i_face) {
      for (unsigned int i_node=0 ; i_node<slave_skin_faces[i_face].size() ; ++i_node) {
        int genesis_mesh_node_id = slave_skin_faces[i_face][i_node];
        slave_skin_faces[i_face][i_node] = genesis_mesh_node_id_to_contact_submodel_id.at(genesis_mesh_node_id);
      }
    }

    // create a list of ghosed nodes in the contact submodel
    std::vector<int> ghosted_contact_node_ids;
    for (auto node_id : ghosted_node_local_ids) {
      std::map<int, int>::iterator it = genesis_mesh_node_id_to_contact_submodel_id.find(node_id);
      if (it != genesis_mesh_node_id_to_contact_submodel_id.end()) {
        ghosted_contact_node_ids.push_back(it->second);
      }
    }

    // allocate data for the contact submodel
    int array_len = 3*node_ids_.size();
    model_coord_.resize(array_len);
    coord_.resize(array_len);
    force_.resize(array_len, 0.0);
    for (unsigned int i_node=0 ; i_node<node_ids_.size() ; i_node++) {
      model_coord_[3*i_node]   = coord_[3*i_node]   = coord_x[node_ids_[i_node]];
      model_coord_[3*i_node+1] = coord_[3*i_node+1] = coord_y[node_ids_[i_node]];
      model_coord_[3*i_node+2] = coord_[3*i_node+2] = coord_z[node_ids_[i_node]];
    }

    // Store nodes in slave faces
    // Create a list of nodes and their characteristic lengths
    const int* genesis_node_global_ids = mesh.GetNodeGlobalIds();
    std::vector<int> slave_node_ids;
    std::vector<int> slave_node_entity_ids;
    std::map<int, double> slave_node_char_lens;
    for (int i_face=0 ; i_face<slave_skin_faces.size(); i_face++) {

      std::vector<int>& face = slave_skin_faces[i_face];

      int num_nodes_in_face = static_cast<int>(face.size());
      // determine a characteristic length based on max edge length
      double max_edge_length = std::numeric_limits<double>::lowest();
      for (int i=0 ; i<num_nodes_in_face ; ++i) {
        int node_id_1 = face[i];
        int node_id_2 = face[0];
        if (i+1 < num_nodes_in_face) {
          node_id_2 = face[i+1];
        }
        double edge_length = sqrt( (coord_[3*node_id_2  ] - coord_[3*node_id_1  ])*(coord_[3*node_id_2  ] - coord_[3*node_id_1  ]) +
                                   (coord_[3*node_id_2+1] - coord_[3*node_id_1+1])*(coord_[3*node_id_2+1] - coord_[3*node_id_1+1]) +
                                   (coord_[3*node_id_2+2] - coord_[3*node_id_1+2])*(coord_[3*node_id_2+2] - coord_[3*node_id_1+2]) );
        if (edge_length > max_edge_length) {
          max_edge_length = edge_length;
        }
      }
      double characteristic_length = max_edge_length;
      for (int i_node=0 ; i_node<num_nodes_in_face; i_node++) {
        int node_id = face[i_node];
        // omit ghosted nodes
        if (std::find(ghosted_contact_node_ids.begin(), ghosted_contact_node_ids.end(), node_id) == ghosted_contact_node_ids.end()) {
          if (std::find(slave_node_ids.begin(), slave_node_ids.end(), node_id) == slave_node_ids.end()) {
            slave_node_ids.push_back(node_id);
            // note the mapping from contact manager local id to real FEM mesh local id to real FEM mesh global id
            int contact_node_entity_id = genesis_node_global_ids[node_ids_[node_id]] + 1;
            slave_node_entity_ids.push_back(contact_node_entity_id);
            slave_node_char_lens[node_id] = characteristic_length;
          }
          else {
            // always use the maximum characteristic length
            // this requires a parallel sync
            if (slave_node_char_lens[node_id] < characteristic_length) {
              slave_node_char_lens[node_id] = characteristic_length;
            }
          }
        }
      }
    }

    contact_nodes_.resize(slave_node_ids.size());
    contact_faces_.resize(4*master_skin_faces.size());
    CreateContactNodesAndFaces(master_skin_faces, master_skin_entity_ids, slave_node_ids, slave_node_entity_ids, slave_node_char_lens, contact_nodes_, contact_faces_);

#ifdef NIMBLE_HAVE_KOKKOS

    nimble_kokkos::HostIntegerArrayView node_ids_h("contact_node_ids_h", node_ids_.size());
    for (unsigned int i_node=0 ; i_node<node_ids_.size() ; i_node++) {
      node_ids_h[i_node] = node_ids_[i_node];
    }

    nimble_kokkos::HostScalarNodeView model_coord_h("contact_model_coord_h", array_len);
    for (unsigned int i_node=0 ; i_node<node_ids_.size() ; i_node++) {
      model_coord_h[3*i_node]   = coord_x[node_ids_[i_node]];
      model_coord_h[3*i_node+1] = coord_y[node_ids_[i_node]];
      model_coord_h[3*i_node+2] = coord_z[node_ids_[i_node]];
    }

    Kokkos::resize(node_ids_d_, node_ids_.size());
    Kokkos::resize(model_coord_d_, array_len);
    Kokkos::resize(coord_d_, array_len);
    Kokkos::resize(force_d_, array_len);

    Kokkos::deep_copy(node_ids_d_, node_ids_h);
    Kokkos::deep_copy(model_coord_d_, model_coord_h);
    Kokkos::deep_copy(coord_d_, model_coord_h);
    Kokkos::deep_copy(force_d_, 0.0);

    Kokkos::resize(contact_nodes_h_, slave_node_ids.size());
    Kokkos::resize(contact_faces_h_, 4*master_skin_faces.size());
    CreateContactNodesAndFaces(master_skin_faces, master_skin_entity_ids, slave_node_ids, slave_node_entity_ids, slave_node_char_lens, contact_nodes_h_, contact_faces_h_);

    Kokkos::resize(contact_nodes_d_, slave_node_ids.size());
    Kokkos::resize(contact_faces_d_, 4*master_skin_faces.size());
    Kokkos::deep_copy(contact_nodes_d_, contact_nodes_h_);
    Kokkos::deep_copy(contact_faces_d_, contact_faces_h_);
#endif

    int num_contact_faces = contact_faces_.size();
    int num_contact_nodes = contact_nodes_.size();
#ifdef NIMBLE_HAVE_MPI
    std::vector<int> input(2);
    std::vector<int> output(2);
    input[0] = num_contact_faces;
    input[1] = num_contact_nodes;
    MPI_Reduce(input.data(), output.data(), input.size(), MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    num_contact_faces = output[0];
    num_contact_nodes = output[1];
#endif
    if (mpi_rank == 0) {
      std::cout << "Contact initialization:" << std::endl;
      std::cout << "  number of triangular contact facets (master blocks): " << num_contact_faces << std::endl;
      std::cout << "  number of contact nodes (slave blocks): " << num_contact_nodes << "\n" << std::endl;
    }
  }

  template <typename ArgT>
  void ContactManager::CreateContactNodesAndFaces(std::vector< std::vector<int> > const & master_skin_faces,
                                                  std::vector<int> const & master_skin_entity_ids,
                                                  std::vector<int> const & slave_node_ids,
                                                  std::vector<int> const & slave_node_entity_ids,
                                                  std::map<int, double> const & slave_node_char_lens,
                                                  ArgT& contact_nodes,
                                                  ArgT& contact_faces) const {

    int index = 0;

    // convert master faces to trangular facets
    for (unsigned int i_face=0 ; i_face < master_skin_faces.size() ; i_face++) {

      auto face = master_skin_faces[i_face];

      int num_nodes_in_face = static_cast<int>(face.size());
      if (num_nodes_in_face != 4) {
        throw std::logic_error("\nError in ContactManager::CreateContactNodesAndFaces(), invalid number of face nodes.\n");
      }

      // determine a characteristic length based on max edge length
      double max_edge_length = std::numeric_limits<double>::lowest();
      for (int i=0 ; i<num_nodes_in_face ; ++i) {
        int node_id_1 = face[i];
        int node_id_2 = face[0];
        if (i+1 < num_nodes_in_face) {
          node_id_2 = face[i+1];
        }
        double edge_length = sqrt( (coord_[3*node_id_2  ] - coord_[3*node_id_1  ])*(coord_[3*node_id_2  ] - coord_[3*node_id_1  ]) +
                                   (coord_[3*node_id_2+1] - coord_[3*node_id_1+1])*(coord_[3*node_id_2+1] - coord_[3*node_id_1+1]) +
                                   (coord_[3*node_id_2+2] - coord_[3*node_id_1+2])*(coord_[3*node_id_2+2] - coord_[3*node_id_1+2]) );
        if (edge_length > max_edge_length) {
          max_edge_length = edge_length;
        }
      }
      double characteristic_length = max_edge_length;

      // create a node at the barycenter of the face
      double fictitious_node[3] = {0.0, 0.0, 0.0};
      for (int i=0 ; i<num_nodes_in_face ; ++i) {
        int node_id = face[i];
        for (int j=0 ; j<3 ; j++) {
          fictitious_node[j] += coord_[3*node_id+j];
        }
      }
      for (int j=0 ; j<3 ; j++) {
        fictitious_node[j] /= num_nodes_in_face;
      }

      // Create a map for transfering displacements and contact forces from the nodes on the
      // triangle patch to the contact manager data structures.  There is a 1-to-1 transfer for the two real nodes,
      // and for the fictitious node the mapping applies an equal fraction of the displacement/force
      // at the fictitious node to each for four real nodes in the original mesh face
      int node_ids_for_fictitious_node[4];
      for (int i=0 ; i<4 ; i++){
        node_ids_for_fictitious_node[i] = face[i];
      }

      double model_coord[9];
      int node_id_1, node_id_2, entity_id;

      // triangle node_0, node_1, fictitious_node
      node_id_1 = face[0];
      node_id_2 = face[1];
      for (int i=0 ; i<3 ; ++i) {
        model_coord[i] = coord_[3*node_id_1+i];
        model_coord[3+i] = coord_[3*node_id_2+i];
      }
      model_coord[6] = fictitious_node[0];
      model_coord[7] = fictitious_node[1];
      model_coord[8] = fictitious_node[2];
      entity_id = master_skin_entity_ids[i_face];
      entity_id |= 0; // triangle ordinal
      contact_faces[index++] = ContactEntity(ContactEntity::TRIANGLE,
                                             entity_id,
                                             model_coord,
                                             characteristic_length,
                                             node_id_1,
                                             node_id_2,
                                             node_ids_for_fictitious_node);

      // triangle node_1, node_2, fictitious_node
      node_id_1 = face[1];
      node_id_2 = face[2];
      for (int i=0 ; i<3 ; ++i) {
        model_coord[i] = coord_[3*node_id_1+i];
        model_coord[3+i] = coord_[3*node_id_2+i];
      }
      model_coord[6] = fictitious_node[0];
      model_coord[7] = fictitious_node[1];
      model_coord[8] = fictitious_node[2];
      entity_id = master_skin_entity_ids[i_face];
      entity_id |= 1; // triangle ordinal
      contact_faces[index++] = ContactEntity(ContactEntity::TRIANGLE,
                                             entity_id,
                                             model_coord,
                                             characteristic_length,
                                             node_id_1,
                                             node_id_2,
                                             node_ids_for_fictitious_node);

      // triangle node_2, node_3, fictitious_node
      node_id_1 = face[2];
      node_id_2 = face[3];
      for (int i=0 ; i<3 ; ++i) {
        model_coord[i] = coord_[3*node_id_1+i];
        model_coord[3+i] = coord_[3*node_id_2+i];
      }
      model_coord[6] = fictitious_node[0];
      model_coord[7] = fictitious_node[1];
      model_coord[8] = fictitious_node[2];
      entity_id = master_skin_entity_ids[i_face];
      entity_id |= 2; // triangle ordinal
      contact_faces[index++] = ContactEntity(ContactEntity::TRIANGLE,
                                             entity_id,
                                             model_coord,
                                             characteristic_length,
                                             node_id_1,
                                             node_id_2,
                                             node_ids_for_fictitious_node);

      // triangle node_3, node_0, fictitious_node
      node_id_1 = face[3];
      node_id_2 = face[0];
      for (int i=0 ; i<3 ; ++i) {
        model_coord[i] = coord_[3*node_id_1+i];
        model_coord[3+i] = coord_[3*node_id_2+i];
      }
      model_coord[6] = fictitious_node[0];
      model_coord[7] = fictitious_node[1];
      model_coord[8] = fictitious_node[2];
      entity_id = master_skin_entity_ids[i_face];
      entity_id |= 3; // triangle ordinal
      contact_faces[index++] = ContactEntity(ContactEntity::TRIANGLE,
                                             entity_id,
                                             model_coord,
                                             characteristic_length,
                                             node_id_1,
                                             node_id_2,
                                             node_ids_for_fictitious_node);
    }

    // Slave node entities
    for (unsigned int i_node=0 ; i_node < slave_node_ids.size() ; ++i_node) {
      int node_id = slave_node_ids.at(i_node);
      int entity_id = slave_node_entity_ids.at(i_node);
      double characteristic_length = slave_node_char_lens.at(node_id);
      double model_coord[3];
      for (int i=0 ; i<3 ; ++i) {
        model_coord[i] = coord_[3*node_id+i];
      }
      contact_nodes[i_node] = ContactEntity(ContactEntity::NODE,
                                            entity_id,
                                            model_coord,
                                            characteristic_length,
                                            node_id);
    }
  }

  void
  ContactManager::BoundingBox(double& x_min,
                              double& x_max,
                              double& y_min,
                              double& y_max,
                              double& z_min,
                              double& z_max) const {

    double big = std::numeric_limits<double>::max();

#ifdef NIMBLE_HAVE_KOKKOS
    nimble_kokkos::DeviceScalarNodeView contact_bounding_box_d("contact_bounding_box_d", 6);
    nimble_kokkos::HostScalarNodeView contact_bounding_box_h("contact_bounding_box_h", 6);
    contact_bounding_box_h(0) = big;       // x_min
    contact_bounding_box_h(1) = -1.0*big;  // x_max
    contact_bounding_box_h(2) = big;       // y_min
    contact_bounding_box_h(3) = -1.0*big;  // y_max
    contact_bounding_box_h(4) = big;       // z_min
    contact_bounding_box_h(5) = -1.0*big;  // z_max
    Kokkos::deep_copy(contact_bounding_box_d, contact_bounding_box_h);

    nimble_kokkos::DeviceScalarNodeView coord_d = coord_d_;
    int contact_vector_size = coord_d.extent(0) / 3;

    Kokkos::parallel_for("Contact Bounding Box",
                         contact_vector_size,
                         KOKKOS_LAMBDA(const int i) {
      double x = coord_d(3*i);
      double y = coord_d(3*i+1);
      double z = coord_d(3*i+2);
      Kokkos::atomic_min_fetch(&contact_bounding_box_d(0), x);
      Kokkos::atomic_max_fetch(&contact_bounding_box_d(1), x);
      Kokkos::atomic_min_fetch(&contact_bounding_box_d(2), y);
      Kokkos::atomic_max_fetch(&contact_bounding_box_d(3), y);
      Kokkos::atomic_min_fetch(&contact_bounding_box_d(4), z);
      Kokkos::atomic_max_fetch(&contact_bounding_box_d(5), z);
    });

    Kokkos::deep_copy(contact_bounding_box_h, contact_bounding_box_d);
    x_min = contact_bounding_box_h(0);
    x_max = contact_bounding_box_h(1);
    y_min = contact_bounding_box_h(2);
    y_max = contact_bounding_box_h(3);
    z_min = contact_bounding_box_h(4);
    z_max = contact_bounding_box_h(5);
#else
    x_min = big;
    x_max = -1.0 * big;
    y_min = big;
    y_max = -1.0 * big;
    z_min = big;
    z_max = -1.0 * big;
    for (unsigned int i=0 ; i<coord_.size()/3 ; i++) {
      double x = coord_[i*3];
      double y = coord_[i*3+1];
      double z = coord_[i*3+2];
      if (x < x_min)
        x_min = x;
      if (x > x_max)
        x_max = x;
      if (y < y_min)
        y_min = y;
      if (y > y_max)
        y_max = y;
      if (z < z_min)
        z_min = z;
      if (z > z_max)
        z_max = z;
    }
#endif
  }

  double
  ContactManager::BoundingBoxAverageCharacteristicLengthOverAllRanks() const {
    double x_min, x_max, y_min, y_max, z_min, z_max;
    BoundingBox(x_min, x_max, y_min, y_max, z_min, z_max);
    double longest_edge = x_max - x_min;
    if ((y_max - y_min) > longest_edge) {
      longest_edge = y_max - y_min;
    }
    if ((z_max - z_min) > longest_edge) {
      longest_edge = z_max - z_min;
    }
    double ave_characteristic_length = longest_edge;
#ifdef NIMBLE_HAVE_MPI
    int num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    MPI_Allreduce(&longest_edge, &ave_characteristic_length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ave_characteristic_length /= num_ranks;
#endif
    return ave_characteristic_length;
  }

  void
  ContactManager::InitializeContactVisualization(std::string const & contact_visualization_exodus_file_name){

    // Exodus id convention for contact visualization:
    //
    // Both node and face contact entities have a unique, parallel-consistent id called contact_entity_global_id_.
    // For faces, the contact_entity_global_id_ a bit-wise combination of the global exodus id of the parent element,
    // plus the face ordinal (1-6), plus the triangle ordinal (1-4).
    // For nodes, the contact_entity_global_id_ is the exodus global node id for the node in the original FEM mesh.
    //
    // For visualization output, we need unique, parallel-consistent node ids and element ids.  For the faces, the
    // contact_entity_global_id_ is used as the element id, and the node ids are constructed here.  For nodes, the
    // contact_entity_global_id_ is used for both the node id and the element id (sphere element containing a single node).
    //
    // For the MPI bounding boxes, both the node ids and the element id are constructed here.
    //
    //   contact faces:
    //     node ids are (3 * contact_entity_global_id_ + max_contact_entity_id + 9,
    //                   3 * contact_entity_global_id_ + max_contact_entity_id + 10,
    //                   3 * contact_entity_global_id_ + max_contact_entity_id + 11)
    //     element id is contact_entity_global_id_
    //   contact nodes
    //     node id is contact_entity_global_id_
    //     element id contact_entity_global_id_
    //   mpi partition bounding box:
    //     nodes id are (3 * max_contact_entity_id + 1,
    //                   3 * max_contact_entity_id + 2,
    //                   3 * max_contact_entity_id + 3,
    //                   3 * max_contact_entity_id + 4,
    //                   3 * max_contact_entity_id + 5,
    //                   3 * max_contact_entity_id + 6,
    //                   3 * max_contact_entity_id + 7,
    //                   3 * max_contact_entity_id + 8)
    //     element id max_contact_entity_id + 1

#ifdef NIMBLE_HAVE_BVH

#elif !defined(NIMBLE_HAVE_KOKKOS)
    throw std::logic_error("\nError in ContactManager::InitializeContactVisualization(), contact visualization currently available only for NimbleSM_Kokkos,\n");
#else

    // determine the maximum contact entity global id over all MPI partitions
    int max_contact_entity_id = 0;
    for (int i_face=0 ; i_face<contact_faces_h_.extent(0) ; i_face++) {
      if (contact_faces_h_[i_face].contact_entity_global_id_ > max_contact_entity_id) {
        max_contact_entity_id = contact_faces_h_[i_face].contact_entity_global_id_;
      }
    }
    for (int i_node=0 ; i_node<contact_nodes_h_.extent(0) ; i_node++) {
      if (contact_nodes_h_[i_node].contact_entity_global_id_ > max_contact_entity_id) {
        max_contact_entity_id = contact_nodes_h_[i_node].contact_entity_global_id_;
      }
    }
#ifdef NIMBLE_HAVE_MPI
    int mpi_rank, num_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    int global_max_contact_entity_id = max_contact_entity_id;
    MPI_Allreduce(&max_contact_entity_id, &global_max_contact_entity_id, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    max_contact_entity_id = global_max_contact_entity_id;
#endif

    std::vector<int> node_global_id;
    std::vector<double> node_x;
    std::vector<double> node_y;
    std::vector<double> node_z;
    std::vector<int> elem_global_id;
    std::vector<int> block_ids;
    std::map<int, std::string> block_names;
    std::map<int, std::vector<int> > block_elem_global_ids;
    std::map<int, int> block_num_nodes_per_elem;
    std::map<int, std::vector<int> > block_elem_connectivity;

    // first block contains the contact faces
    int block_id = 1;
    block_ids.push_back(block_id);
    block_names[block_id] = "contact_faces";
    block_elem_global_ids[block_id] = std::vector<int>();
    block_num_nodes_per_elem[block_id] = 3;
    block_elem_connectivity[block_id] = std::vector<int>();

    int node_index(0);
    for (int i_face=0 ; i_face<contact_faces_h_.extent(0) ; i_face++) {
      ContactEntity const & face = contact_faces_h_(i_face);
      int contact_entity_global_id = face.contact_entity_global_id_;
      node_global_id.push_back(3 * contact_entity_global_id + max_contact_entity_id + 9);
      node_x.push_back(face.coord_1_x_);
      node_y.push_back(face.coord_1_y_);
      node_z.push_back(face.coord_1_z_);
      block_elem_connectivity[block_id].push_back(node_index++);
      node_global_id.push_back(3 * contact_entity_global_id + max_contact_entity_id + 10);
      node_x.push_back(face.coord_2_x_);
      node_y.push_back(face.coord_2_y_);
      node_z.push_back(face.coord_2_z_);
      block_elem_connectivity[block_id].push_back(node_index++);
      node_global_id.push_back(3 * contact_entity_global_id + max_contact_entity_id + 11);
      node_x.push_back(face.coord_3_x_);
      node_y.push_back(face.coord_3_y_);
      node_z.push_back(face.coord_3_z_);
      block_elem_connectivity[block_id].push_back(node_index++);
      elem_global_id.push_back(contact_entity_global_id);
    }

    // second block contains the contact nodes
    block_id = 2;
    block_ids.push_back(block_id);
    block_names[block_id] = "contact_nodes";
    block_elem_global_ids[block_id] = std::vector<int>();
    block_num_nodes_per_elem[block_id] = 1;
    block_elem_connectivity[block_id] = std::vector<int>();

    for (int i_node=0 ; i_node<contact_nodes_h_.extent(0) ; i_node++) {
      ContactEntity const & node = contact_nodes_h_(i_node);
      int contact_entity_global_id = node.contact_entity_global_id_;
      node_global_id.push_back(contact_entity_global_id);
      node_x.push_back(node.coord_1_x_);
      node_y.push_back(node.coord_1_y_);
      node_z.push_back(node.coord_1_z_);
      block_elem_connectivity[block_id].push_back(node_index++);
      elem_global_id.push_back(contact_entity_global_id);
    }

    // third block is the bounding box for this mpi rank
    /**
    block_id = 3;
    block_ids.push_back(block_id);
    block_names[block_id] = "contact_mpi_rank_bounding_box";
    block_elem_global_ids[block_id] = std::vector<int>();
    block_num_nodes_per_elem[block_id] = 8;
    block_elem_connectivity[block_id] = std::vector<int>();

    double x_min, x_max, y_min, y_max, z_min, z_max;
    BoundingBox(x_min, x_max, y_min, y_max, z_min, z_max);
    node_global_id.push_back(3 * max_contact_entity_id + 1);
    node_x.push_back(x_min); node_y.push_back(y_min); node_z.push_back(z_max);
    block_elem_connectivity[block_id].push_back(node_index++);
    node_global_id.push_back(3 * max_contact_entity_id + 2);
    node_x.push_back(x_max); node_y.push_back(y_min); node_z.push_back(z_max);
    block_elem_connectivity[block_id].push_back(node_index++);
    node_global_id.push_back(3 * max_contact_entity_id + 3);
    node_x.push_back(x_max); node_y.push_back(y_min); node_z.push_back(z_min);
    block_elem_connectivity[block_id].push_back(node_index++);
    node_global_id.push_back(3 * max_contact_entity_id + 4);
    node_x.push_back(x_min); node_y.push_back(y_min); node_z.push_back(z_min);
    block_elem_connectivity[block_id].push_back(node_index++);
    node_global_id.push_back(3 * max_contact_entity_id + 5);
    node_x.push_back(x_min); node_y.push_back(y_max); node_z.push_back(z_max);
    block_elem_connectivity[block_id].push_back(node_index++);
    node_global_id.push_back(3 * max_contact_entity_id + 6);
    node_x.push_back(x_max); node_y.push_back(y_max); node_z.push_back(z_max);
    block_elem_connectivity[block_id].push_back(node_index++);
    node_global_id.push_back(3 * max_contact_entity_id + 7);
    node_x.push_back(x_max); node_y.push_back(y_max); node_z.push_back(z_min);
    block_elem_connectivity[block_id].push_back(node_index++);
    node_global_id.push_back(3 * max_contact_entity_id + 8);
    node_x.push_back(x_min); node_y.push_back(y_max); node_z.push_back(z_min);
    block_elem_connectivity[block_id].push_back(node_index++);
    elem_global_id.push_back(max_contact_entity_id + 1);

    // store the model coordinate bounding box
    contact_visualization_model_coord_bounding_box_[0] = x_min;
    contact_visualization_model_coord_bounding_box_[1] = x_max;
    contact_visualization_model_coord_bounding_box_[2] = y_min;
    contact_visualization_model_coord_bounding_box_[3] = y_max;
    contact_visualization_model_coord_bounding_box_[4] = z_min;
    contact_visualization_model_coord_bounding_box_[5] = z_max;
    **/
    genesis_mesh_for_contact_visualization_.Initialize("contact_visualization",
                                                       node_global_id,
                                                       node_x,
                                                       node_y,
                                                       node_z,
                                                       elem_global_id,
                                                       block_ids,
                                                       block_names,
                                                       block_elem_global_ids,
                                                       block_num_nodes_per_elem,
                                                       block_elem_connectivity);

    exodus_output_for_contact_visualization_.Initialize(contact_visualization_exodus_file_name,
                                                        genesis_mesh_for_contact_visualization_);

    std::vector<std::string> global_data_labels;
    std::vector<std::string> node_data_labels_for_output;
    node_data_labels_for_output.push_back("displacement_x");
    node_data_labels_for_output.push_back("displacement_y");
    node_data_labels_for_output.push_back("displacement_z");
    std::map<int, std::vector<std::string> > elem_data_labels_for_output;
    std::map<int, std::vector<std::string> > derived_elem_data_labels;
    for (auto & block_id : block_ids) {
      elem_data_labels_for_output[block_id] = std::vector<std::string>();
      derived_elem_data_labels[block_id] = std::vector<std::string>();
    }
    exodus_output_for_contact_visualization_.InitializeDatabase(genesis_mesh_for_contact_visualization_,
                                                                global_data_labels,
                                                                node_data_labels_for_output,
                                                                elem_data_labels_for_output,
                                                                derived_elem_data_labels);
#endif
  }

  void
  ContactManager::ContactVisualizationWriteStep(double time_current){
#ifdef NIMBLE_HAVE_BVH

#elif !defined(NIMBLE_HAVE_KOKKOS)
    throw std::logic_error("\nError in ContactManager::ContactVisualizationWriteStep(), contact visualization currently available only for NimbleSM_Kokkos,\n");
#else
    // copy contact entities from host to device
    Kokkos::deep_copy(contact_nodes_h_, contact_nodes_d_);
    Kokkos::deep_copy(contact_faces_h_, contact_faces_d_);

    std::vector<double> global_data;
    std::vector< std::vector<double> > node_data_for_output(3);
    std::map<int, std::vector<std::string> > elem_data_labels_for_output;
    std::map<int, std::vector< std::vector<double> > > elem_data_for_output;
    std::map<int, std::vector<std::string> > derived_elem_data_labels;
    std::map<int, std::vector< std::vector<double> > > derived_elem_data;

    std::vector<int> const & block_ids = genesis_mesh_for_contact_visualization_.GetBlockIds();
    for (auto & block_id : block_ids) {
      elem_data_labels_for_output[block_id] = std::vector<std::string>();
      derived_elem_data_labels[block_id] = std::vector<std::string>();
    }

    // node_data_for_output contains displacement_x, displacement_y, displacement_z
    int num_nodes = genesis_mesh_for_contact_visualization_.GetNumNodes();
    node_data_for_output[0].resize(num_nodes);
    node_data_for_output[1].resize(num_nodes);
    node_data_for_output[2].resize(num_nodes);
    const double * model_coord_x = genesis_mesh_for_contact_visualization_.GetCoordinatesX();
    const double * model_coord_y = genesis_mesh_for_contact_visualization_.GetCoordinatesY();
    const double * model_coord_z = genesis_mesh_for_contact_visualization_.GetCoordinatesZ();

    int node_index(0);
    for (int i_face=0 ; i_face<contact_faces_h_.extent(0) ; i_face++) {
      ContactEntity const & face = contact_faces_h_[i_face];
      node_data_for_output[0][node_index] = face.coord_1_x_ - model_coord_x[node_index];
      node_data_for_output[1][node_index] = face.coord_1_y_ - model_coord_y[node_index];
      node_data_for_output[2][node_index] = face.coord_1_z_ - model_coord_z[node_index];
      node_index += 1;
      node_data_for_output[0][node_index] = face.coord_2_x_ - model_coord_x[node_index];
      node_data_for_output[1][node_index] = face.coord_2_y_ - model_coord_y[node_index];
      node_data_for_output[2][node_index] = face.coord_2_z_ - model_coord_z[node_index];
      node_index += 1;
      node_data_for_output[0][node_index] = face.coord_3_x_ - model_coord_x[node_index];
      node_data_for_output[1][node_index] = face.coord_3_y_ - model_coord_y[node_index];
      node_data_for_output[2][node_index] = face.coord_3_z_ - model_coord_z[node_index];
      node_index += 1;
    }
    for (int i_node=0 ; i_node<contact_nodes_h_.extent(0) ; i_node++) {
      ContactEntity const & node = contact_nodes_h_(i_node);
      node_data_for_output[0][node_index] = node.coord_1_x_ - model_coord_x[node_index];
      node_data_for_output[1][node_index] = node.coord_1_y_ - model_coord_y[node_index];
      node_data_for_output[2][node_index] = node.coord_1_z_ - model_coord_z[node_index];
      node_index += 1;
    }
    /*
    double x_min_model_coord = contact_visualization_model_coord_bounding_box_[0];
    double x_max_model_coord = contact_visualization_model_coord_bounding_box_[1];
    double y_min_model_coord = contact_visualization_model_coord_bounding_box_[2];
    double y_max_model_coord = contact_visualization_model_coord_bounding_box_[3];
    double z_min_model_coord = contact_visualization_model_coord_bounding_box_[4];
    double z_max_model_coord = contact_visualization_model_coord_bounding_box_[5];
    double x_min, x_max, y_min, y_max, z_min, z_max;
    BoundingBox(x_min, x_max, y_min, y_max, z_min, z_max);
    // node_x.push_back(x_min); node_y.push_back(y_min); node_z.push_back(z_max);
    node_data_for_output[0][node_index] = x_min - x_min_model_coord;
    node_data_for_output[1][node_index] = y_min - y_min_model_coord;
    node_data_for_output[2][node_index] = z_min - z_min_model_coord;
    node_index += 1;
    // node_x.push_back(x_max); node_y.push_back(y_min); node_z.push_back(z_max);
    node_data_for_output[0][node_index] = x_max - x_max_model_coord;
    node_data_for_output[1][node_index] = y_min - y_min_model_coord;
    node_data_for_output[2][node_index] = z_max - z_max_model_coord;
    node_index += 1;
    // node_x.push_back(x_max); node_y.push_back(y_min); node_z.push_back(z_min);
    node_data_for_output[0][node_index] = x_max - x_max_model_coord;
    node_data_for_output[1][node_index] = y_min - y_min_model_coord;
    node_data_for_output[2][node_index] = z_min - z_min_model_coord;
    node_index += 1;
    // node_x.push_back(x_min); node_y.push_back(y_min); node_z.push_back(z_min);
    node_data_for_output[0][node_index] = x_min - x_min_model_coord;
    node_data_for_output[1][node_index] = y_min - y_min_model_coord;
    node_data_for_output[2][node_index] = z_min - z_min_model_coord;
    node_index += 1;
    // node_x.push_back(x_min); node_y.push_back(y_max); node_z.push_back(z_max);
    node_data_for_output[0][node_index] = x_min - x_min_model_coord;
    node_data_for_output[1][node_index] = y_max - y_max_model_coord;
    node_data_for_output[2][node_index] = z_max - z_max_model_coord;
    node_index += 1;
    // node_x.push_back(x_max); node_y.push_back(y_max); node_z.push_back(z_max);
    node_data_for_output[0][node_index] = x_max - x_max_model_coord;
    node_data_for_output[1][node_index] = y_max - y_max_model_coord;
    node_data_for_output[2][node_index] = z_max - z_max_model_coord;
    node_index += 1;
    // node_x.push_back(x_max); node_y.push_back(y_max); node_z.push_back(z_min);
    node_data_for_output[0][node_index] = x_max - x_max_model_coord;
    node_data_for_output[1][node_index] = y_max - y_max_model_coord;
    node_data_for_output[2][node_index] = z_min - z_min_model_coord;
    node_index += 1;
    // node_x.push_back(x_min); node_y.push_back(y_max); node_z.push_back(z_min);
    node_data_for_output[0][node_index] = x_min - x_min_model_coord;
    node_data_for_output[1][node_index] = y_max - y_max_model_coord;
    node_data_for_output[2][node_index] = z_min - z_min_model_coord;
    node_index += 1;
    */
    exodus_output_for_contact_visualization_.WriteStep(time_current,
                                                       global_data,
                                                       node_data_for_output,
                                                       elem_data_labels_for_output,
                                                       elem_data_for_output,
                                                       derived_elem_data_labels,
                                                       derived_elem_data);
#endif
  }

#ifdef NIMBLE_HAVE_BVH
namespace
{
  void
  WriteContactEntitiesToVTKFile(const std::vector<ContactEntity> &faces,
                                const std::vector<ContactEntity> &nodes,
                                const std::string &prefix,
                                int step) {

    std::stringstream file_name_ss;
    file_name_ss << prefix << step << ".vtk";
    std::ofstream vis_file;
    vis_file.open(file_name_ss.str().c_str());

    // file version and identifier
    vis_file << "# vtk DataFile Version 3.0" << std::endl;

    // header
    vis_file << "Contact entity visualization" << std::endl;

    // file format ASCII | BINARY
    vis_file << "ASCII" << std::endl;

    // dataset structure STRUCTURED_POINTS | STRUCTURED_GRID | UNSTRUCTURED_GRID | POLYDATA | RECTILINEAR_GRID | FIELD
    vis_file << "DATASET UNSTRUCTURED_GRID" << std::endl;

    std::vector<int> vtk_vertices(nodes.size());
    std::vector< std::vector<int> > vtk_triangles(faces.size());

    vis_file << "POINTS " << nodes.size() + 3*faces.size() << " float" << std::endl;

    for (int i=0 ; i<nodes.size() ; i++) {
      vis_file << nodes[i].coord_1_x_ << " " << nodes[i].coord_1_y_ << " " << nodes[i].coord_1_z_ << std::endl;
      vtk_vertices[i] = i;
    }
    int offset = static_cast<int>(vtk_vertices.size());
    for (int i=0 ; i<faces.size() ; i++) {
      ContactEntity const & face = faces[i];
      vis_file << face.coord_1_x_ << " " << face.coord_1_y_ << " " << face.coord_1_z_ << std::endl;
      vis_file << face.coord_2_x_ << " " << face.coord_2_y_ << " " << face.coord_2_z_ << std::endl;
      vis_file << face.coord_3_x_ << " " << face.coord_3_y_ << " " << face.coord_3_z_ << std::endl;
      vtk_triangles[i].push_back(offset + 3*i);
      vtk_triangles[i].push_back(offset + 3*i + 1);
      vtk_triangles[i].push_back(offset + 3*i + 2);
    }

    vis_file << "CELLS " << vtk_vertices.size() + vtk_triangles.size() << " " << 2*vtk_vertices.size() + 4*vtk_triangles.size() << std::endl;
    for(const auto & vtk_vertex : vtk_vertices) {
      vis_file << "1 " << vtk_vertex << std::endl;
    }
    for(const auto & vtk_triangle : vtk_triangles) {
      vis_file << "3 " << vtk_triangle[0] << " " << vtk_triangle[1] << " " << vtk_triangle[2] << std::endl;
    }

    vis_file << "CELL_TYPES " << nodes.size() + faces.size() << std::endl;
    for(unsigned int i=0 ; i<vtk_vertices.size() ; ++i) {
      // cell type 1 is VTK_VERTEX
      vis_file << 1 << std::endl;
    }
    for(unsigned int i=0 ; i<vtk_triangles.size() ; ++i) {
      // cell type 6 is VTK_TRIANGLE
      vis_file << 6 << std::endl;
    }

    vis_file.close();
  }
}

  void
  ContactManager::VisualizeCollisionInfo( const bvh::bvh_tree_26d &faces_tree, const bvh::bvh_tree_26d &nodes_tree,
                                          const bvh::bvh_tree_26d::collision_query_result_type &collision_result,
                                          int step )
  {
    std::vector<ContactEntity> colliding_faces;
    std::vector<ContactEntity> noncolliding_faces;
    for ( auto &&face : contact_faces_ )
    {
      if ( std::count_if( collision_result.begin(), collision_result.end(),
                          [&face]( auto &&pair ){ return pair.first == face.contact_entity_global_id(); } ) )
      {
        colliding_faces.push_back(face);
      } else {
        noncolliding_faces.push_back(face);
      }
    }

    std::vector<ContactEntity> colliding_nodes;
    std::vector<ContactEntity> noncolliding_nodes;
    for ( auto &&node : contact_nodes_ )
    {
      if ( std::count_if( collision_result.begin(), collision_result.end(),
                          [&node]( auto &&pair ){ return pair.first == node.contact_entity_global_id(); } ) )
      {
        colliding_nodes.push_back(node);
      } else {
        noncolliding_nodes.push_back(node);
      }
    }

#ifdef BVH_USE_VTK
    // Visualize bvh trees
    std::stringstream tree_faces_file_name_ss;
    tree_faces_file_name_ss << "bvh_tree_faces_" << step << ".vtp";
    std::ofstream tree_faces_vis_file(tree_faces_file_name_ss.str().c_str());

    bvh::vis::write_bvh( tree_faces_vis_file, faces_tree );

    tree_faces_vis_file.close();

    std::stringstream tree_nodes_file_name_ss;
    tree_nodes_file_name_ss << "bvh_tree_nodes_" << step << ".vtp";
    std::ofstream tree_nodes_vis_file(tree_nodes_file_name_ss.str().c_str());

    bvh::vis::write_bvh( tree_nodes_vis_file, nodes_tree );

    tree_nodes_vis_file.close();
#endif
  }
#endif

  void
  ContactManager::BruteForceBoxIntersectionSearch(std::vector<ContactEntity> const & nodes,
                                                  std::vector<ContactEntity> const & triangles) {

    for (int i_node = 0 ; i_node < nodes.size() ; i_node++) {
      for (int i_tri = 0 ; i_tri < triangles.size() ; i_tri++) {

        //ContactEntity.get_x_min()

      }
    }
  }

  void
  ContactManager::ClosestPointProjection(std::vector<ContactEntity> const & nodes,
                                         std::vector<ContactEntity> const & triangles,
                                         std::vector<ContactEntity::vertex>& closest_points,
                                         std::vector<ContactManager::PROJECTION_TYPE>& projection_types) {

    // Wolfgang Heidrich, 2005, Computing the Barycentric Coordinates of a Projected Point, Journal of Graphics Tools, pp 9-12, 10(3).

    double tol = 1.0e-16;

    for (unsigned int i_proj=0 ; i_proj<nodes.size() ; i_proj++) {

      ContactEntity const & node = nodes[i_proj];
      ContactEntity const & tri = triangles[i_proj];

      double p[3];
      p[0] = node.coord_1_x_;
      p[1] = node.coord_1_y_;
      p[2] = node.coord_1_z_;

      double p1[3];
      p1[0] = tri.coord_1_x_;
      p1[1] = tri.coord_1_y_;
      p1[2] = tri.coord_1_z_;

      double p2[3];
      p2[0] = tri.coord_2_x_;
      p2[1] = tri.coord_2_y_;
      p2[2] = tri.coord_2_z_;

      double p3[3];
      p3[0] = tri.coord_3_x_;
      p3[1] = tri.coord_3_y_;
      p3[2] = tri.coord_3_z_;

      double u[3], v[3], w[3];
      for (int i=0 ; i<3 ; i++) {
        u[i] = p2[i] - p1[i];
        v[i] = p3[i] - p1[i];
        w[i] = p[i]  - p1[i];
      }

      double n[3];
      CrossProduct(u, v, n);

      double n_squared = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];

      double cross[3];
      CrossProduct(u, w, cross);
      double gamma = (cross[0]*n[0] + cross[1]*n[1] + cross[2]*n[2]) / n_squared;

      CrossProduct(w, v, cross);
      double beta = (cross[0]*n[0] + cross[1]*n[1] + cross[2]*n[2]) / n_squared;

      double alpha = 1 - gamma - beta;

      bool alpha_is_zero, alpha_in_range;
      (alpha > -tol && alpha < tol) ? alpha_is_zero = true : alpha_is_zero = false;
      (alpha > -tol && alpha < 1.0 + tol) ? alpha_in_range = true : alpha_in_range = false;

      bool beta_is_zero, beta_in_range;
      (beta > -tol && beta < tol) ? beta_is_zero = true : beta_is_zero = false;
      (beta > -tol && beta < 1.0 + tol) ? beta_in_range = true : beta_in_range = false;

      bool gamma_is_zero, gamma_in_range;
      (gamma > -tol && gamma < tol) ? gamma_is_zero = true : gamma_is_zero = false;
      (gamma > -tol && gamma < 1.0 + tol) ? gamma_in_range = true : gamma_in_range = false;

      if (alpha_in_range && beta_in_range && gamma_in_range) {
        closest_points[i_proj].coords_[0] = alpha*p1[0] + beta*p2[0] + gamma*p3[0];
        closest_points[i_proj].coords_[1] = alpha*p1[1] + beta*p2[1] + gamma*p3[1];
        closest_points[i_proj].coords_[2] = alpha*p1[2] + beta*p2[2] + gamma*p3[2];
        if (alpha_is_zero || beta_is_zero || gamma_is_zero) {
          projection_types[i_proj] = PROJECTION_TYPE::NODE_OR_EDGE;
        }
        else {
          projection_types[i_proj] = PROJECTION_TYPE::FACE;
        }
      }
      else {

        projection_types[i_proj] = PROJECTION_TYPE::NODE_OR_EDGE;

        double x = p1[0] - p[0];
        double y = p1[1] - p[1];
        double z = p1[2] - p[2];
        double distance_squared = x*x + y*y + z*z;

        double min_distance_squared = distance_squared;
        double min_distance_squared_t;
        int min_case = 1;

        x = p2[0] - p[0];
        y = p2[1] - p[1];
        z = p2[2] - p[2];
        distance_squared = x*x + y*y + z*z;
        if (distance_squared < min_distance_squared) {
          min_distance_squared = distance_squared;
          min_case = 2;
        }

        x = p3[0] - p[0];
        y = p3[1] - p[1];
        z = p3[2] - p[2];
        distance_squared = x*x + y*y + z*z;
        if (distance_squared < min_distance_squared) {
          min_distance_squared = distance_squared;
          min_case = 3;
        }

        double t = PointEdgeClosestPointFindT(p1, p2, p);
        if (t > 0.0 && t < 1.0) {
          distance_squared = PointEdgeClosestPointFindDistanceSquared(p1, p2, p, t);
          if (distance_squared < min_distance_squared) {
            min_distance_squared = distance_squared;
            min_distance_squared_t = t;
            min_case = 4;
          }
        }

        t = PointEdgeClosestPointFindT(p2, p3, p);
        if (t > 0.0 && t < 1.0) {
          distance_squared = PointEdgeClosestPointFindDistanceSquared(p2, p3, p, t);
          if (distance_squared < min_distance_squared) {
            min_distance_squared = distance_squared;
            min_distance_squared_t = t;
            min_case = 5;
          }
        }

        t = PointEdgeClosestPointFindT(p3, p1, p);
        if (t > 0.0 && t < 1.0) {
          distance_squared = PointEdgeClosestPointFindDistanceSquared(p3, p1, p, t);
          if (distance_squared < min_distance_squared) {
            min_distance_squared = distance_squared;
            min_distance_squared_t = t;
            min_case = 6;
          }
        }

        switch (min_case) {
        case 1:
          closest_points[i_proj].coords_[0] = p1[0];
          closest_points[i_proj].coords_[1] = p1[1];
          closest_points[i_proj].coords_[2] = p1[2];
          break;
        case 2:
          closest_points[i_proj].coords_[0] = p2[0];
          closest_points[i_proj].coords_[1] = p2[1];
          closest_points[i_proj].coords_[2] = p2[2];
          break;
        case 3:
          closest_points[i_proj].coords_[0] = p3[0];
          closest_points[i_proj].coords_[1] = p3[1];
          closest_points[i_proj].coords_[2] = p3[2];
          break;
        case 4:
          closest_points[i_proj].coords_[0] = p1[0] + (p2[0] - p1[0])*min_distance_squared_t;
          closest_points[i_proj].coords_[1] = p1[1] + (p2[1] - p1[1])*min_distance_squared_t;
          closest_points[i_proj].coords_[2] = p1[2] + (p2[2] - p1[2])*min_distance_squared_t;
          break;
        case 5:
          closest_points[i_proj].coords_[0] = p2[0] + (p3[0] - p2[0])*min_distance_squared_t;
          closest_points[i_proj].coords_[1] = p2[1] + (p3[1] - p2[1])*min_distance_squared_t;
          closest_points[i_proj].coords_[2] = p2[2] + (p3[2] - p2[2])*min_distance_squared_t;
          break;
        case 6:
          closest_points[i_proj].coords_[0] = p3[0] + (p1[0] - p3[0])*min_distance_squared_t;
          closest_points[i_proj].coords_[1] = p3[1] + (p1[1] - p3[1])*min_distance_squared_t;
          closest_points[i_proj].coords_[2] = p3[2] + (p1[2] - p3[2])*min_distance_squared_t;
          break;
        }
      }
    }
  }

  void ContactManager::ComputeContactForce(int step, bool debug_output) {

    if (penalty_parameter_ <= 0.0) {
      throw std::logic_error("\nError in ComputeContactForce(), invalid penalty_parameter.\n");
    }

    double background_grid_cell_size = BoundingBoxAverageCharacteristicLengthOverAllRanks();
    // std::cout << "DEBUGGING background_grid_cell_size " << background_grid_cell_size << std::endl;
    // DJL PARALLEL CONTACT  processCollision(coord_.data(), coord_.size(), background_grid_cell_size);

    // kokkos_view::extent(), accessor kokkos_view::(i, 0) x coordinte of ith entry in list

    // ANTONIO, HERE IS THE COMMENTED OUT LINE
    // std::vector<int> exchange_members = getExchangeMembers(coord_d_,
    //                                                        background_grid_cell_size,
    //                                                        MPI_COMM_WORLD,
    //                                                        0,
    //                                                        1,
    //                                                        2);

#ifdef NIMBLE_HAVE_BVH

    // Force kdop recomputation
    for ( auto &&ent : contact_faces_ )
    {
      ent.RecomputeKdop();
    }
    for ( auto &&ent : contact_nodes_ )
    {
      ent.RecomputeKdop();
    }

    // Construct the BVH trees for narrowphase
    // Can also be done by bvh internally, but using thsese for visualization
    auto faces_tree = bvh::bvh_tree_26d{ contact_faces_.begin(), contact_faces_.end() };
    auto nodes_tree = bvh::bvh_tree_26d{ contact_nodes_.begin(), contact_nodes_.end() };

    // For now, if leaves collide say the primitves also collide
    auto collision_results = bvh::narrowphase_tree( faces_tree, nodes_tree,
                                                    contact_faces_, contact_nodes_,
                                                    []( const auto &, const auto & ){ return true; } );

    if ( debug_output )
    {
      VisualizeCollisionInfo(faces_tree, nodes_tree, collision_results, step);
    }

#endif

    // DJL
    // entities are stored in contact_nodes_d, contact_nodes_h
    // 1) box box search
    // 2) node-face projection
    // 3) culling
    // 4) enforcement

#ifdef NIMBLE_HAVE_KOKKOS
  contact_interface->ComputeContact(contact_nodes_d_, contact_faces_d_, force_d_, penalty_parameter_);
#endif
  }
#if defined(NIMBLE_HAVE_MPI) && defined(NIMBLE_HAVE_BVH)
  namespace
  {
    bvh::span< ContactEntity, bvh::dynamic_extent() >
    span_for_decomp(int od, int od_factor, std::vector<ContactEntity> &_vec)
    {
      auto size = _vec.size();

      auto count = size / static_cast<std::size_t>(od_factor);
      auto rem = size % static_cast<std::size_t>(od_factor);

      bvh::span< ContactEntity, bvh::dynamic_extent() > sp( &_vec[0], _vec.size() );

      if (od < rem) {
        return sp.subspan(od * (count + 1), count + 1);
      } else {
        return sp.subspan(rem + od * count, count);
      }
    }

    auto
    build_patch(ContactManager::patch_collection &patch_collection,
                               std::vector<ContactEntity> &entities,
                               int od_factor) {
      // Build patches (averages centroids and builds kdops)
      // od_factor is used for overdecomposition
      // TODO: Use tree splitting metric for splitting patches rather than indices

      // Force kdop recomputation
      for ( auto &&ent : entities )
      {
        ent.RecomputeKdop();
      }

      int rank = bvh::vt::context::current()->rank();

      std::vector< bvh::patch< ContactEntity > >  patches_vec;
      patches_vec.reserve( static_cast< std::size_t >( od_factor ) );
      for ( int i = 0; i < od_factor; ++i )
      {
        patches_vec.emplace_back(i + rank * od_factor, span_for_decomp(i, od_factor, entities));
      }

      return bvh::vt::vt_collection_data(patches_vec, patch_collection);
    }
  }

  void
  ContactManager::ComputeParallelContactForce(int step, bool is_output_step, bool visualize) {

    if (penalty_parameter_ <= 0.0) {
      throw std::logic_error("\nError in ComputeParallelContactForce(), invalid penalty_parameter.\n");
    }
    auto od_factor = static_cast<int>(dicing_factor_);

    auto face_patches_future = build_patch(face_patch_collection_, contact_faces_, od_factor);
    auto node_patches_future = build_patch(node_patch_collection_, contact_nodes_, od_factor);

    auto &world = collision_world_;

    // As soon as the patches have been transferred to a VT collection, start using them to build the trees
    auto trees_future = node_patches_future.then([&world](auto &&tree_nodes){
      return world.build_trees(std::forward<decltype(tree_nodes)>(tree_nodes));
    } );

    // As soon as the trees have completed and the face patches are available, starting to collision
    // tests between the trees and patches
    auto bpr = bvh::vt::when_all(trees_future, face_patches_future).then([&world](auto &&tup) {
      auto face_patches = std::get<1>(std::forward<decltype(tup)>(tup));

      return world.find_collisions(face_patches);
    });

    // When collisions have been processed, transfer back to a standard vector per rank from a VT collection
    auto collision_result_future = bpr.then([](auto &&results_collection){
      return bvh::vt::vt_collection_to_mpi<std::vector<bvh::bvh_tree_26d::collision_query_result_type>>(results_collection);
    });

    std::vector< bvh::bvh_tree_26d::collision_query_result_type > results_vec(dicing_factor_);
    bvh::vt::debug( "{}: ============begin get results\n", ::vt::theContext()->getNode() );
    results_vec = collision_result_future.get();
    bvh::vt::debug( "{}: ============end get results\n", ::vt::theContext()->getNode() );

#if 1

    bvh::bvh_tree_26d::collision_query_result_type collision_result;

    for ( const auto &r : results_vec )
      collision_result.pairs.insert( collision_result.end(), r.begin(), r.end() );

    unsigned long long total_num_collisions = 0;
    for ( auto &&v : results_vec)
      total_num_collisions += v.size();

    int rank = bvh::vt::context::current()->rank();

    if ( is_output_step )
    {
      unsigned long long ncollisions = total_num_collisions;
      MPI_Reduce(&ncollisions, &total_num_collisions, 1, MPI_UNSIGNED_LONG_LONG,
          MPI_SUM, 0, MPI_COMM_WORLD );

      if (rank == 0)
        bvh::vt::print( "num collisions: {}\n", total_num_collisions );
      //VisualizeCollisionInfo(faces_tree, nodes_tree, collision_results, step);
      std::vector<ContactEntity> colliding_faces;
      std::vector<ContactEntity> noncolliding_faces;
      for ( auto &&face : contact_faces_ )
      {
        if ( std::count_if( collision_result.begin(), collision_result.end(),
                            [&face]( auto &&pair ){ return pair.first == face.contact_entity_global_id(); } ) )
        {
          colliding_faces.push_back(face);
        } else {
          noncolliding_faces.push_back(face);
        }
      }

      std::vector<ContactEntity> colliding_nodes;
      std::vector<ContactEntity> noncolliding_nodes;
      for ( auto &&node : contact_nodes_ )
      {
        if ( std::count_if( collision_result.begin(), collision_result.end(),
                            [&node]( auto &&pair ){ return pair.first == node.contact_entity_global_id(); } ) )
        {
          colliding_nodes.push_back(node);
        } else {
          noncolliding_nodes.push_back(node);
        }
      }

      std::stringstream colliding_out_name;
      colliding_out_name << "contact_entities_colliding." << rank << '.';

      std::stringstream noncolliding_out_name;
      noncolliding_out_name << "contact_entities_noncolliding_." << rank << '.';

      WriteContactEntitiesToVTKFile(colliding_faces, colliding_nodes, colliding_out_name.str(), step);
      WriteContactEntitiesToVTKFile(noncolliding_faces, noncolliding_nodes, noncolliding_out_name.str(), step);
    }
#endif
  }
#endif

}
