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

#include "nimble_contact.h"

#ifdef NIMBLE_HAVE_BVH
  #include <bvh/broadphase.hpp>
  #include <bvh/narrowphase.hpp>
  #include <bvh/tree.hpp>
  #include <bvh/vis/vis_bvh.hpp>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>

namespace nimble {

#ifdef NIMBLE_HAVE_EXTRAS
  std::string GTKProjectionTypeToString(short int val) {
    gtk::ProjectionType type = static_cast<gtk::ProjectionType>(val);
    std::string type_str("Unknown");
    switch (type) {
    case gtk::NODE_PROJECTION:
      type_str = "NODE_PROJECTION";
      break;
    case gtk::EDGE_PROJECTION:
      type_str = "EDGE_PROJECTION";
      break;
    case gtk::FACE_PROJECTION:
      type_str = "FACE_PROJECTION";
      break;
    case gtk::NUM_PROJ_TYPE:
      type_str = "NUM_PROJ_TYPE";
      break;
    case gtk::NULL_PROJECTION:
      type_str = "NULL_PROJECTION";
      break;
    default:
      throw std::logic_error("\nError in GTKProjectionTypeToString(), unrecognized ProjectionType.\n");
      break;
    return type_str;
    }
  }
#endif

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

  std::vector< std::vector<int> >
  ContactManager::SkinBlocks(GenesisMesh const & mesh,
                             std::vector<int> block_ids) {

    std::map< std::vector<int>, std::vector<int> > faces;
    std::map< std::vector<int>, std::vector<int> >::iterator face_it;

    for (auto & block_id : block_ids) {

      int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      int num_node_per_elem = mesh.GetNumNodesPerElement(block_id);
      const int * const conn = mesh.GetConnectivity(block_id);
      int conn_index = 0;

      // key is sorted node list for a face
      int num_node_per_face = 4;
      std::vector<int> key(num_node_per_face);
      // value is count and unsorted node list for a face
      std::vector<int> value(num_node_per_face + 1);

      for (int i_elem = 0 ; i_elem < num_elem_in_block ; i_elem++) {

        // Examine each face, following the Exodus node-ordering convention

        // face 0: 0 1 5 4
        value[0] = 1;
        key[0] = value[1] = conn[conn_index + 0];
        key[1] = value[2] = conn[conn_index + 1];
        key[2] = value[3] = conn[conn_index + 5];
        key[3] = value[4] = conn[conn_index + 4];
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
        std::sort(key.begin(), key.end());
        face_it = faces.find(key);
        if (face_it == faces.end())
          faces[key] = value;
        else
          face_it->second[0] += 1;

        conn_index += num_node_per_elem;
      }
    }

    std::vector< std::vector<int> > skin_faces;
    for (auto face : faces) {
      if (face.second[0] == 1) {
        std::vector<int> skin_face;
        for (int i=0 ; i<face.second.size()-1 ; i++) {
          skin_face.push_back(face.second.at(i+1));
        }
        skin_faces.push_back(skin_face);
      }
      else if (face.second[0] != 2) {
        throw std::logic_error("Error in mesh skinning routine, face found more than two times!\n");
      }
    }

    return skin_faces;
  }

  void
  ContactManager::CreateContactEntities(GenesisMesh const & mesh,
                                        std::vector<int> master_block_ids,
                                        std::vector<int> slave_block_ids) {
    contact_enabled_ = true;

    const double* coord_x = mesh.GetCoordinatesX();
    const double* coord_y = mesh.GetCoordinatesY();
    const double* coord_z = mesh.GetCoordinatesZ();

    // find all the element faces on the master and slave contact blocks
    std::vector< std::vector<int> > master_skin_faces = SkinBlocks(mesh, master_block_ids);
    std::vector< std::vector<int> > slave_skin_faces = SkinBlocks(mesh, slave_block_ids);

    // construct containers for the subset of the model that is involved with contact
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

    std::map<int, int> genesis_mesh_node_id_to_contact_vector_id;
    for (unsigned int i_node=0 ; i_node<node_ids_.size() ; ++i_node) {
      genesis_mesh_node_id_to_contact_vector_id[node_ids_[i_node]] = i_node;
    }

    // replace the node ids that correspond to the genesis mesh
    // with node ids that correspond to the contact submodel
    for (unsigned int i_face=0 ; i_face<master_skin_faces.size() ; ++i_face) {
      for (unsigned int i_node=0 ; i_node<master_skin_faces[i_face].size() ; ++i_node) {
        int genesis_mesh_node_id = master_skin_faces[i_face][i_node];
        master_skin_faces[i_face][i_node] = genesis_mesh_node_id_to_contact_vector_id.at(genesis_mesh_node_id);
      }
    }
    for (unsigned int i_face=0 ; i_face<slave_skin_faces.size() ; ++i_face) {
      for (unsigned int i_node=0 ; i_node<slave_skin_faces[i_face].size() ; ++i_node) {
        int genesis_mesh_node_id = slave_skin_faces[i_face][i_node];
        slave_skin_faces[i_face][i_node] = genesis_mesh_node_id_to_contact_vector_id.at(genesis_mesh_node_id);
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
    std::vector<int> slave_node_ids;
    std::map<int, double> slave_node_char_lens;
    for (auto const & face : slave_skin_faces) {

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

      for (auto const & node_id : face) {
        std::vector<int>::iterator it = std::find(slave_node_ids.begin(), slave_node_ids.end(), node_id);
        if (it == slave_node_ids.end()) {
          slave_node_ids.push_back(node_id);
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

    contact_nodes_.resize(slave_node_ids.size());
    contact_faces_.resize(4*master_skin_faces.size());
    CreateContactNodesAndFaces(master_skin_faces, slave_node_ids, slave_node_char_lens, contact_nodes_, contact_faces_);

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
    CreateContactNodesAndFaces(master_skin_faces, slave_node_ids, slave_node_char_lens, contact_nodes_h_, contact_faces_h_);

    Kokkos::resize(contact_nodes_d_, slave_node_ids.size());
    Kokkos::resize(contact_faces_d_, 4*master_skin_faces.size());
    Kokkos::deep_copy(contact_nodes_d_, contact_nodes_h_);
    Kokkos::deep_copy(contact_faces_d_, contact_faces_h_);
#endif

    std::cout << "Contact initialization:" << std::endl;
    std::cout << "  number of triangular contact facets (master blocks): " << contact_faces_.size() << std::endl;
    std::cout << "  number of contact nodes (slave blocks): " << contact_nodes_.size() << "\n" << std::endl;
  }

  template <typename ArgT>
  void ContactManager::CreateContactNodesAndFaces(std::vector< std::vector<int> > const & master_skin_faces,
                                                  std::vector<int> const & slave_node_ids,
                                                  std::map<int, double> const & slave_node_char_lens,
                                                  ArgT& contact_nodes,
                                                  ArgT& contact_faces) const {

    int contact_entity_id = 0;
    int index = 0;

    // convert master faces to trangular facets
    for (auto face : master_skin_faces) {

      int num_nodes_in_face = static_cast<int>(face.size());
      if (num_nodes_in_face != 4) {
        throw std::logic_error("\nError in ContactManager::CreateContactEntities(), invalid number of face nodes.\n");
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
      int node_id_1, node_id_2;

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
      contact_faces[index++] = ContactEntity(ContactEntity::TRIANGLE,
                                             contact_entity_id++,
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
      contact_faces[index++] = ContactEntity(ContactEntity::TRIANGLE,
                                             contact_entity_id++,
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
      contact_faces[index++] = ContactEntity(ContactEntity::TRIANGLE,
                                             contact_entity_id++,
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
      contact_faces[index++] = ContactEntity(ContactEntity::TRIANGLE,
                                             contact_entity_id++,
                                             model_coord,
                                             characteristic_length,
                                             node_id_1,
                                             node_id_2,
                                             node_ids_for_fictitious_node);
    }

    // Slave node entities
    index = 0;
    for (auto const & node_id : slave_node_ids) {
      double model_coord[3];
      for (int i=0 ; i<3 ; ++i) {
        model_coord[i] = coord_[3*node_id+i];
      }
      double characteristic_length = slave_node_char_lens.at(node_id);
      contact_nodes[index++] = ContactEntity(ContactEntity::NODE,
                                             contact_entity_id++,
                                             model_coord,
                                             characteristic_length,
                                             node_id);
    }
  }

  void
  ContactManager::WriteContactEntitiesToVTKFile(int step) {
    WriteContactEntitiesToVTKFile(contact_faces_, contact_nodes_, "contact_entities_", step);
  }

  void
  ContactManager::WriteContactEntitiesToVTKFile(const std::vector<ContactEntity> &faces,
                                                const std::vector<ContactEntity> &nodes,
                                                const std::string &_prefix, int step) {

    std::stringstream file_name_ss;
    file_name_ss << _prefix << step << ".vtk";
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

#ifdef NIMBLE_HAVE_BVH
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
                          [&face]( auto &&pair ){ return pair.first == face.global_id(); } ) )
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
                          [&node]( auto &&pair ){ return pair.first == node.global_id(); } ) )
      {
        colliding_nodes.push_back(node);
      } else {
        noncolliding_nodes.push_back(node);
      }
    }

    WriteContactEntitiesToVTKFile(colliding_faces, colliding_nodes, "contact_entities_colliding_", step);
    WriteContactEntitiesToVTKFile(noncolliding_faces, noncolliding_nodes, "contact_entities_noncolliding_", step);

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
  ContactManager::ComputeContactForce(int step, bool debug_output) {

    if (penalty_parameter_ <= 0.0) {
      throw std::logic_error("\nError in ComputeContactForce(), invalid penalty_parameter.\n");
    }

#ifdef NIMBLE_HAVE_BVH
    // Construct the BVH trees for narrowphase
    // Can also be done by bvh internall,y but using thsese for visualization
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

#ifdef NIMBLE_HAVE_EXTRAS

    stk::search::CollisionList<nimble_kokkos::kokkos_device_execution_space> collision_list("contact_collision_list");
    stk::search::MortonLBVHSearch_Timers timers;

    stk::search::mas_aabb_tree_loader<double, nimble_kokkos::kokkos_device_execution_space>
    contact_nodes_tree_loader(contact_nodes_search_tree_, contact_nodes_.size());

    int num_contact_nodes = contact_nodes_.size();
    DeviceContactEntityArrayView contact_nodes = contact_nodes_d_;
    Kokkos::parallel_for("Load contact nodes search tree",
                         num_contact_nodes,
                         KOKKOS_LAMBDA(const int i_contact_node) {
      double min_x = contact_nodes(i_contact_node).get_x_min();
      double max_x = contact_nodes(i_contact_node).get_x_max();
      double min_y = contact_nodes(i_contact_node).get_y_min();
      double max_y = contact_nodes(i_contact_node).get_y_max();
      double min_z = contact_nodes(i_contact_node).get_z_min();
      double max_z = contact_nodes(i_contact_node).get_z_max();
      contact_nodes_tree_loader.set_box(i_contact_node, min_x, max_x, min_y, max_y, min_z, max_z);
    });

    stk::search::mas_aabb_tree_loader<double, nimble_kokkos::kokkos_device_execution_space>
    contact_faces_tree_loader(contact_faces_search_tree_, contact_faces_.size());

    int num_contact_faces = contact_faces_.size();
    DeviceContactEntityArrayView contact_faces = contact_faces_d_;
    Kokkos::parallel_for("Load contact faces search tree",
                         num_contact_faces,
                         KOKKOS_LAMBDA(const int i_contact_face) {
      double min_x = contact_faces(i_contact_face).get_x_min();
      double max_x = contact_faces(i_contact_face).get_x_max();
      double min_y = contact_faces(i_contact_face).get_y_min();
      double max_y = contact_faces(i_contact_face).get_y_max();
      double min_z = contact_faces(i_contact_face).get_z_min();
      double max_z = contact_faces(i_contact_face).get_z_max();
      contact_faces_tree_loader.set_box(i_contact_face, min_x, max_x, min_y, max_y, min_z, max_z);
    });

    stk::search::TimedMortonLBVHSearch<double, nimble_kokkos::kokkos_device_execution_space>(contact_nodes_search_tree_,
                                                                                             contact_faces_search_tree_,
                                                                                             collision_list,
                                                                                             timers);

    auto num_collisions = collision_list.get_num_collisions();
    gtk::PointsView<nimble_kokkos::kokkos_device_execution_space> points("points", num_collisions);
    gtk::TrianglesView<nimble_kokkos::kokkos_device_execution_space> triangles("triangles", num_collisions);
    gtk::PointsView<nimble_kokkos::kokkos_device_execution_space> closest_points("closest_points", num_collisions);

    // TODO CARRY OUT ENFORCEMENT ON THE DEVICE
    // For now, copy everything to the host (ugh)
    collision_list.sync_from_device();
    Kokkos::deep_copy(contact_nodes_h_, contact_nodes_d_);
    Kokkos::deep_copy(contact_faces_h_, contact_faces_d_);

    std::map<int, std::vector<int> > collision_indices_for_each_contact_node;
    for (int i=0 ; i<num_collisions ; i++) {
      int node_index = collision_list.hm_data(i,0);
      collision_indices_for_each_contact_node[node_index] = std::vector<int>();
    }

    for (int i_collision=0 ; i_collision<num_collisions ; i_collision++) {
      int contact_node_index = collision_list.hm_data(i_collision,0);
      int contact_face_index = collision_list.hm_data(i_collision,1);
      ContactEntity const & node = contact_nodes_h_[contact_node_index];
      ContactEntity const & face = contact_faces_h_[contact_face_index];
      points.setHostPointValue(i_collision, node.coord_1_x_, node.coord_1_y_, node.coord_1_z_);
      triangles.setHostVertexValues(i_collision,
                                    mtk::Vec3<double>(face.coord_1_x_, face.coord_1_y_, face.coord_1_z_),
                                    mtk::Vec3<double>(face.coord_2_x_, face.coord_2_y_, face.coord_2_z_),
                                    mtk::Vec3<double>(face.coord_3_x_, face.coord_3_y_, face.coord_3_z_));
      collision_indices_for_each_contact_node[contact_node_index].push_back(i_collision);
    }

    // TEMP COPY TO DEVICE
    points.inhaleHostData();
    triangles.inhaleHostData();

    constexpr bool save_projection_types_computed = true;
    Kokkos::View<short *, nimble_kokkos::kokkos_device_execution_space> proj_types_returned_d;
    if (save_projection_types_computed) {
      Kokkos::resize(proj_types_returned_d, num_collisions);
    }

    Kokkos::Impl::Timer timer;
    gtk::ComputeProjections<nimble_kokkos::kokkos_device_execution_space, save_projection_types_computed> projection(points,
                                                                                                                     triangles,
                                                                                                                     closest_points,
                                                                                                                     proj_types_returned_d);

    // TEMP COPY TO HOST
    Kokkos::deep_copy(closest_points.m_hostData, closest_points.m_data);

    Kokkos::View<short *, nimble_kokkos::kokkos_device_execution_space>::HostMirror proj_types_returned_h;
    Kokkos::resize(proj_types_returned_h, proj_types_returned_d.extent(0));
    Kokkos::deep_copy(proj_types_returned_h, proj_types_returned_d);

    std::vector<int> contact_collisions;

    for (auto const & entry : collision_indices_for_each_contact_node) {

      std::vector<int> const & collision_indices = entry.second;
      int first_collision_index = collision_indices[0];
      double pt_value[3];
      //points.getPointValue(first_collision_index, pt_value);
#ifdef NIMBLE_NVIDIA_BUILD
      pt_value[0] = points.m_hostData(first_collision_index, 0);
      pt_value[1] = points.m_hostData(first_collision_index, 1);
      pt_value[2] = points.m_hostData(first_collision_index, 2);
#else
      pt_value[0] = points.m_hostData(first_collision_index).X();
      pt_value[1] = points.m_hostData(first_collision_index).Y();
      pt_value[2] = points.m_hostData(first_collision_index).Z();
#endif

      int contact_node_index = collision_list.hm_data(first_collision_index,0);

      std::vector<int> face_projections;
      std::vector<int> node_and_edge_projections;

      for (auto const & collision_index : collision_indices) {

        int contact_node_index_check = collision_list.hm_data(collision_index,0);
        int contact_face_index = collision_list.hm_data(collision_index,1);
        if (contact_node_index_check != contact_node_index) {
          throw std::logic_error("\nError in ComputeContactForce(), bad node index.\n");
        }

        double closest_pt_value[3], proj_vector[3];
        //closest_points.getPointValue(collision_index, closest_pt_value);
#ifdef NIMBLE_NVIDIA_BUILD
        closest_pt_value[0] = closest_points.m_hostData(collision_index, 0);
        closest_pt_value[1] = closest_points.m_hostData(collision_index, 2);
        closest_pt_value[2] = closest_points.m_hostData(collision_index, 3);
#else
        closest_pt_value[0] = closest_points.m_hostData(collision_index).X();
        closest_pt_value[1] = closest_points.m_hostData(collision_index).Y();
        closest_pt_value[2] = closest_points.m_hostData(collision_index).Z();
#endif

        for (int i=0 ; i<3 ; i++) {
          proj_vector[i] = closest_pt_value[i] - pt_value[i];
        }

        double tri_node_1[3], tri_node_2[3], tri_node_3[3];
#ifdef NIMBLE_NVIDIA_BUILD
        //triangles.getVertexValue(collision_index, 0, tri_node_1);
        tri_node_1[0] = triangles.m_hostData(collision_index, 0, 0);
        tri_node_1[1] = triangles.m_hostData(collision_index, 0, 1);
        tri_node_1[2] = triangles.m_hostData(collision_index, 0, 2);
        //triangles.getVertexValue(collision_index, 1, tri_node_2);
        tri_node_2[0] = triangles.m_hostData(collision_index, 1, 0);
        tri_node_2[1] = triangles.m_hostData(collision_index, 1, 1);
        tri_node_2[2] = triangles.m_hostData(collision_index, 1, 2);
        //triangles.getVertexValue(collision_index, 2, tri_node_3);
        tri_node_3[0] = triangles.m_hostData(collision_index, 2, 0);
        tri_node_3[1] = triangles.m_hostData(collision_index, 2, 1);
        tri_node_3[2] = triangles.m_hostData(collision_index, 2, 2);
#else
        //triangles.getVertexValue(collision_index, 0, tri_node_1);
        tri_node_1[0] = triangles.m_hostData(collision_index).GetNode(0)[0];
        tri_node_1[1] = triangles.m_hostData(collision_index).GetNode(0)[1];
        tri_node_1[2] = triangles.m_hostData(collision_index).GetNode(0)[2];
        //triangles.getVertexValue(collision_index, 1, tri_node_2);
        tri_node_2[0] = triangles.m_hostData(collision_index).GetNode(1)[0];
        tri_node_2[1] = triangles.m_hostData(collision_index).GetNode(1)[1];
        tri_node_2[2] = triangles.m_hostData(collision_index).GetNode(1)[2];
        //triangles.getVertexValue(collision_index, 2, tri_node_3);
        tri_node_3[0] = triangles.m_hostData(collision_index).GetNode(2)[0];
        tri_node_3[1] = triangles.m_hostData(collision_index).GetNode(2)[1];
        tri_node_3[2] = triangles.m_hostData(collision_index).GetNode(2)[2];
#endif
        double tri_edge_1[3], tri_edge_2[3], tri_normal[3];
        for (int i=0 ; i<3 ; i++) {
          tri_edge_1[i] = tri_node_2[i] - tri_node_1[i];
          tri_edge_2[i] = tri_node_3[i] - tri_node_2[i];
        }
        tri_normal[0] = tri_edge_1[1]*tri_edge_2[2] - tri_edge_1[2]*tri_edge_2[1];
        tri_normal[1] = tri_edge_1[0]*tri_edge_2[2] - tri_edge_1[2]*tri_edge_2[0];
        tri_normal[2] = tri_edge_1[0]*tri_edge_2[1] - tri_edge_1[1]*tri_edge_2[0];

        double dot = proj_vector[0]*tri_normal[0] + proj_vector[1]*tri_normal[1] + proj_vector[2]*tri_normal[2];

        bool penetration = (dot >= 0.0);

        if (penetration) {
          gtk::ProjectionType type = static_cast<gtk::ProjectionType>(proj_types_returned_h[collision_index]);
          switch (type) {
          case gtk::NODE_PROJECTION:
            node_and_edge_projections.push_back(collision_index);
            break;
          case gtk::EDGE_PROJECTION:
            node_and_edge_projections.push_back(collision_index);
            break;
          case gtk::FACE_PROJECTION:
            face_projections.push_back(collision_index);
            break;
          default:
            throw std::logic_error("\nError in ComputeContactForce(), unrecognized ProjectionType.\n");
            break;
          }
        }
      }

      bool has_contact = false;
      int contact_index(-1);

      if (face_projections.size() == 0 && node_and_edge_projections.size() > 0) {
        // if there are no face projections, use the edge or node projection with the shortest projection distance
        has_contact = true;
        double min_distance(std::numeric_limits<double>::max());
        double closest_pt_value[3], proj_vector[3], distance;
        for (auto & collision_index : node_and_edge_projections) {
          //closest_points.getPointValue(collision_index, closest_pt_value);
#ifdef NIMBLE_NVIDIA_BUILD
          closest_pt_value[0] = closest_points.m_hostData(collision_index, 0);
          closest_pt_value[1] = closest_points.m_hostData(collision_index, 1);
          closest_pt_value[2] = closest_points.m_hostData(collision_index, 2);
#else
          closest_pt_value[0] = closest_points.m_hostData(collision_index).X();
          closest_pt_value[1] = closest_points.m_hostData(collision_index).Y();
          closest_pt_value[2] = closest_points.m_hostData(collision_index).Z();
#endif
          for (int i=0 ; i<3 ; i++) {
            proj_vector[i] = closest_pt_value[i] - pt_value[i];
          }
          distance = std::sqrt(proj_vector[0]*proj_vector[0] + proj_vector[1]*proj_vector[1] + proj_vector[2]*proj_vector[2]);
          if (distance < min_distance) {
            contact_index = collision_index;
            min_distance = distance;
          }
        }
      }
      else if (face_projections.size() == 1) {
        // if there is exactly one face projection, then it is a contact interaction
        has_contact = true;
        contact_index = face_projections[0];

        double closest_pt_value[3], proj_vector[3], distance;
        //closest_points.getPointValue(face_projections[0], closest_pt_value);
#ifdef NIMBLE_NVIDIA_BUILD
        closest_pt_value[0] = closest_points.m_hostData(face_projections[0], 0);
        closest_pt_value[1] = closest_points.m_hostData(face_projections[0], 1);
        closest_pt_value[2] = closest_points.m_hostData(face_projections[0], 2);
#else
        closest_pt_value[0] = closest_points.m_hostData(face_projections[0]).X();
        closest_pt_value[1] = closest_points.m_hostData(face_projections[0]).Y();
        closest_pt_value[2] = closest_points.m_hostData(face_projections[0]).Z();
#endif
        for (int i=0 ; i<3 ; i++) {
          proj_vector[i] = closest_pt_value[i] - pt_value[i];
        }
        distance = std::sqrt(proj_vector[0]*proj_vector[0] + proj_vector[1]*proj_vector[1] + proj_vector[2]*proj_vector[2]);
      }
      else if (face_projections.size() > 1) {
        // if there are multiple face projections, use the one with the minimum projection distance
        has_contact = true;
        double min_distance(std::numeric_limits<double>::max());
        double closest_pt_value[3], proj_vector[3], distance;
        for (auto & collision_index : face_projections) {
          //closest_points.getPointValue(collision_index, closest_pt_value);
#ifdef NIMBLE_NVIDIA_BUILD
          closest_pt_value[0] = closest_points.m_hostData(collision_index, 0);
          closest_pt_value[1] = closest_points.m_hostData(collision_index, 1);
          closest_pt_value[2] = closest_points.m_hostData(collision_index, 2);
#else
          closest_pt_value[0] = closest_points.m_hostData(collision_index, 0);
          closest_pt_value[1] = closest_points.m_hostData(collision_index, 1);
          closest_pt_value[2] = closest_points.m_hostData(collision_index, 2);
#endif
          for (int i=0 ; i<3 ; i++) {
            proj_vector[i] = closest_pt_value[i] - pt_value[i];
          }
          distance = std::sqrt(proj_vector[0]*proj_vector[0] + proj_vector[1]*proj_vector[1] + proj_vector[2]*proj_vector[2]);
          if (distance < min_distance) {
            contact_index = collision_index;
            min_distance = distance;
          }
        }
      }

      if (has_contact) {
        if (contact_index == -1) {
          throw std::logic_error("\nError in ComputeContactForce(), contact_index == -1.\n");
        }
        contact_collisions.push_back(contact_index);
      }
    }

    // compute penalty contact force
    for (unsigned int i=0 ; i<force_.size() ; ++i) {
      force_[i] = 0.0;
    }
    for (auto const & contact_index : contact_collisions) {
        int contact_node_index = collision_list.hm_data(contact_index,0);
        int contact_face_index = collision_list.hm_data(contact_index,1);
        ContactEntity & node = contact_nodes_h_[contact_node_index];
        ContactEntity & face = contact_faces_h_[contact_face_index];

        double pt_value[3], closest_pt_value[3];
#ifdef NIMBLE_NVIDIA_BUILD
        //points.getPointValue(contact_index, pt_value);
        pt_value[0] = points.m_hostData(contact_index, 0);
        pt_value[1] = points.m_hostData(contact_index, 1);
        pt_value[2] = points.m_hostData(contact_index, 2);
        //closest_points.getPointValue(contact_index, closest_pt_value);
        closest_pt_value[0] = closest_points.m_hostData(contact_index, 0);
        closest_pt_value[1] = closest_points.m_hostData(contact_index, 1);
        closest_pt_value[2] = closest_points.m_hostData(contact_index, 2);
#else
        //points.getPointValue(contact_index, pt_value);
        pt_value[0] = points.m_hostData(contact_index).X();
        pt_value[1] = points.m_hostData(contact_index).Y();
        pt_value[2] = points.m_hostData(contact_index).Z();
        //closest_points.getPointValue(contact_index, closest_pt_value);
        closest_pt_value[0] = closest_points.m_hostData(contact_index).X();
        closest_pt_value[1] = closest_points.m_hostData(contact_index).Y();
        closest_pt_value[2] = closest_points.m_hostData(contact_index).Z();
#endif

        double tri_node_1[3], tri_node_2[3], tri_node_3[3];
#ifdef NIMBLE_NVIDIA_BUILD
        //triangles.getVertexValue(contact_index, 0, tri_node_1);
        tri_node_1[0] = triangles.m_hostData(contact_index, 0, 0);
        tri_node_1[1] = triangles.m_hostData(contact_index, 0, 1);
        tri_node_1[2] = triangles.m_hostData(contact_index, 0, 2);
        //triangles.getVertexValue(contact_index, 1, tri_node_2);
        tri_node_2[0] = triangles.m_hostData(contact_index, 1, 0);
        tri_node_2[1] = triangles.m_hostData(contact_index, 1, 1);
        tri_node_2[2] = triangles.m_hostData(contact_index, 1, 2);
        //triangles.getVertexValue(contact_index, 2, tri_node_3);
        tri_node_3[0] = triangles.m_hostData(contact_index, 2, 0);
        tri_node_3[1] = triangles.m_hostData(contact_index, 2, 1);
        tri_node_3[2] = triangles.m_hostData(contact_index, 2, 2);
#else
        //triangles.getVertexValue(contact_index, 0, tri_node_1);
        tri_node_1[0] = triangles.m_hostData(contact_index).GetNode(0)[0];
        tri_node_1[1] = triangles.m_hostData(contact_index).GetNode(0)[1];
        tri_node_1[2] = triangles.m_hostData(contact_index).GetNode(0)[2];
        //triangles.getVertexValue(contact_index, 1, tri_node_2);
        tri_node_2[0] = triangles.m_hostData(contact_index).GetNode(1)[0];
        tri_node_2[1] = triangles.m_hostData(contact_index).GetNode(1)[1];
        tri_node_2[2] = triangles.m_hostData(contact_index).GetNode(1)[2];
        //triangles.getVertexValue(contact_index, 2, tri_node_3);
        tri_node_3[0] = triangles.m_hostData(contact_index).GetNode(2)[0];
        tri_node_3[1] = triangles.m_hostData(contact_index).GetNode(2)[1];
        tri_node_3[2] = triangles.m_hostData(contact_index).GetNode(2)[2];
#endif

        double tri_edge_1[3], tri_edge_2[3], tri_normal[3];
        for (int i=0 ; i<3 ; i++) {
          tri_edge_1[i] = tri_node_2[i] - tri_node_1[i];
          tri_edge_2[i] = tri_node_3[i] - tri_node_2[i];
        }
        tri_normal[0] = tri_edge_1[1]*tri_edge_2[2] - tri_edge_1[2]*tri_edge_2[1];
        tri_normal[1] = tri_edge_1[0]*tri_edge_2[2] - tri_edge_1[2]*tri_edge_2[0];
        tri_normal[2] = tri_edge_1[0]*tri_edge_2[1] - tri_edge_1[1]*tri_edge_2[0];

        double gap =
          (pt_value[0] - closest_pt_value[0])*tri_normal[0] +
          (pt_value[1] - closest_pt_value[1])*tri_normal[1] +
          (pt_value[2] - closest_pt_value[2])*tri_normal[2] ;

        double tri_normal_magnitude = std::sqrt(tri_normal[0]*tri_normal[0] + tri_normal[1]*tri_normal[1] + tri_normal[2]*tri_normal[2]);

        double contact_force[3];
        for (int i=0 ; i<3 ; ++i) {
          contact_force[i] = -1.0 * penalty_parameter_ * gap * tri_normal[i] / tri_normal_magnitude;
        }

        node.ComputeNodalContactForces(contact_force,
                                       closest_pt_value);

        for (int i=0 ; i<3 ; ++i) {
          contact_force[i] *= -1.0;
        }
        face.ComputeNodalContactForces(contact_force,
                                       closest_pt_value);

        node.GetForces(force_.data());
        face.GetForces(force_.data());
    }

    nimble_kokkos::HostScalarNodeView force_h("contact force_h", force_d_.extent(0));
    for (int i=0 ; i<force_.size() ; i++) {
      force_h(i) = force_[i];
    }
    Kokkos::deep_copy(force_d_, force_h);
#endif
  }

}
