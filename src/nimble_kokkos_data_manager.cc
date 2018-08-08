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

#include "nimble_kokkos_data_manager.h"
#include <stdexcept>

#include <Kokkos_ScatterView.hpp>

namespace nimble_kokkos {

  int ModelData::AllocateNodeData(nimble::Length length,
                                  std::string label,
                                  int num_objects) {

    int field_id;
    auto it = field_label_to_field_id_map_.find(label);
    if (it == field_label_to_field_id_map_.end()) {
      field_id = field_label_to_field_id_map_.size();
      field_label_to_field_id_map_[label] = field_id;
    }
    else {
      field_id = it->second;
    }

    if (length == nimble::SCALAR) {
      // device_node_data_ is of type std::vector< std::unique_ptr< FieldBase > >
      device_node_data_.emplace_back( new Field< FieldType::DeviceScalarNode >( label, num_objects ) );
    }
    else if (length == nimble::VECTOR) {
      device_node_data_.emplace_back( new Field< FieldType::DeviceVectorNode >( label, num_objects ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid device data length in nimble_kokkos::ModelData::AllocateNodeData().\n");
    }

    field_id_to_device_node_data_index_[field_id] = device_node_data_.size() - 1;

    FieldBase * d_field = device_node_data_.back().get();

    if (d_field->type() == FieldType::DeviceScalarNode) {
      Field< FieldType::DeviceScalarNode > * field = dynamic_cast< Field< FieldType::DeviceScalarNode> * >( d_field );
      Field< FieldType::DeviceScalarNode >::View d_view = field->data();
      auto h_view = Kokkos::create_mirror_view( d_view );
      host_node_data_.emplace_back( new Field< FieldType::HostScalarNode >( h_view ) );
    }
    else if (d_field->type() == FieldType::DeviceVectorNode) {
      Field< FieldType::DeviceVectorNode > * field = dynamic_cast< Field< FieldType::DeviceVectorNode> * >( d_field );
      Field< FieldType::DeviceVectorNode >::View d_view = field->data();
      auto h_view = Kokkos::create_mirror_view( d_view );
      host_node_data_.emplace_back( new Field< FieldType::HostVectorNode >( h_view ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid host data length in nimble_kokkos::ModelData::AllocateNodeData().\n");
    }

    field_id_to_host_node_data_index_[field_id] = host_node_data_.size() - 1;

    return field_id;
  }

  int ModelData::AllocateElementData(int block_id,
                                     nimble::Length length,
                                     std::string label,
                                     int num_objects) {
    int field_id;
    auto it = field_label_to_field_id_map_.find(label);
    if (it == field_label_to_field_id_map_.end()) {
      field_id = field_label_to_field_id_map_.size();
      field_label_to_field_id_map_[label] = field_id;
    }
    else {
      field_id = it->second;
    }

    if (block_id_to_element_data_index_.find(block_id) == block_id_to_element_data_index_.end()) {
      block_id_to_element_data_index_[block_id] = host_element_data_.size();
      host_element_data_.push_back(std::vector< Data >());
      device_element_data_.push_back(std::vector< Data >());
      field_id_to_host_element_data_index_.push_back(std::map<int, int>());
      field_id_to_device_element_data_index_.push_back(std::map<int, int>());
    }
    int block_index = block_id_to_element_data_index_.at(block_id);

    if (length == nimble::SCALAR) {
      device_element_data_.at(block_index).emplace_back( new Field< FieldType::DeviceScalarElem >( label, num_objects ) );
    }
    else if (length == nimble::SYMMETRIC_TENSOR) {
      device_element_data_.at(block_index).emplace_back( new Field< FieldType::DeviceSymTensorElem >( label, num_objects ) );
    }
    else if (length == nimble::FULL_TENSOR) {
      device_element_data_.at(block_index).emplace_back( new Field< FieldType::DeviceFullTensorElem >( label, num_objects ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid device data length in nimble_kokkos::ModelData::AllocateElementData().\n");
    }

    field_id_to_device_element_data_index_.at(block_index)[field_id] = device_element_data_.at(block_index).size() - 1;

    FieldBase * d_field = device_element_data_.at(block_index).back().get();

    if (d_field->type() == FieldType::DeviceScalarElem) {
      Field< FieldType::DeviceScalarElem > * field = dynamic_cast< Field< FieldType::DeviceScalarElem> * >( d_field );
      Field< FieldType::DeviceScalarElem >::View d_view = field->data();
      auto h_view = Kokkos::create_mirror_view( d_view );
      host_element_data_.at(block_index).emplace_back( new Field< FieldType::HostScalarElem >( h_view ) );
    }
    else if (d_field->type() == FieldType::DeviceSymTensorElem) {
      Field< FieldType::DeviceSymTensorElem > * field = dynamic_cast< Field< FieldType::DeviceSymTensorElem> * >( d_field );
      Field< FieldType::DeviceSymTensorElem >::View d_view = field->data();
      auto h_view = Kokkos::create_mirror_view( d_view );
      host_element_data_.at(block_index).emplace_back( new Field< FieldType::HostSymTensorElem >( h_view ) );
    }
    else if (d_field->type() == FieldType::DeviceFullTensorElem) {
      Field< FieldType::DeviceFullTensorElem > * field = dynamic_cast< Field< FieldType::DeviceFullTensorElem> * >( d_field );
      Field< FieldType::DeviceFullTensorElem >::View d_view = field->data();
      auto h_view = Kokkos::create_mirror_view( d_view );
      host_element_data_.at(block_index).emplace_back( new Field< FieldType::HostFullTensorElem >( h_view ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid host data length in nimble_kokkos::ModelData::AllocateElementData().\n");
    }

    field_id_to_host_element_data_index_.at(block_index)[field_id] = host_element_data_.at(block_index).size() - 1;

    return field_id;
  }

  int ModelData::AllocateIntegrationPointData(int block_id,
                                              nimble::Length length,
                                              std::string label,
                                              int num_objects,
                                              std::vector<double> initial_value) {
    bool set_initial_value = false;
    if (initial_value.size() > 0) {
      set_initial_value = true;
    }

    int field_id;
    auto it = field_label_to_field_id_map_.find(label);
    if (it == field_label_to_field_id_map_.end()) {
      field_id = field_label_to_field_id_map_.size();
      field_label_to_field_id_map_[label] = field_id;
    }
    else {
      field_id = it->second;
    }

    if (block_id_to_integration_point_data_index_.find(block_id) == block_id_to_integration_point_data_index_.end()) {
      block_id_to_integration_point_data_index_[block_id] = host_integration_point_data_step_n_.size();
      host_integration_point_data_step_n_.push_back(std::vector< Data >());
      host_integration_point_data_step_np1_.push_back(std::vector< Data >());
      device_integration_point_data_step_n_.push_back(std::vector< Data >());
      device_integration_point_data_step_np1_.push_back(std::vector< Data >());
      field_id_to_host_integration_point_data_index_.push_back(std::map<int, int>());
      field_id_to_device_integration_point_data_index_.push_back(std::map<int, int>());
    }
    int block_index = block_id_to_integration_point_data_index_.at(block_id);

    if (length == nimble::SCALAR) {
      device_integration_point_data_step_n_.at(block_index).emplace_back( new Field< FieldType::DeviceScalarNode >( label, num_objects ) );
      device_integration_point_data_step_np1_.at(block_index).emplace_back( new Field< FieldType::DeviceScalarNode >( label, num_objects ) );
    }
    else if (length == nimble::VECTOR) {
      device_integration_point_data_step_n_.at(block_index).emplace_back( new Field< FieldType::DeviceVectorNode >( label, num_objects ) );
      device_integration_point_data_step_np1_.at(block_index).emplace_back( new Field< FieldType::DeviceVectorNode >( label, num_objects ) );
    }
    else if (length == nimble::SYMMETRIC_TENSOR) {
      device_integration_point_data_step_n_.at(block_index).emplace_back( new Field< FieldType::DeviceSymTensorIntPt >( label, num_objects ) );
      device_integration_point_data_step_np1_.at(block_index).emplace_back( new Field< FieldType::DeviceSymTensorIntPt >( label, num_objects ) );
    }
    else if (length == nimble::FULL_TENSOR) {
      device_integration_point_data_step_n_.at(block_index).emplace_back( new Field< FieldType::DeviceFullTensorIntPt >( label, num_objects ) );
      device_integration_point_data_step_np1_.at(block_index).emplace_back( new Field< FieldType::DeviceFullTensorIntPt >( label, num_objects ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid device data length in nimble_kokkos::ModelData::AllocateIntegrationPointData().\n");
    }

    field_id_to_device_integration_point_data_index_.at(block_index)[field_id] = device_integration_point_data_step_n_.at(block_index).size() - 1;

    FieldBase * d_field_step_n = device_integration_point_data_step_n_.at(block_index).back().get();
    FieldBase * d_field_step_np1 = device_integration_point_data_step_np1_.at(block_index).back().get();

    if (d_field_step_n->type() == FieldType::DeviceScalarNode) {
      Field< FieldType::DeviceScalarNode > * field_step_n = dynamic_cast< Field< FieldType::DeviceScalarNode> * >( d_field_step_n );
      Field< FieldType::DeviceScalarNode >::View d_view_step_n = field_step_n->data();
      auto h_view_step_n = Kokkos::create_mirror_view( d_view_step_n );
      host_integration_point_data_step_n_.at(block_index).emplace_back( new Field< FieldType::HostScalarNode >( h_view_step_n ) );

      Field< FieldType::DeviceScalarNode > * field_step_np1 = dynamic_cast< Field< FieldType::DeviceScalarNode> * >( d_field_step_np1 );
      Field< FieldType::DeviceScalarNode >::View d_view_step_np1 = field_step_np1->data();
      auto h_view_step_np1 = Kokkos::create_mirror_view( d_view_step_np1 );
      host_integration_point_data_step_np1_.at(block_index).emplace_back( new Field< FieldType::HostScalarNode >( h_view_step_np1 ) );

      if (set_initial_value) {
        int num_elem = h_view_step_n.extent(0);
        for (int i_elem=0 ; i_elem<num_elem ; ++i_elem) {
          h_view_step_n(i_elem) = initial_value.at(0);
          h_view_step_np1(i_elem) = initial_value.at(0);
        }
        Kokkos::deep_copy(d_view_step_n, h_view_step_n);
        Kokkos::deep_copy(d_view_step_np1, h_view_step_np1);
      }
    }
    else if (d_field_step_n->type() == FieldType::DeviceVectorNode) {
      Field< FieldType::DeviceVectorNode > * field_step_n = dynamic_cast< Field< FieldType::DeviceVectorNode> * >( d_field_step_n );
      Field< FieldType::DeviceVectorNode >::View d_view_step_n = field_step_n->data();
      auto h_view_step_n = Kokkos::create_mirror_view( d_view_step_n );
      host_integration_point_data_step_n_.at(block_index).emplace_back( new Field< FieldType::HostVectorNode >( h_view_step_n ) );

      Field< FieldType::DeviceVectorNode > * field_step_np1 = dynamic_cast< Field< FieldType::DeviceVectorNode> * >( d_field_step_np1 );
      Field< FieldType::DeviceVectorNode >::View d_view_step_np1 = field_step_np1->data();
      auto h_view_step_np1 = Kokkos::create_mirror_view( d_view_step_np1 );
      host_integration_point_data_step_np1_.at(block_index).emplace_back( new Field< FieldType::HostVectorNode >( h_view_step_np1 ) );

      if (set_initial_value) {
        int num_elem = h_view_step_n.extent(0);
        int num_entries = 3;
        for (int i_elem=0 ; i_elem<num_elem ; ++i_elem) {
          for (int i_entry=0 ; i_entry<num_entries ; ++i_entry) {
            h_view_step_n(i_elem, i_entry) = initial_value.at(i_entry);
            h_view_step_np1(i_elem, i_entry) = initial_value.at(i_entry);
          }
        }
        Kokkos::deep_copy(d_view_step_n, h_view_step_n);
        Kokkos::deep_copy(d_view_step_np1, h_view_step_np1);
      }
    }
    else if (d_field_step_n->type() == FieldType::DeviceSymTensorIntPt) {
      Field< FieldType::DeviceSymTensorIntPt > * field_step_n = dynamic_cast< Field< FieldType::DeviceSymTensorIntPt> * >( d_field_step_n );
      Field< FieldType::DeviceSymTensorIntPt >::View d_view_step_n = field_step_n->data();
      auto h_view_step_n = Kokkos::create_mirror_view( d_view_step_n );
      host_integration_point_data_step_n_.at(block_index).emplace_back( new Field< FieldType::HostSymTensorIntPt >( h_view_step_n ) );

      Field< FieldType::DeviceSymTensorIntPt > * field_step_np1 = dynamic_cast< Field< FieldType::DeviceSymTensorIntPt> * >( d_field_step_np1 );
      Field< FieldType::DeviceSymTensorIntPt >::View d_view_step_np1 = field_step_np1->data();
      auto h_view_step_np1 = Kokkos::create_mirror_view( d_view_step_np1 );
      host_integration_point_data_step_np1_.at(block_index).emplace_back( new Field< FieldType::HostSymTensorIntPt >( h_view_step_np1 ) );

      if (set_initial_value) {
        int num_elem = h_view_step_n.extent(0);
        int num_int_pt = h_view_step_n.extent(1);
        int num_entries = 6;
        for (int i_elem=0 ; i_elem<num_elem ; ++i_elem) {
          for (int i_int_pt=0 ; i_int_pt<num_int_pt ; ++i_int_pt) {
            for (int i_entry=0 ; i_entry<num_entries ; ++i_entry) {
              h_view_step_n(i_elem, i_int_pt, i_entry) = initial_value.at(i_entry);
              h_view_step_np1(i_elem, i_int_pt, i_entry) = initial_value.at(i_entry);
            }
          }
        }
        Kokkos::deep_copy(d_view_step_n, h_view_step_n);
        Kokkos::deep_copy(d_view_step_np1, h_view_step_np1);
      }
    }
    else if (d_field_step_n->type() == FieldType::DeviceFullTensorIntPt) {
      Field< FieldType::DeviceFullTensorIntPt > * field_step_n = dynamic_cast< Field< FieldType::DeviceFullTensorIntPt> * >( d_field_step_n );
      Field< FieldType::DeviceFullTensorIntPt >::View d_view_step_n = field_step_n->data();
      auto h_view_step_n = Kokkos::create_mirror_view( d_view_step_n );
      host_integration_point_data_step_n_.at(block_index).emplace_back( new Field< FieldType::HostFullTensorIntPt >( h_view_step_n ) );

      Field< FieldType::DeviceFullTensorIntPt > * field_step_np1 = dynamic_cast< Field< FieldType::DeviceFullTensorIntPt> * >( d_field_step_np1 );
      Field< FieldType::DeviceFullTensorIntPt >::View d_view_step_np1 = field_step_np1->data();
      auto h_view_step_np1 = Kokkos::create_mirror_view( d_view_step_np1 );
      host_integration_point_data_step_np1_.at(block_index).emplace_back( new Field< FieldType::HostFullTensorIntPt >( h_view_step_np1 ) );

      if (set_initial_value) {
        int num_elem = h_view_step_n.extent(0);
        int num_int_pt = h_view_step_n.extent(1);
        int num_entries = 9;
        for (int i_elem=0 ; i_elem<num_elem ; ++i_elem) {
          for (int i_int_pt=0 ; i_int_pt<num_int_pt ; ++i_int_pt) {
            for (int i_entry=0 ; i_entry<num_entries ; ++i_entry) {
              h_view_step_n(i_elem, i_int_pt, i_entry) = initial_value.at(i_entry);
              h_view_step_np1(i_elem, i_int_pt, i_entry) = initial_value.at(i_entry);
            }
          }
        }
        Kokkos::deep_copy(d_view_step_n, h_view_step_n);
        Kokkos::deep_copy(d_view_step_np1, h_view_step_np1);
      }
    }
    else {
      throw std::logic_error("\nError:  Invalid host data length in nimble_kokkos::ModelData::AllocateElementData().\n");
    }

    field_id_to_host_integration_point_data_index_.at(block_index)[field_id] = host_integration_point_data_step_n_.at(block_index).size() - 1;

    return field_id;
  }

  std::vector<int> ModelData::GetBlockIds() const {
    std::vector<int> block_ids;
    for (auto const & entry : block_id_to_integration_point_data_index_) {
      block_ids.push_back(entry.first);
    }
    return block_ids;
  }

  std::vector<std::string> ModelData::GetScalarNodeDataLabels() const {
    std::vector<std::string> node_data_labels;
    for (auto const & entry : field_label_to_field_id_map_) {
      std::string const & field_label = entry.first;
      int field_id = entry.second;
      for (auto const & node_entry : field_id_to_host_node_data_index_) {
        int node_data_field_id = node_entry.first;
        if (field_id == node_data_field_id) {
          int node_data_index = node_entry.second;
          if (host_node_data_.at(node_data_index)->type() == FieldType::HostScalarNode)
            node_data_labels.push_back(field_label);
        }
      }
    }
    return node_data_labels;
  }

  std::vector<std::string> ModelData::GetVectorNodeDataLabels() const {
    std::vector<std::string> node_data_labels;
    for (auto const & entry : field_label_to_field_id_map_) {
      std::string const & field_label = entry.first;
      int field_id = entry.second;
      for (auto const & node_entry : field_id_to_host_node_data_index_) {
        int node_data_field_id = node_entry.first;
        if (field_id == node_data_field_id) {
          int node_data_index = node_entry.second;
          if (host_node_data_.at(node_data_index)->type() == FieldType::HostVectorNode)
            node_data_labels.push_back(field_label);
        }
      }
    }
    return node_data_labels;
  }

  std::vector<std::string> ModelData::GetSymmetricTensorIntegrationPointDataLabels(int block_id) const {
    int block_index = block_id_to_integration_point_data_index_.at(block_id);
    int num_blocks = static_cast<int>(block_id_to_integration_point_data_index_.size());
    std::vector<std::string> ipt_data_labels;
    for (auto const & entry : field_label_to_field_id_map_) {
      std::string const & field_label = entry.first;
      int field_id = entry.second;
      for (auto const & ipt_entry : field_id_to_device_integration_point_data_index_.at(block_index)) {
        int ipt_data_field_id = ipt_entry.first;
        if (field_id == ipt_data_field_id) {
          int ipt_data_index = ipt_entry.second;
          if (device_integration_point_data_step_np1_.at(block_index).at(ipt_data_index)->type() == FieldType::DeviceSymTensorIntPt) {
            if (std::find(ipt_data_labels.begin(), ipt_data_labels.end(), field_label) == ipt_data_labels.end()) {
              ipt_data_labels.push_back(field_label);
            }
          }
        }
      }
    }
    return ipt_data_labels;
  }

  std::vector<std::string> ModelData::GetFullTensorIntegrationPointDataLabels(int block_id) const {
    int block_index = block_id_to_integration_point_data_index_.at(block_id);
    std::vector<std::string> ipt_data_labels;
    for (auto const & entry : field_label_to_field_id_map_) {
      std::string const & field_label = entry.first;
      int field_id = entry.second;
      for (auto const & ipt_entry : field_id_to_device_integration_point_data_index_.at(block_index)) {
        int ipt_data_field_id = ipt_entry.first;
        if (field_id == ipt_data_field_id) {
          int ipt_data_index = ipt_entry.second;
          if (device_integration_point_data_step_np1_.at(block_index).at(ipt_data_index)->type() == FieldType::DeviceFullTensorIntPt) {
            if (std::find(ipt_data_labels.begin(), ipt_data_labels.end(), field_label) == ipt_data_labels.end()) {
              ipt_data_labels.push_back(field_label);
            }
          }
        }
      }
    }
    return ipt_data_labels;
  }

  HostScalarNodeView ModelData::GetHostScalarNodeData(int field_id) {
    int index = field_id_to_host_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = host_node_data_.at(index).get();
    Field< FieldType::HostScalarNode >* derived_field_ptr = static_cast< Field< FieldType::HostScalarNode >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  HostVectorNodeView ModelData::GetHostVectorNodeData(int field_id) {
    int index = field_id_to_host_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = host_node_data_.at(index).get();
    Field< FieldType::HostVectorNode >* derived_field_ptr = static_cast< Field< FieldType::HostVectorNode >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  HostScalarElemView ModelData::GetHostScalarElementData(int block_id,
                                                         int field_id) {
    int block_index = block_id_to_element_data_index_.at(block_id);
    int data_index = field_id_to_host_element_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    base_field_ptr = host_element_data_.at(block_index).at(data_index).get();
    Field< FieldType::HostScalarElem >* derived_field_ptr = static_cast< Field< FieldType::HostScalarElem >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  HostSymTensorIntPtView ModelData::GetHostSymTensorIntegrationPointData(int block_id,
                                                                         int field_id,
                                                                         nimble::Step step) {
    int block_index = block_id_to_integration_point_data_index_.at(block_id);
    int data_index = field_id_to_host_integration_point_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    if (step == nimble::STEP_N) {
      base_field_ptr = host_integration_point_data_step_n_.at(block_index).at(data_index).get();
    }
    else if (step == nimble::STEP_NP1) {
      base_field_ptr = host_integration_point_data_step_np1_.at(block_index).at(data_index).get();
    }
    Field< FieldType::HostSymTensorIntPt >* derived_field_ptr = static_cast< Field< FieldType::HostSymTensorIntPt >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  HostFullTensorIntPtView ModelData::GetHostFullTensorIntegrationPointData(int block_id,
                                                                           int field_id,
                                                                           nimble::Step step) {
    int block_index = block_id_to_integration_point_data_index_.at(block_id);
    int data_index = field_id_to_host_integration_point_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    if (step == nimble::STEP_N) {
      base_field_ptr = host_integration_point_data_step_n_.at(block_index).at(data_index).get();
    }
    else if (step == nimble::STEP_NP1) {
      base_field_ptr = host_integration_point_data_step_np1_.at(block_index).at(data_index).get();
    }
    Field< FieldType::HostFullTensorIntPt >* derived_field_ptr = static_cast< Field< FieldType::HostFullTensorIntPt >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  HostSymTensorElemView ModelData::GetHostSymTensorElementData(int block_id,
                                                               int field_id) {
    int block_index = block_id_to_element_data_index_.at(block_id);
    int data_index = field_id_to_host_element_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    base_field_ptr = host_element_data_.at(block_index).at(data_index).get();
    Field< FieldType::HostSymTensorElem >* derived_field_ptr = static_cast< Field< FieldType::HostSymTensorElem >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  HostFullTensorElemView ModelData::GetHostFullTensorElementData(int block_id,
                                                                 int field_id) {
    int block_index = block_id_to_element_data_index_.at(block_id);
    int data_index = field_id_to_host_element_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    base_field_ptr = host_element_data_.at(block_index).at(data_index).get();
    Field< FieldType::HostFullTensorElem >* derived_field_ptr = static_cast< Field< FieldType::HostFullTensorElem >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceScalarNodeView ModelData::GetDeviceScalarNodeData(int field_id) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceScalarNode >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalarNode >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceVectorNodeView ModelData::GetDeviceVectorNodeData(int field_id) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceVectorNode >* derived_field_ptr = static_cast< Field< FieldType::DeviceVectorNode >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceSymTensorIntPtView ModelData::GetDeviceSymTensorIntegrationPointData(int block_id,
                                                                             int field_id,
                                                                             nimble::Step step) {
    int block_index = block_id_to_integration_point_data_index_.at(block_id);
    int data_index = field_id_to_device_integration_point_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    if (step == nimble::STEP_N) {
      base_field_ptr = device_integration_point_data_step_n_.at(block_index).at(data_index).get();
    }
    else if (step == nimble::STEP_NP1) {
      base_field_ptr = device_integration_point_data_step_np1_.at(block_index).at(data_index).get();
    }
    Field< FieldType::DeviceSymTensorIntPt >* derived_field_ptr = static_cast< Field< FieldType::DeviceSymTensorIntPt >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceFullTensorIntPtView ModelData::GetDeviceFullTensorIntegrationPointData(int block_id,
                                                                               int field_id,
                                                                               nimble::Step step) {
    int block_index = block_id_to_integration_point_data_index_.at(block_id);
    int data_index = field_id_to_device_integration_point_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    if (step == nimble::STEP_N) {
      base_field_ptr = device_integration_point_data_step_n_.at(block_index).at(data_index).get();
    }
    else if (step == nimble::STEP_NP1) {
      base_field_ptr = device_integration_point_data_step_np1_.at(block_index).at(data_index).get();
    }
    Field< FieldType::DeviceFullTensorIntPt >* derived_field_ptr = static_cast< Field< FieldType::DeviceFullTensorIntPt >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceScalarElemView ModelData::GetDeviceScalarElementData(int block_id,
                                                             int field_id) {
    int block_index = block_id_to_element_data_index_.at(block_id);
    int data_index = field_id_to_device_element_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    base_field_ptr = device_element_data_.at(block_index).at(data_index).get();
    Field< FieldType::DeviceScalarElem >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalarElem >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceSymTensorElemView ModelData::GetDeviceSymTensorElementData(int block_id,
                                                                   int field_id) {
    int block_index = block_id_to_element_data_index_.at(block_id);
    int data_index = field_id_to_device_element_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    base_field_ptr = device_element_data_.at(block_index).at(data_index).get();
    Field< FieldType::DeviceSymTensorElem >* derived_field_ptr = static_cast< Field< FieldType::DeviceSymTensorElem >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceFullTensorElemView ModelData::GetDeviceFullTensorElementData(int block_id,
                                                                     int field_id) {
    int block_index = block_id_to_element_data_index_.at(block_id);
    int data_index = field_id_to_device_element_data_index_.at(block_index).at(field_id);
    FieldBase* base_field_ptr(0);
    base_field_ptr = device_element_data_.at(block_index).at(data_index).get();
    Field< FieldType::DeviceFullTensorElem >* derived_field_ptr = static_cast< Field< FieldType::DeviceFullTensorElem >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceScalarNodeGatheredView ModelData::GatherScalarNodeData(int field_id,
                                                               int num_elements,
                                                               int num_nodes_per_element,
                                                               DeviceElementConnectivityView elem_conn_d,
                                                               DeviceScalarNodeGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceScalarNode >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalarNode >* >(base_field_ptr);
    auto data = derived_field_ptr->data();
    Kokkos::parallel_for("GatherScalarNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
        for (int i_node=0 ; i_node < num_nodes_per_element ; i_node++) {
          gathered_view_d(i_elem, i_node) = data(elem_conn_d(num_nodes_per_element*i_elem + i_node));
        }
      });
    return gathered_view_d;
  }

  DeviceVectorNodeGatheredView ModelData::GatherVectorNodeData(int field_id,
                                                               int num_elements,
                                                               int num_nodes_per_element,
                                                               DeviceElementConnectivityView elem_conn_d,
                                                               DeviceVectorNodeGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceVectorNode >* derived_field_ptr = static_cast< Field< FieldType::DeviceVectorNode >* >(base_field_ptr);
    auto data = derived_field_ptr->data();
    Kokkos::parallel_for("GatherVectorNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
        for (int i_node=0 ; i_node < num_nodes_per_element ; i_node++) {
          int node_index = elem_conn_d(num_nodes_per_element*i_elem + i_node);
          for (int i_coord=0 ; i_coord < 3 ; i_coord++) {
            gathered_view_d(i_elem, i_node, i_coord) = data(node_index, i_coord);
          }
        }
      });
    return gathered_view_d;
  }

  void ModelData::ScatterScalarNodeData(int field_id,
                                        int num_elements,
                                        int num_nodes_per_element,
                                        DeviceElementConnectivityView elem_conn_d,
                                        DeviceScalarNodeGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceScalarNode >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalarNode >* >(base_field_ptr);
    Field< FieldType::DeviceScalarNode >::AtomicView data = derived_field_ptr->data();
    Kokkos::deep_copy(data, (double)(0.0));
    Kokkos::parallel_for("ScatterScalarNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
        for (int i_node=0 ; i_node < num_nodes_per_element ; i_node++) {
          data(elem_conn_d(num_nodes_per_element*i_elem + i_node)) += gathered_view_d(i_elem, i_node);
        }
      });
    return;
  }

  void ModelData::ScatterVectorNodeData(int field_id,
                                        int num_elements,
                                        int num_nodes_per_element,
                                        DeviceElementConnectivityView elem_conn_d,
                                        DeviceVectorNodeGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceVectorNode >* derived_field_ptr = static_cast< Field< FieldType::DeviceVectorNode >* >(base_field_ptr);
    Field< FieldType::DeviceVectorNode >::AtomicView data = derived_field_ptr->data();
    Kokkos::deep_copy(data, (double)(0.0));
    Kokkos::parallel_for("ScatterVectorNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
        for (int i_node=0 ; i_node < num_nodes_per_element ; i_node++) {
          int node_index = elem_conn_d(num_nodes_per_element*i_elem + i_node);
          for (int i_coord=0 ; i_coord < 3 ; i_coord++) {
            data(node_index, i_coord) += gathered_view_d(i_elem, i_node, i_coord);
          }
        }
      });
    return;
  }

  void ModelData::ScatterScalarNodeDataUsingKokkosScatterView(int field_id,
                                                              int num_elements,
                                                              int num_nodes_per_element,
                                                              DeviceElementConnectivityView elem_conn_d,
                                                              DeviceScalarNodeGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceScalarNode >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalarNode >* >(base_field_ptr);
    auto data = derived_field_ptr->data();
    Kokkos::deep_copy(data, (double)(0.0));
    auto scatter_view = Kokkos::Experimental::create_scatter_view(data); // DJL it is a terrible idea to allocate this here
    scatter_view.reset();
    Kokkos::parallel_for("GatherVectorNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
        auto scattered_access = scatter_view.access();
        for (int i_node=0 ; i_node < num_nodes_per_element ; i_node++) {
          scattered_access(elem_conn_d(num_nodes_per_element*i_elem + i_node)) += gathered_view_d(i_elem, i_node);
        }
      });
    Kokkos::Experimental::contribute(data, scatter_view);
    return;
  }

  ModelData& DataManager::GetMacroScaleData() {
    return macroscale_data_;
  }

} // namespace nimble
