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
      device_node_data_.emplace_back( new Field< FieldType::DeviceScalar >( label, num_objects ) );
    }
    else if (length == nimble::VECTOR) {
      device_node_data_.emplace_back( new Field< FieldType::DeviceVector >( label, num_objects ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid device data length in nimble_kokkos::ModelData::AllocateNodeData().\n");
    }

    field_id_to_device_node_data_index_[field_id] = device_node_data_.size() - 1;

    FieldBase * d_field = device_node_data_.back().get();

    if (d_field->type() == FieldType::DeviceScalar) {
      Field< FieldType::DeviceScalar > * field = dynamic_cast< Field< FieldType::DeviceScalar> * >( d_field );
      Field< FieldType::DeviceScalar >::View d_view = field->data();
      auto h_view = Kokkos::create_mirror_view( d_view );
      host_node_data_.emplace_back( new Field< FieldType::HostScalar >( h_view ) );
    }
    else if (d_field->type() == FieldType::DeviceVector) {
      Field< FieldType::DeviceVector > * field = dynamic_cast< Field< FieldType::DeviceVector> * >( d_field );
      Field< FieldType::DeviceVector >::View d_view = field->data();
      auto h_view = Kokkos::create_mirror_view( d_view );
      host_node_data_.emplace_back( new Field< FieldType::HostVector >( h_view ) );
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
    int block_index = block_id_to_element_data_index_[block_id];

    if (length == nimble::SCALAR) {
      // m_host_element_block_data is of type std::vector< std::vector< std::unique_ptr< FieldBase > > >
      device_element_data_[block_index].emplace_back( new Field< FieldType::DeviceScalar >( label, num_objects ) );
    }
    else if (length == nimble::VECTOR) {
      device_element_data_[block_index].emplace_back( new Field< FieldType::DeviceVector >( label, num_objects ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid device data length in nimble_kokkos::ModelData::AllocateElementData().\n");
    }

    field_id_to_device_element_data_index_[block_index][field_id] = device_element_data_[block_index].size() - 1;

    FieldBase * d_field = device_element_data_[block_index].back().get();

    if (d_field->type() == FieldType::DeviceScalar) {
      Field< FieldType::DeviceScalar > * field = dynamic_cast< Field< FieldType::DeviceScalar> * >( d_field );
      Field< FieldType::DeviceScalar >::View d_view = field->data();
      auto h_view = Kokkos::create_mirror_view( d_view );
      host_element_data_[block_index].emplace_back( new Field< FieldType::HostScalar >( h_view ) );
    }
    else if (d_field->type() == FieldType::DeviceVector) {
      Field< FieldType::DeviceVector > * field = dynamic_cast< Field< FieldType::DeviceVector> * >( d_field );
      Field< FieldType::DeviceVector >::View d_view = field->data();
      auto h_view = Kokkos::create_mirror_view( d_view );
      host_element_data_[block_index].emplace_back( new Field< FieldType::HostVector >( h_view ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid host data length in nimble_kokkos::ModelData::AllocateElementData().\n");
    }

    field_id_to_host_element_data_index_[block_index][field_id] = host_element_data_.size() - 1;

    return field_id;
  }

  int ModelData::AllocateIntegrationPointData(int block_id,
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

    if (block_id_to_integration_point_data_index_.find(block_id) == block_id_to_integration_point_data_index_.end()) {
      block_id_to_integration_point_data_index_[block_id] = host_integration_point_data_step_n_.size();
      host_integration_point_data_step_n_.push_back(std::vector< Data >());
      host_integration_point_data_step_np1_.push_back(std::vector< Data >());
      device_integration_point_data_step_n_.push_back(std::vector< Data >());
      device_integration_point_data_step_np1_.push_back(std::vector< Data >());
      field_id_to_host_integration_point_data_index_.push_back(std::map<int, int>());
      field_id_to_device_integration_point_data_index_.push_back(std::map<int, int>());
    }
    int block_index = block_id_to_integration_point_data_index_[block_id];

    if (length == nimble::SCALAR) {
      device_integration_point_data_step_n_[block_index].emplace_back( new Field< FieldType::DeviceScalar >( label, num_objects ) );
      device_integration_point_data_step_np1_[block_index].emplace_back( new Field< FieldType::DeviceScalar >( label, num_objects ) );
    }
    else if (length == nimble::VECTOR) {
      device_integration_point_data_step_n_[block_index].emplace_back( new Field< FieldType::DeviceVector >( label, num_objects ) );
      device_integration_point_data_step_np1_[block_index].emplace_back( new Field< FieldType::DeviceVector >( label, num_objects ) );
    }
    else if (length == nimble::SYMMETRIC_TENSOR) {
      device_integration_point_data_step_n_[block_index].emplace_back( new Field< FieldType::DeviceSymTensor >( label, num_objects ) );
      device_integration_point_data_step_np1_[block_index].emplace_back( new Field< FieldType::DeviceSymTensor >( label, num_objects ) );
    }
    else if (length == nimble::FULL_TENSOR) {
      device_integration_point_data_step_n_[block_index].emplace_back( new Field< FieldType::DeviceFullTensor >( label, num_objects ) );
      device_integration_point_data_step_np1_[block_index].emplace_back( new Field< FieldType::DeviceFullTensor >( label, num_objects ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid device data length in nimble_kokkos::ModelData::AllocateElementData().\n");
    }

    field_id_to_device_integration_point_data_index_[block_index][field_id] = device_integration_point_data_step_n_[block_index].size() - 1;

    FieldBase * d_field_step_n = device_integration_point_data_step_n_[block_index].back().get();
    FieldBase * d_field_step_np1 = device_integration_point_data_step_np1_[block_index].back().get();

    if (d_field_step_n->type() == FieldType::DeviceScalar) {
      Field< FieldType::DeviceScalar > * field_step_n = dynamic_cast< Field< FieldType::DeviceScalar> * >( d_field_step_n );
      Field< FieldType::DeviceScalar >::View d_view_step_n = field_step_n->data();
      auto h_view_step_n = Kokkos::create_mirror_view( d_view_step_n );
      host_integration_point_data_step_n_[block_index].emplace_back( new Field< FieldType::HostScalar >( h_view_step_n ) );

      Field< FieldType::DeviceScalar > * field_step_np1 = dynamic_cast< Field< FieldType::DeviceScalar> * >( d_field_step_np1 );
      Field< FieldType::DeviceScalar >::View d_view_step_np1 = field_step_np1->data();
      auto h_view_step_np1 = Kokkos::create_mirror_view( d_view_step_np1 );
      host_integration_point_data_step_np1_[block_index].emplace_back( new Field< FieldType::HostScalar >( h_view_step_np1 ) );
    }
    else if (d_field_step_n->type() == FieldType::DeviceVector) {
      Field< FieldType::DeviceVector > * field_step_n = dynamic_cast< Field< FieldType::DeviceVector> * >( d_field_step_n );
      Field< FieldType::DeviceVector >::View d_view_step_n = field_step_n->data();
      auto h_view_step_n = Kokkos::create_mirror_view( d_view_step_n );
      host_integration_point_data_step_n_[block_index].emplace_back( new Field< FieldType::HostVector >( h_view_step_n ) );

      Field< FieldType::DeviceVector > * field_step_np1 = dynamic_cast< Field< FieldType::DeviceVector> * >( d_field_step_np1 );
      Field< FieldType::DeviceVector >::View d_view_step_np1 = field_step_np1->data();
      auto h_view_step_np1 = Kokkos::create_mirror_view( d_view_step_np1 );
      host_integration_point_data_step_np1_[block_index].emplace_back( new Field< FieldType::HostVector >( h_view_step_np1 ) );
    }
    else if (d_field_step_n->type() == FieldType::DeviceSymTensor) {
      Field< FieldType::DeviceSymTensor > * field_step_n = dynamic_cast< Field< FieldType::DeviceSymTensor> * >( d_field_step_n );
      Field< FieldType::DeviceSymTensor >::View d_view_step_n = field_step_n->data();
      auto h_view_step_n = Kokkos::create_mirror_view( d_view_step_n );
      host_integration_point_data_step_n_[block_index].emplace_back( new Field< FieldType::HostSymTensor >( h_view_step_n ) );

      Field< FieldType::DeviceSymTensor > * field_step_np1 = dynamic_cast< Field< FieldType::DeviceSymTensor> * >( d_field_step_np1 );
      Field< FieldType::DeviceSymTensor >::View d_view_step_np1 = field_step_np1->data();
      auto h_view_step_np1 = Kokkos::create_mirror_view( d_view_step_np1 );
      host_integration_point_data_step_np1_[block_index].emplace_back( new Field< FieldType::HostSymTensor >( h_view_step_np1 ) );
    }
    else if (d_field_step_n->type() == FieldType::DeviceFullTensor) {
      Field< FieldType::DeviceFullTensor > * field_step_n = dynamic_cast< Field< FieldType::DeviceFullTensor> * >( d_field_step_n );
      Field< FieldType::DeviceFullTensor >::View d_view_step_n = field_step_n->data();
      auto h_view_step_n = Kokkos::create_mirror_view( d_view_step_n );
      host_integration_point_data_step_n_[block_index].emplace_back( new Field< FieldType::HostFullTensor >( h_view_step_n ) );

      Field< FieldType::DeviceFullTensor > * field_step_np1 = dynamic_cast< Field< FieldType::DeviceFullTensor> * >( d_field_step_np1 );
      Field< FieldType::DeviceFullTensor >::View d_view_step_np1 = field_step_np1->data();
      auto h_view_step_np1 = Kokkos::create_mirror_view( d_view_step_np1 );
      host_integration_point_data_step_np1_[block_index].emplace_back( new Field< FieldType::HostFullTensor >( h_view_step_np1 ) );
    }
    else {
      throw std::logic_error("\nError:  Invalid host data length in nimble_kokkos::ModelData::AllocateElementData().\n");
    }

    field_id_to_host_integration_point_data_index_[block_index][field_id] = host_integration_point_data_step_n_.size() - 1;

    return field_id;
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
          if (host_node_data_.at(node_data_index)->type() == FieldType::HostScalar)
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
          if (host_node_data_.at(node_data_index)->type() == FieldType::HostVector)
            node_data_labels.push_back(field_label);
        }
      }
    }
    return node_data_labels;
  }

  HostScalarView ModelData::GetHostScalarNodeData(int field_id) {
    int index = field_id_to_host_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = host_node_data_.at(index).get();
    Field< FieldType::HostScalar >* derived_field_ptr = static_cast< Field< FieldType::HostScalar >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  HostVectorView ModelData::GetHostVectorNodeData(int field_id) {
    int index = field_id_to_host_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = host_node_data_.at(index).get();
    Field< FieldType::HostVector >* derived_field_ptr = static_cast< Field< FieldType::HostVector >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceScalarView ModelData::GetDeviceScalarNodeData(int field_id) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceScalar >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalar >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceVectorView ModelData::GetDeviceVectorNodeData(int field_id) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceVector >* derived_field_ptr = static_cast< Field< FieldType::DeviceVector >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceFullTensorView ModelData::GetDeviceFullTensorIntegrationPointData(int block_id,
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
    Field< FieldType::DeviceFullTensor >* derived_field_ptr = static_cast< Field< FieldType::DeviceFullTensor >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceSymTensorView ModelData::GetDeviceSymTensorIntegrationPointData(int block_id,
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
    Field< FieldType::DeviceSymTensor >* derived_field_ptr = static_cast< Field< FieldType::DeviceSymTensor >* >(base_field_ptr);
    return derived_field_ptr->data();
  }

  DeviceScalarGatheredView ModelData::GatherScalarNodeData(int field_id,
                                                           int num_elements,
                                                           int num_nodes_per_element,
                                                           DeviceElementConnectivityView elem_conn_d,
                                                           DeviceScalarGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceScalar >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalar >* >(base_field_ptr);
    auto data = derived_field_ptr->data();
    Kokkos::parallel_for("GatherScalarNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
        for (int i_node=0 ; i_node < num_nodes_per_element ; i_node++) {
          gathered_view_d(i_elem, i_node) = data(elem_conn_d(num_nodes_per_element*i_elem + i_node));
        }
      });
    return gathered_view_d;
  }

  DeviceVectorGatheredView ModelData::GatherVectorNodeData(int field_id,
                                                           int num_elements,
                                                           int num_nodes_per_element,
                                                           DeviceElementConnectivityView elem_conn_d,
                                                           DeviceVectorGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceVector >* derived_field_ptr = static_cast< Field< FieldType::DeviceVector >* >(base_field_ptr);
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
                                        DeviceScalarGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceScalar >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalar >* >(base_field_ptr);
    Field< FieldType::DeviceScalar >::AtomicView data = derived_field_ptr->data();
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
                                        DeviceVectorGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceVector >* derived_field_ptr = static_cast< Field< FieldType::DeviceVector >* >(base_field_ptr);
    Field< FieldType::DeviceVector >::AtomicView data = derived_field_ptr->data();
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
                                                              DeviceScalarGatheredView gathered_view_d) {
    int index = field_id_to_device_node_data_index_.at(field_id);
    FieldBase* base_field_ptr = device_node_data_.at(index).get();
    Field< FieldType::DeviceScalar >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalar >* >(base_field_ptr);
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
