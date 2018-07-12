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

#include "nimble_exodus_output_manager.h"
#include "nimble_utils.h"
#include <sstream>

#ifndef NIMBLE_HAVE_DARMA
  #include <vector>
#endif

namespace nimble_kokkos {

  void ExodusOutputManager::SpecifyOutputFields(nimble_kokkos::ModelData& model_data,
                                                std::string const & output_command_string) {

    // Parse the string into individual field labels
    std::vector<std::string> requested_labels;
    std::stringstream ss(output_command_string);
    std::string entry;
    while (std::getline(ss, entry, ' ')) {
      requested_labels.push_back(entry);
    }

    std::vector<std::string> scalar_node_data_labels = model_data.GetScalarNodeDataLabels();
    for (auto const & requested_label : requested_labels) {
      for (auto& node_label : scalar_node_data_labels) {
        if (requested_label == node_label) {
          int field_id = model_data.GetFieldId(node_label);
          int num_nodes = model_data.GetHostScalarNodeData(field_id).extent(0);
          node_data_labels_.push_back(node_label);
          node_data_field_ids_.push_back(field_id);
          node_data_types_.push_back(FieldType::HostScalar);
          node_data_components_.push_back(0);
          node_data_.push_back(std::vector<double>(num_nodes, 0.0));
        }
      }
    }

    std::vector<std::string> vector_node_data_labels = model_data.GetVectorNodeDataLabels();
    for (auto const & requested_label : requested_labels) {
      for (auto& node_label : vector_node_data_labels) {
        if (requested_label == node_label) {
          int field_id = model_data.GetFieldId(node_label);
          int num_nodes = model_data.GetHostVectorNodeData(field_id).extent(0);
          // x component
          node_data_labels_.push_back(node_label + "_x");
          node_data_field_ids_.push_back(field_id);
          node_data_types_.push_back(FieldType::HostVector);
          node_data_components_.push_back(K_X);
          node_data_.push_back(std::vector<double>(num_nodes, 0.0));
          // y component
          node_data_labels_.push_back(node_label + "_y");
          node_data_field_ids_.push_back(field_id);
          node_data_types_.push_back(FieldType::HostVector);
          node_data_components_.push_back(K_Y);
          node_data_.push_back(std::vector<double>(num_nodes, 0.0));
          // z component
          node_data_labels_.push_back(node_label + "_z");
          node_data_field_ids_.push_back(field_id);
          node_data_types_.push_back(FieldType::HostVector);
          node_data_components_.push_back(K_Z);
          node_data_.push_back(std::vector<double>(num_nodes, 0.0));
        }
      }
    }

    std::vector<int> block_ids = model_data.GetBlockIds();
    for (auto const & block_id : block_ids) {
      elem_data_labels_[block_id] = std::vector<std::string>();
      elem_data_field_ids_[block_id] = std::vector<int>();
      elem_data_types_[block_id] = std::vector<FieldType>();
      elem_data_components_[block_id] = std::vector<int>();
      elem_data_[block_id] = std::vector< std::vector<double> >();
    }

    for (unsigned int i_block=0 ; i_block<block_ids.size() ; ++i_block) {

      int block_id = block_ids[i_block];

      std::vector<std::string> full_tensor_integration_point_data_labels = model_data.GetFullTensorIntegrationPointDataLabels(block_id);
      for (auto const & requested_label : requested_labels) {
        for (auto& ipt_label : full_tensor_integration_point_data_labels) {
          if (requested_label == ipt_label) {
            int field_id = model_data.GetFieldId(ipt_label);
            int num_elem = model_data.GetDeviceFullTensorIntegrationPointData(block_id, field_id, nimble::STEP_NP1).extent(0);
            // xx component
            elem_data_labels_[block_id].push_back(ipt_label + "_xx");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceFullTensor);
            elem_data_components_[block_id].push_back(K_F_XX);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // yy component
            elem_data_labels_[block_id].push_back(ipt_label + "_yy");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceFullTensor);
            elem_data_components_[block_id].push_back(K_F_YY);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // zz component
            elem_data_labels_[block_id].push_back(ipt_label + "_zz");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceFullTensor);
            elem_data_components_[block_id].push_back(K_F_ZZ);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            //  xy component
            elem_data_labels_[block_id].push_back(ipt_label + "_xy");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceFullTensor);
            elem_data_components_[block_id].push_back(K_F_XY);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // yz component
            elem_data_labels_[block_id].push_back(ipt_label + "_yz");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceFullTensor);
            elem_data_components_[block_id].push_back(K_F_YZ);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // zx component
            elem_data_labels_[block_id].push_back(ipt_label + "_zx");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceFullTensor);
            elem_data_components_[block_id].push_back(K_F_ZX);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // yx component
            elem_data_labels_[block_id].push_back(ipt_label + "_yx");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceFullTensor);
            elem_data_components_[block_id].push_back(K_F_YX);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // zy component
            elem_data_labels_[block_id].push_back(ipt_label + "_zy");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceFullTensor);
            elem_data_components_[block_id].push_back(K_F_ZY);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // xz component
            elem_data_labels_[block_id].push_back(ipt_label + "_xz");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceFullTensor);
            elem_data_components_[block_id].push_back(K_F_XZ);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          }
        }
      }

      std::vector<std::string> symmetric_tensor_integration_point_data_labels = model_data.GetSymmetricTensorIntegrationPointDataLabels(block_id);
      for (auto const & requested_label : requested_labels) {
        for (auto& ipt_label : symmetric_tensor_integration_point_data_labels) {
          if (requested_label == ipt_label) {
            int field_id = model_data.GetFieldId(ipt_label);
            int num_elem = model_data.GetDeviceSymTensorIntegrationPointData(block_id, field_id, nimble::STEP_NP1).extent(0);
            // xx component
            elem_data_labels_[block_id].push_back(ipt_label + "_xx");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceSymTensor);
            elem_data_components_[block_id].push_back(K_S_XX);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // yy component
            elem_data_labels_[block_id].push_back(ipt_label + "_yy");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceSymTensor);
            elem_data_components_[block_id].push_back(K_S_YY);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // zz component
            elem_data_labels_[block_id].push_back(ipt_label + "_zz");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceSymTensor);
            elem_data_components_[block_id].push_back(K_S_ZZ);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            //  xy component
            elem_data_labels_[block_id].push_back(ipt_label + "_xy");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceSymTensor);
            elem_data_components_[block_id].push_back(K_S_XY);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // yz component
            elem_data_labels_[block_id].push_back(ipt_label + "_yz");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceSymTensor);
            elem_data_components_[block_id].push_back(K_S_YZ);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
            // zx component
            elem_data_labels_[block_id].push_back(ipt_label + "_zx");
            elem_data_field_ids_[block_id].push_back(field_id);
            elem_data_types_[block_id].push_back(FieldType::DeviceSymTensor);
            elem_data_components_[block_id].push_back(K_S_ZX);
            elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          }
        }
      }
    }
  }

  std::vector< std::vector<double> > ExodusOutputManager::GetNodeDataForOutput(nimble_kokkos::ModelData& model_data) {
    for (unsigned int i_data=0 ; i_data < node_data_labels_.size() ; ++i_data) {
      int field_id = node_data_field_ids_.at(i_data);
      FieldType field_type = node_data_types_.at(i_data);
      if (field_type == FieldType::HostScalar) {
        HostScalarView data = model_data.GetHostScalarNodeData(field_id);
        for (unsigned int i=0 ; i<node_data_[i_data].size() ; i++) {
          node_data_[i_data][i] = data(i);
        }
      }
      else if (field_type == FieldType::HostVector) {
        HostVectorView data = model_data.GetHostVectorNodeData(field_id);
        int component = node_data_components_[i_data];
        for (unsigned int i=0 ; i<node_data_[i_data].size() ; i++) {
          node_data_[i_data][i] = data(i, component);
        }
      }
    }
    return node_data_;
  }


  std::map<int, std::vector< std::vector<double> > > ExodusOutputManager::GetElementDataForOutput(nimble_kokkos::ModelData& model_data) {

    return elem_data_;
  }

}
