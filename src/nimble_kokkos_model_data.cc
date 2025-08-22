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

#include "nimble_kokkos_model_data.h"

#include <Kokkos_ScatterView.hpp>
#include <memory>
#include <stdexcept>

#include "nimble_data_manager.h"
#include "nimble_kokkos_block.h"
#include "nimble_kokkos_material_factory.h"
#include "nimble_parser.h"
#include "nimble_vector_communicator.h"

namespace nimble_kokkos {

////// Definition of private class
class ExodusOutputManager
{
 public:
  ExodusOutputManager() : output_element_volume_(false), volume_field_id_(0){};

  void
  SpecifyOutputFields(nimble_kokkos::ModelData& model_data, std::string const& output_command_string);

  void
  ComputeElementData(
      nimble::GenesisMesh&                                      mesh,
      nimble_kokkos::ModelData&                                 model_data,
      std::map<int, nimble_kokkos::Block>&                      blocks,
      std::vector<nimble_kokkos::DeviceVectorNodeGatheredView>& gathered_reference_coordinate_d,
      std::vector<nimble_kokkos::DeviceVectorNodeGatheredView>& gathered_displacement_d);

  std::vector<std::string>
  GetNodeDataLabelsForOutput()
  {
    return node_data_labels_;
  }

  std::vector<std::vector<double>>
  GetNodeDataForOutput(nimble_kokkos::ModelData& model_data);

  std::map<int, std::vector<std::string>>
  GetElementDataLabelsForOutput()
  {
    return elem_data_labels_;
  }

  std::map<int, std::vector<std::vector<double>>>
  GetElementDataForOutput(nimble_kokkos::ModelData& model_data);

 private:
  bool output_element_volume_;

  std::vector<std::string>         node_data_labels_;
  std::vector<int>                 node_data_field_ids_;
  std::vector<FieldType>           node_data_types_;
  std::vector<int>                 node_data_components_;
  std::vector<std::vector<double>> node_data_;

  std::map<int, std::vector<std::string>>         elem_data_labels_;
  std::map<int, std::vector<int>>                 elem_data_iptdata_field_ids_;
  std::map<int, std::vector<int>>                 elem_data_edata_field_ids_;
  std::map<int, std::vector<FieldType>>           elem_data_types_;
  std::map<int, std::vector<int>>                 elem_data_integration_point_index_;
  std::map<int, std::vector<int>>                 elem_data_components_;
  std::map<int, std::vector<std::vector<double>>> elem_data_;

  int                             volume_field_id_;
  std::map<int, std::vector<int>> sym_tensor_field_ids_requiring_volume_average_;
  std::map<int, std::vector<int>> full_tensor_field_ids_requiring_volume_average_;
};

void
ExodusOutputManager::SpecifyOutputFields(nimble_kokkos::ModelData& model_data, std::string const& output_command_string)
{
  // Parse the string into individual field labels
  std::vector<std::string> requested_labels;
  std::stringstream        ss(output_command_string);
  std::string              entry;
  while (std::getline(ss, entry, ' ')) { requested_labels.push_back(entry); }

  std::vector<std::string> scalar_node_data_labels = model_data.GetScalarNodeDataLabels();
  for (auto const& requested_label : requested_labels) {
    for (auto& node_label : scalar_node_data_labels) {
      if (requested_label == node_label) {
        int field_id  = model_data.GetFieldId(node_label);
        int num_nodes = model_data.GetHostScalarNodeData(field_id).extent(0);
        node_data_labels_.push_back(node_label);
        node_data_field_ids_.push_back(field_id);
        node_data_types_.push_back(FieldType::HostScalarNode);
        node_data_components_.push_back(0);
        node_data_.push_back(std::vector<double>(num_nodes, 0.0));
      }
    }
  }

  std::vector<std::string> vector_node_data_labels = model_data.GetVectorNodeDataLabels();
  for (auto const& requested_label : requested_labels) {
    for (auto& node_label : vector_node_data_labels) {
      if (requested_label == node_label) {
        int field_id  = model_data.GetFieldId(node_label);
        int num_nodes = model_data.GetHostVectorNodeData(field_id).extent(0);
        // x component
        node_data_labels_.push_back(node_label + "_x");
        node_data_field_ids_.push_back(field_id);
        node_data_types_.push_back(FieldType::HostVectorNode);
        node_data_components_.push_back(K_X);
        node_data_.push_back(std::vector<double>(num_nodes, 0.0));
        // y component
        node_data_labels_.push_back(node_label + "_y");
        node_data_field_ids_.push_back(field_id);
        node_data_types_.push_back(FieldType::HostVectorNode);
        node_data_components_.push_back(K_Y);
        node_data_.push_back(std::vector<double>(num_nodes, 0.0));
        // z component
        node_data_labels_.push_back(node_label + "_z");
        node_data_field_ids_.push_back(field_id);
        node_data_types_.push_back(FieldType::HostVectorNode);
        node_data_components_.push_back(K_Z);
        node_data_.push_back(std::vector<double>(num_nodes, 0.0));
      }
    }
  }

  std::vector<int> block_ids = model_data.GetBlockIds();
  for (auto const& block_id : block_ids) {
    elem_data_labels_[block_id]                               = std::vector<std::string>();
    elem_data_iptdata_field_ids_[block_id]                    = std::vector<int>();
    elem_data_edata_field_ids_[block_id]                      = std::vector<int>();
    elem_data_types_[block_id]                                = std::vector<FieldType>();
    elem_data_integration_point_index_[block_id]              = std::vector<int>();
    elem_data_components_[block_id]                           = std::vector<int>();
    elem_data_[block_id]                                      = std::vector<std::vector<double>>();
    sym_tensor_field_ids_requiring_volume_average_[block_id]  = std::vector<int>();
    full_tensor_field_ids_requiring_volume_average_[block_id] = std::vector<int>();
  }

  for (unsigned int i_block = 0; i_block < block_ids.size(); ++i_block) {
    int block_id = block_ids[i_block];

    // volume is special, it is not tied to integration point data
    for (auto const& requested_label : requested_labels) {
      if (requested_label == "volume") {
        output_element_volume_ = true;
        volume_field_id_       = model_data.GetFieldId(requested_label);
        int num_elem           = model_data.GetDeviceScalarElementData(block_id, volume_field_id_).extent(0);
        elem_data_labels_[block_id].push_back(requested_label);
        elem_data_iptdata_field_ids_[block_id].push_back(volume_field_id_);
        elem_data_edata_field_ids_[block_id].push_back(volume_field_id_);
        elem_data_types_[block_id].push_back(FieldType::HostScalarElem);
        elem_data_integration_point_index_[block_id].push_back(-1);
        elem_data_components_[block_id].push_back(0);
        elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
      }
    }

    std::vector<std::string> symmetric_tensor_integration_point_data_labels =
        model_data.GetSymmetricTensorIntegrationPointDataLabels(block_id);
    for (auto const& requested_label : requested_labels) {
      for (auto& ipt_label : symmetric_tensor_integration_point_data_labels) {
        std::string field_label              = requested_label;
        bool        single_integration_point = false;
        bool        require_volume_average   = true;
        int         integration_point_index  = -1;
        if (requested_label.substr(0, 3) == "ipt") {
          field_label              = requested_label.substr(6, requested_label.size() - 6);
          single_integration_point = true;
          require_volume_average   = false;
          integration_point_index  = atoi(requested_label.substr(4, 2).c_str()) - 1;
        }
        if (field_label == ipt_label) {
          int iptdata_field_id = model_data.GetFieldId(ipt_label);
          int edata_field_id   = iptdata_field_id;
          if (require_volume_average) {
            sym_tensor_field_ids_requiring_volume_average_[block_id].push_back(iptdata_field_id);
          }
          int num_elem =
              model_data.GetDeviceSymTensorIntegrationPointData(block_id, iptdata_field_id, nimble::STEP_NP1).extent(0);
          if (single_integration_point) {
            edata_field_id =
                model_data.AllocateElementData(block_id, nimble::SYMMETRIC_TENSOR, requested_label, num_elem);
          }
          // xx component
          elem_data_labels_[block_id].push_back(requested_label + "_xx");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostSymTensorElem);
          elem_data_components_[block_id].push_back(K_S_XX);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // yy component
          elem_data_labels_[block_id].push_back(requested_label + "_yy");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostSymTensorElem);
          elem_data_components_[block_id].push_back(K_S_YY);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // zz component
          elem_data_labels_[block_id].push_back(requested_label + "_zz");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostSymTensorElem);
          elem_data_components_[block_id].push_back(K_S_ZZ);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          //  xy component
          elem_data_labels_[block_id].push_back(requested_label + "_xy");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostSymTensorElem);
          elem_data_components_[block_id].push_back(K_S_XY);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // yz component
          elem_data_labels_[block_id].push_back(requested_label + "_yz");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostSymTensorElem);
          elem_data_components_[block_id].push_back(K_S_YZ);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // zx component
          elem_data_labels_[block_id].push_back(requested_label + "_zx");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostSymTensorElem);
          elem_data_components_[block_id].push_back(K_S_ZX);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
        }
      }
    }

    std::vector<std::string> full_tensor_integration_point_data_labels =
        model_data.GetFullTensorIntegrationPointDataLabels(block_id);
    for (auto const& requested_label : requested_labels) {
      // Handle case where user has requested output for a specific integration
      // point
      std::string field_label              = requested_label;
      bool        single_integration_point = false;
      bool        require_volume_average   = true;
      int         integration_point_index  = -1;
      if (requested_label.substr(0, 3) == "ipt") {
        field_label              = requested_label.substr(6, requested_label.size() - 6);
        single_integration_point = true;
        require_volume_average   = false;
        integration_point_index  = atoi(requested_label.substr(4, 2).c_str()) - 1;
      }
      for (auto& ipt_label : full_tensor_integration_point_data_labels) {
        if (field_label == ipt_label) {
          int iptdata_field_id = model_data.GetFieldId(ipt_label);
          int edata_field_id   = iptdata_field_id;
          if (require_volume_average) {
            full_tensor_field_ids_requiring_volume_average_[block_id].push_back(iptdata_field_id);
          }
          int num_elem =
              model_data.GetDeviceFullTensorIntegrationPointData(block_id, iptdata_field_id, nimble::STEP_NP1)
                  .extent(0);
          if (single_integration_point) {
            edata_field_id = model_data.AllocateElementData(block_id, nimble::FULL_TENSOR, requested_label, num_elem);
          }
          // xx component
          elem_data_labels_[block_id].push_back(requested_label + "_xx");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostFullTensorElem);
          elem_data_components_[block_id].push_back(K_F_XX);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // yy component
          elem_data_labels_[block_id].push_back(requested_label + "_yy");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostFullTensorElem);
          elem_data_components_[block_id].push_back(K_F_YY);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // zz component
          elem_data_labels_[block_id].push_back(requested_label + "_zz");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostFullTensorElem);
          elem_data_components_[block_id].push_back(K_F_ZZ);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          //  xy component
          elem_data_labels_[block_id].push_back(requested_label + "_xy");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostFullTensorElem);
          elem_data_components_[block_id].push_back(K_F_XY);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // yz component
          elem_data_labels_[block_id].push_back(requested_label + "_yz");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostFullTensorElem);
          elem_data_components_[block_id].push_back(K_F_YZ);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // zx component
          elem_data_labels_[block_id].push_back(requested_label + "_zx");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostFullTensorElem);
          elem_data_components_[block_id].push_back(K_F_ZX);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // yx component
          elem_data_labels_[block_id].push_back(requested_label + "_yx");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostFullTensorElem);
          elem_data_components_[block_id].push_back(K_F_YX);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // zy component
          elem_data_labels_[block_id].push_back(requested_label + "_zy");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostFullTensorElem);
          elem_data_components_[block_id].push_back(K_F_ZY);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
          // xz component
          elem_data_labels_[block_id].push_back(requested_label + "_xz");
          elem_data_iptdata_field_ids_[block_id].push_back(iptdata_field_id);
          elem_data_edata_field_ids_[block_id].push_back(edata_field_id);
          elem_data_types_[block_id].push_back(FieldType::HostFullTensorElem);
          elem_data_components_[block_id].push_back(K_F_XZ);
          elem_data_integration_point_index_[block_id].push_back(integration_point_index);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
        }
      }
    }
  }
}

void
ExodusOutputManager::ComputeElementData(
    nimble::GenesisMesh&                                      mesh,
    nimble_kokkos::ModelData&                                 model_data,
    std::map<int, nimble_kokkos::Block>&                      blocks,
    std::vector<nimble_kokkos::DeviceVectorNodeGatheredView>& gathered_reference_coordinate_d,
    std::vector<nimble_kokkos::DeviceVectorNodeGatheredView>& gathered_displacement_d)
{
  int                                           block_index;
  std::map<int, nimble_kokkos::Block>::iterator block_it;

  // Element volume
  if (output_element_volume_) {
    for (block_index = 0, block_it = blocks.begin(); block_it != blocks.end(); block_index++, block_it++) {
      int                                          block_id           = block_it->first;
      nimble_kokkos::Block&                        block              = block_it->second;
      nimble::Element*                             element_d          = block.GetDeviceElement();
      int                                          num_elem_in_block  = mesh.GetNumElementsInBlock(block_id);
      int                                          num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);
      nimble_kokkos::DeviceElementConnectivityView elem_conn_d        = block.GetDeviceElementConnectivityView();
      nimble_kokkos::DeviceVectorNodeGatheredView  gathered_reference_coordinate_block_d =
          gathered_reference_coordinate_d.at(block_index);
      nimble_kokkos::DeviceVectorNodeGatheredView gathered_displacement_block_d =
          gathered_displacement_d.at(block_index);
      nimble_kokkos::DeviceScalarElemView volume_d = model_data.GetDeviceScalarElementData(block_id, volume_field_id_);
      Kokkos::parallel_for(
          "Element Volume", num_elem_in_block, KOKKOS_LAMBDA(const int i_elem) {
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d =
                Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d =
                Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceScalarElemSingleEntryView element_volume_d = Kokkos::subview(volume_d, i_elem);
            element_d->ComputeVolume(element_reference_coordinate_d, element_displacement_d, element_volume_d);
          });
      nimble_kokkos::HostScalarElemView volume_h = model_data.GetHostScalarElementData(block_id, volume_field_id_);
      deep_copy(volume_h, volume_d);
    }
  }

  // Extract data for specific integration points
  for (block_index = 0, block_it = blocks.begin(); block_it != blocks.end(); block_index++, block_it++) {
    int block_id          = block_it->first;
    int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    for (unsigned int i_data = 0; i_data < elem_data_labels_.at(block_id).size(); ++i_data) {
      int integration_point_index = elem_data_integration_point_index_.at(block_id).at(i_data);
      if (integration_point_index != -1) {
        int       iptdata_field_id = elem_data_iptdata_field_ids_.at(block_id).at(i_data);
        int       edata_field_id   = elem_data_edata_field_ids_.at(block_id).at(i_data);
        FieldType field_type       = elem_data_types_.at(block_id).at(i_data);
        if (field_type == FieldType::HostSymTensorElem) {
          nimble_kokkos::DeviceSymTensorIntPtView sym_tensor_data_step_np1_d =
              model_data.GetDeviceSymTensorIntegrationPointData(block_id, iptdata_field_id, nimble::STEP_NP1);
          nimble_kokkos::DeviceSymTensorElemView sym_tensor_data_single_int_pt_d =
              model_data.GetDeviceSymTensorElementData(block_id, edata_field_id);
          Kokkos::parallel_for(
              "Extract Symmetric Tensor Integration Point Data for Output",
              num_elem_in_block,
              KOKKOS_LAMBDA(const int i_elem) {
                nimble_kokkos::DeviceSymTensorIntPtSubView element_sym_tensor_step_np1_d =
                    Kokkos::subview(sym_tensor_data_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
                nimble_kokkos::DeviceSymTensorElemSingleEntryView element_sym_tensor_single_int_pt_d =
                    Kokkos::subview(sym_tensor_data_single_int_pt_d, i_elem, Kokkos::ALL);
                for (int i = 0; i < element_sym_tensor_single_int_pt_d.extent(0); i++) {
                  element_sym_tensor_single_int_pt_d(i) = element_sym_tensor_step_np1_d(integration_point_index, i);
                }
              });
          nimble_kokkos::HostSymTensorElemView sym_tensor_data_single_int_pt_h =
              model_data.GetHostSymTensorElementData(block_id, edata_field_id);
          deep_copy(sym_tensor_data_single_int_pt_h, sym_tensor_data_single_int_pt_d);

        } else if (field_type == FieldType::HostFullTensorElem) {
          nimble_kokkos::DeviceFullTensorIntPtView full_tensor_data_step_np1_d =
              model_data.GetDeviceFullTensorIntegrationPointData(block_id, iptdata_field_id, nimble::STEP_NP1);
          nimble_kokkos::DeviceFullTensorElemView full_tensor_data_single_int_pt_d =
              model_data.GetDeviceFullTensorElementData(block_id, edata_field_id);
          Kokkos::parallel_for(
              "Extract Full Tensor Integration Point Data for Output",
              num_elem_in_block,
              KOKKOS_LAMBDA(const int i_elem) {
                nimble_kokkos::DeviceFullTensorIntPtSubView element_full_tensor_step_np1_d =
                    Kokkos::subview(full_tensor_data_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
                nimble_kokkos::DeviceFullTensorElemSingleEntryView element_full_tensor_single_int_pt_d =
                    Kokkos::subview(full_tensor_data_single_int_pt_d, i_elem, Kokkos::ALL);
                for (int i = 0; i < element_full_tensor_single_int_pt_d.extent(0); i++) {
                  element_full_tensor_single_int_pt_d(i) = element_full_tensor_step_np1_d(integration_point_index, i);
                }
              });
          nimble_kokkos::HostFullTensorElemView full_tensor_data_single_int_pt_h =
              model_data.GetHostFullTensorElementData(block_id, edata_field_id);
          deep_copy(full_tensor_data_single_int_pt_h, full_tensor_data_single_int_pt_d);
        }
      }
    }
  }

  // Volume averaging of symmetric and full tensors stored at integration points
  for (block_index = 0, block_it = blocks.begin(); block_it != blocks.end(); block_index++, block_it++) {
    int                                          block_id           = block_it->first;
    nimble_kokkos::Block&                        block              = block_it->second;
    nimble::Element*                             element_d          = block.GetDeviceElement();
    int                                          num_elem_in_block  = mesh.GetNumElementsInBlock(block_id);
    int                                          num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);
    nimble_kokkos::DeviceElementConnectivityView elem_conn_d        = block.GetDeviceElementConnectivityView();
    nimble_kokkos::DeviceVectorNodeGatheredView  gathered_reference_coordinate_block_d =
        gathered_reference_coordinate_d.at(block_index);
    nimble_kokkos::DeviceVectorNodeGatheredView gathered_displacement_block_d = gathered_displacement_d.at(block_index);
    for (auto& field_id : sym_tensor_field_ids_requiring_volume_average_.at(block_id)) {
      nimble_kokkos::DeviceSymTensorIntPtView sym_tensor_data_step_np1_d =
          model_data.GetDeviceSymTensorIntegrationPointData(block_id, field_id, nimble::STEP_NP1);
      nimble_kokkos::DeviceSymTensorElemView sym_tensor_data_vol_ave_d =
          model_data.GetDeviceSymTensorElementData(block_id, field_id);
      Kokkos::parallel_for(
          "Volume Averaging Sym Tensor", num_elem_in_block, KOKKOS_LAMBDA(const int i_elem) {
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d =
                Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d =
                Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceSymTensorIntPtSubView element_sym_tensor_step_np1_d =
                Kokkos::subview(sym_tensor_data_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceSymTensorElemSingleEntryView element_sym_tensor_vol_ave_d =
                Kokkos::subview(sym_tensor_data_vol_ave_d, i_elem, Kokkos::ALL);
            element_d->ComputeVolumeAverageSymTensor(
                element_reference_coordinate_d,
                element_displacement_d,
                element_sym_tensor_step_np1_d,
                element_sym_tensor_vol_ave_d);
          });
      nimble_kokkos::HostSymTensorElemView sym_tensor_data_vol_ave_h =
          model_data.GetHostSymTensorElementData(block_id, field_id);
      deep_copy(sym_tensor_data_vol_ave_h, sym_tensor_data_vol_ave_d);
    }
    for (auto& field_id : full_tensor_field_ids_requiring_volume_average_.at(block_id)) {
      nimble_kokkos::DeviceFullTensorIntPtView full_tensor_data_step_np1_d =
          model_data.GetDeviceFullTensorIntegrationPointData(block_id, field_id, nimble::STEP_NP1);
      nimble_kokkos::DeviceFullTensorElemView full_tensor_data_vol_ave_d =
          model_data.GetDeviceFullTensorElementData(block_id, field_id);
      Kokkos::parallel_for(
          "Volume Averaging Full Tensor", num_elem_in_block, KOKKOS_LAMBDA(const int i_elem) {
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d =
                Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d =
                Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceFullTensorIntPtSubView element_full_tensor_step_np1_d =
                Kokkos::subview(full_tensor_data_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceFullTensorElemSingleEntryView element_full_tensor_vol_ave_d =
                Kokkos::subview(full_tensor_data_vol_ave_d, i_elem, Kokkos::ALL);
            element_d->ComputeVolumeAverageFullTensor(
                element_reference_coordinate_d,
                element_displacement_d,
                element_full_tensor_step_np1_d,
                element_full_tensor_vol_ave_d);
          });
      nimble_kokkos::HostFullTensorElemView full_tensor_data_vol_ave_h =
          model_data.GetHostFullTensorElementData(block_id, field_id);
      deep_copy(full_tensor_data_vol_ave_h, full_tensor_data_vol_ave_d);
    }
  }
}

std::vector<std::vector<double>>
ExodusOutputManager::GetNodeDataForOutput(nimble_kokkos::ModelData& model_data)
{
  for (unsigned int i_data = 0; i_data < node_data_labels_.size(); ++i_data) {
    int       field_id   = node_data_field_ids_.at(i_data);
    FieldType field_type = node_data_types_.at(i_data);
    if (field_type == FieldType::HostScalarNode) {
      HostScalarNodeView data = model_data.GetHostScalarNodeData(field_id);
      for (unsigned int i = 0; i < node_data_[i_data].size(); i++) { node_data_[i_data][i] = data(i); }
    } else if (field_type == FieldType::HostVectorNode) {
      HostVectorNodeView data      = model_data.GetHostVectorNodeData(field_id);
      int                component = node_data_components_[i_data];
      for (unsigned int i = 0; i < node_data_[i_data].size(); i++) { node_data_[i_data][i] = data(i, component); }
    }
  }
  return node_data_;
}

std::map<int, std::vector<std::vector<double>>>
ExodusOutputManager::GetElementDataForOutput(nimble_kokkos::ModelData& model_data)
{
  std::vector<int> block_ids = model_data.GetBlockIds();
  for (auto const& block_id : block_ids) {
    for (unsigned int i_data = 0; i_data < elem_data_labels_.at(block_id).size(); ++i_data) {
      int       edata_field_id = elem_data_edata_field_ids_.at(block_id).at(i_data);
      FieldType field_type     = elem_data_types_.at(block_id).at(i_data);
      if (field_type == FieldType::HostScalarElem) {
        HostScalarElemView data = model_data.GetHostScalarElementData(block_id, edata_field_id);
        for (unsigned int i = 0; i < elem_data_.at(block_id)[i_data].size(); i++) {
          elem_data_.at(block_id)[i_data][i] = data(i);
        }
      } else if (field_type == FieldType::HostSymTensorElem) {
        HostSymTensorElemView data      = model_data.GetHostSymTensorElementData(block_id, edata_field_id);
        int                   component = elem_data_components_.at(block_id).at(i_data);
        for (unsigned int i = 0; i < elem_data_.at(block_id)[i_data].size(); i++) {
          elem_data_.at(block_id)[i_data][i] = data(i, component);
        }
      } else if (field_type == FieldType::HostFullTensorElem) {
        HostFullTensorElemView data      = model_data.GetHostFullTensorElementData(block_id, edata_field_id);
        int                    component = elem_data_components_.at(block_id).at(i_data);
        for (unsigned int i = 0; i < elem_data_.at(block_id)[i_data].size(); i++) {
          elem_data_.at(block_id)[i_data][i] = data(i, component);
        }
      }
    }
  }
  return elem_data_;
}
//////////////////////////////////////

ModelData::ModelData() : exodus_output_manager_(std::unique_ptr<ExodusOutputManager>(new ExodusOutputManager())) {}

ModelData::~ModelData() {}

int
ModelData::AllocateNodeData(nimble::Length length, std::string label, int num_objects)
{
  int  field_id;
  auto it = field_label_to_field_id_map_.find(label);
  if (it == field_label_to_field_id_map_.end()) {
    field_id                            = field_label_to_field_id_map_.size();
    field_label_to_field_id_map_[label] = field_id;
  } else {
    field_id = it->second;
  }

  if (length == nimble::SCALAR) {
    // device_node_data_ is of type std::vector< std::unique_ptr< FieldBase > >
    device_node_data_.emplace_back(new Field<FieldType::DeviceScalarNode>(label, num_objects));
  } else if (length == nimble::VECTOR) {
    device_node_data_.emplace_back(new Field<FieldType::DeviceVectorNode>(label, num_objects));
  } else {
    throw std::invalid_argument(
        "\nError:  Invalid device data length in "
        "nimble_kokkos::ModelData::AllocateNodeData().\n");
  }

  field_id_to_device_node_data_index_[field_id] = device_node_data_.size() - 1;

  FieldBase* d_field = device_node_data_.back().get();

  if (d_field->type() == FieldType::DeviceScalarNode) {
    auto                                     field  = dynamic_cast<Field<FieldType::DeviceScalarNode>*>(d_field);
    Field<FieldType::DeviceScalarNode>::View d_view = field->data();
    auto                                     h_view = Kokkos::create_mirror_view(d_view);
    host_node_data_.emplace_back(new Field<FieldType::HostScalarNode>(h_view));
  } else if (d_field->type() == FieldType::DeviceVectorNode) {
    auto                                     field  = dynamic_cast<Field<FieldType::DeviceVectorNode>*>(d_field);
    Field<FieldType::DeviceVectorNode>::View d_view = field->data();
    auto                                     h_view = Kokkos::create_mirror_view(d_view);
    host_node_data_.emplace_back(new Field<FieldType::HostVectorNode>(h_view));
  } else {
    throw std::invalid_argument(
        "\nError:  Invalid host data length in "
        "nimble_kokkos::ModelData::AllocateNodeData().\n");
  }

  field_id_to_host_node_data_index_[field_id] = host_node_data_.size() - 1;

  return field_id;
}

int
ModelData::AllocateElementData(int block_id, nimble::Length length, std::string label, int num_objects)
{
  int  field_id;
  auto it = field_label_to_field_id_map_.find(label);
  if (it == field_label_to_field_id_map_.end()) {
    field_id                            = field_label_to_field_id_map_.size();
    field_label_to_field_id_map_[label] = field_id;
  } else {
    field_id = it->second;
  }

  if (block_id_to_element_data_index_.find(block_id) == block_id_to_element_data_index_.end()) {
    block_id_to_element_data_index_[block_id] = host_element_data_.size();
    host_element_data_.emplace_back();
    device_element_data_.emplace_back();
    field_id_to_host_element_data_index_.emplace_back();
    field_id_to_device_element_data_index_.emplace_back();
  }
  int block_index = block_id_to_element_data_index_.at(block_id);

  if (length == nimble::SCALAR) {
    device_element_data_.at(block_index).emplace_back(new Field<FieldType::DeviceScalarElem>(label, num_objects));
  } else if (length == nimble::SYMMETRIC_TENSOR) {
    device_element_data_.at(block_index).emplace_back(new Field<FieldType::DeviceSymTensorElem>(label, num_objects));
  } else if (length == nimble::FULL_TENSOR) {
    device_element_data_.at(block_index).emplace_back(new Field<FieldType::DeviceFullTensorElem>(label, num_objects));
  } else {
    throw std::invalid_argument(
        "\nError:  Invalid device data length in "
        "nimble_kokkos::ModelData::AllocateElementData().\n");
  }

  field_id_to_device_element_data_index_.at(block_index)[field_id] = device_element_data_.at(block_index).size() - 1;

  FieldBase* d_field = device_element_data_.at(block_index).back().get();

  if (d_field->type() == FieldType::DeviceScalarElem) {
    auto                                     field  = dynamic_cast<Field<FieldType::DeviceScalarElem>*>(d_field);
    Field<FieldType::DeviceScalarElem>::View d_view = field->data();
    auto                                     h_view = Kokkos::create_mirror_view(d_view);
    host_element_data_.at(block_index).emplace_back(new Field<FieldType::HostScalarElem>(h_view));
  } else if (d_field->type() == FieldType::DeviceSymTensorElem) {
    auto                                        field  = dynamic_cast<Field<FieldType::DeviceSymTensorElem>*>(d_field);
    Field<FieldType::DeviceSymTensorElem>::View d_view = field->data();
    auto                                        h_view = Kokkos::create_mirror_view(d_view);
    host_element_data_.at(block_index).emplace_back(new Field<FieldType::HostSymTensorElem>(h_view));
  } else if (d_field->type() == FieldType::DeviceFullTensorElem) {
    auto                                         field = dynamic_cast<Field<FieldType::DeviceFullTensorElem>*>(d_field);
    Field<FieldType::DeviceFullTensorElem>::View d_view = field->data();
    auto                                         h_view = Kokkos::create_mirror_view(d_view);
    host_element_data_.at(block_index).emplace_back(new Field<FieldType::HostFullTensorElem>(h_view));
  } else {
    throw std::invalid_argument(
        "\nError:  Invalid host data length in "
        "nimble_kokkos::ModelData::AllocateElementData().\n");
  }

  field_id_to_host_element_data_index_.at(block_index)[field_id] = host_element_data_.at(block_index).size() - 1;

  return field_id;
}

int
ModelData::AllocateIntegrationPointData(
    int                 block_id,
    nimble::Length      length,
    std::string         label,
    int                 num_objects,
    std::vector<double> initial_value)
{
  bool set_initial_value = false;
  if (!initial_value.empty()) set_initial_value = true;

  int  field_id;
  auto it = field_label_to_field_id_map_.find(label);
  if (it == field_label_to_field_id_map_.end()) {
    field_id                            = field_label_to_field_id_map_.size();
    field_label_to_field_id_map_[label] = field_id;
  } else {
    field_id = it->second;
  }

  if (block_id_to_integration_point_data_index_.find(block_id) == block_id_to_integration_point_data_index_.end()) {
    block_id_to_integration_point_data_index_[block_id] = host_integration_point_data_step_n_.size();
    host_integration_point_data_step_n_.emplace_back();
    host_integration_point_data_step_np1_.emplace_back();
    device_integration_point_data_step_n_.emplace_back();
    device_integration_point_data_step_np1_.emplace_back();
    field_id_to_host_integration_point_data_index_.emplace_back();
    field_id_to_device_integration_point_data_index_.emplace_back();
  }
  int block_index = block_id_to_integration_point_data_index_.at(block_id);

  if (length == nimble::SCALAR) {
    device_integration_point_data_step_n_.at(block_index)
        .emplace_back(new Field<FieldType::DeviceScalarNode>(label, num_objects));
    device_integration_point_data_step_np1_.at(block_index)
        .emplace_back(new Field<FieldType::DeviceScalarNode>(label, num_objects));
  } else if (length == nimble::VECTOR) {
    device_integration_point_data_step_n_.at(block_index)
        .emplace_back(new Field<FieldType::DeviceVectorNode>(label, num_objects));
    device_integration_point_data_step_np1_.at(block_index)
        .emplace_back(new Field<FieldType::DeviceVectorNode>(label, num_objects));
  } else if (length == nimble::SYMMETRIC_TENSOR) {
    device_integration_point_data_step_n_.at(block_index)
        .emplace_back(new Field<FieldType::DeviceSymTensorIntPt>(label, num_objects));
    device_integration_point_data_step_np1_.at(block_index)
        .emplace_back(new Field<FieldType::DeviceSymTensorIntPt>(label, num_objects));
  } else if (length == nimble::FULL_TENSOR) {
    device_integration_point_data_step_n_.at(block_index)
        .emplace_back(new Field<FieldType::DeviceFullTensorIntPt>(label, num_objects));
    device_integration_point_data_step_np1_.at(block_index)
        .emplace_back(new Field<FieldType::DeviceFullTensorIntPt>(label, num_objects));
  } else {
    throw std::invalid_argument(
        "\nError:  Invalid device data length in "
        "nimble_kokkos::ModelData::AllocateIntegrationPointData().\n");
  }

  field_id_to_device_integration_point_data_index_.at(block_index)[field_id] =
      device_integration_point_data_step_n_.at(block_index).size() - 1;

  FieldBase* d_field_step_n   = device_integration_point_data_step_n_.at(block_index).back().get();
  FieldBase* d_field_step_np1 = device_integration_point_data_step_np1_.at(block_index).back().get();

  if (d_field_step_n->type() == FieldType::DeviceScalarNode) {
    auto field_step_n = dynamic_cast<Field<FieldType::DeviceScalarNode>*>(d_field_step_n);
    Field<FieldType::DeviceScalarNode>::View d_view_step_n = field_step_n->data();
    auto                                     h_view_step_n = Kokkos::create_mirror_view(d_view_step_n);
    host_integration_point_data_step_n_.at(block_index)
        .emplace_back(new Field<FieldType::HostScalarNode>(h_view_step_n));

    auto field_step_np1 = dynamic_cast<Field<FieldType::DeviceScalarNode>*>(d_field_step_np1);
    Field<FieldType::DeviceScalarNode>::View d_view_step_np1 = field_step_np1->data();
    auto                                     h_view_step_np1 = Kokkos::create_mirror_view(d_view_step_np1);
    host_integration_point_data_step_np1_.at(block_index)
        .emplace_back(new Field<FieldType::HostScalarNode>(h_view_step_np1));

    if (set_initial_value) {
      int num_elem = h_view_step_n.extent(0);
      for (int i_elem = 0; i_elem < num_elem; ++i_elem) {
        h_view_step_n(i_elem)   = initial_value.at(0);
        h_view_step_np1(i_elem) = initial_value.at(0);
      }
      Kokkos::deep_copy(d_view_step_n, h_view_step_n);
      Kokkos::deep_copy(d_view_step_np1, h_view_step_np1);
    }
  } else if (d_field_step_n->type() == FieldType::DeviceVectorNode) {
    auto field_step_n = dynamic_cast<Field<FieldType::DeviceVectorNode>*>(d_field_step_n);
    Field<FieldType::DeviceVectorNode>::View d_view_step_n = field_step_n->data();
    auto                                     h_view_step_n = Kokkos::create_mirror_view(d_view_step_n);
    host_integration_point_data_step_n_.at(block_index)
        .emplace_back(new Field<FieldType::HostVectorNode>(h_view_step_n));

    auto field_step_np1 = dynamic_cast<Field<FieldType::DeviceVectorNode>*>(d_field_step_np1);
    Field<FieldType::DeviceVectorNode>::View d_view_step_np1 = field_step_np1->data();
    auto                                     h_view_step_np1 = Kokkos::create_mirror_view(d_view_step_np1);
    host_integration_point_data_step_np1_.at(block_index)
        .emplace_back(new Field<FieldType::HostVectorNode>(h_view_step_np1));

    if (set_initial_value) {
      int num_elem    = h_view_step_n.extent(0);
      int num_entries = 3;
      for (int i_elem = 0; i_elem < num_elem; ++i_elem) {
        for (int i_entry = 0; i_entry < num_entries; ++i_entry) {
          h_view_step_n(i_elem, i_entry)   = initial_value.at(i_entry);
          h_view_step_np1(i_elem, i_entry) = initial_value.at(i_entry);
        }
      }
      Kokkos::deep_copy(d_view_step_n, h_view_step_n);
      Kokkos::deep_copy(d_view_step_np1, h_view_step_np1);
    }
  } else if (d_field_step_n->type() == FieldType::DeviceSymTensorIntPt) {
    auto field_step_n = dynamic_cast<Field<FieldType::DeviceSymTensorIntPt>*>(d_field_step_n);
    Field<FieldType::DeviceSymTensorIntPt>::View d_view_step_n = field_step_n->data();
    auto                                         h_view_step_n = Kokkos::create_mirror_view(d_view_step_n);
    host_integration_point_data_step_n_.at(block_index)
        .emplace_back(new Field<FieldType::HostSymTensorIntPt>(h_view_step_n));

    auto field_step_np1 = dynamic_cast<Field<FieldType::DeviceSymTensorIntPt>*>(d_field_step_np1);
    Field<FieldType::DeviceSymTensorIntPt>::View d_view_step_np1 = field_step_np1->data();
    auto                                         h_view_step_np1 = Kokkos::create_mirror_view(d_view_step_np1);
    host_integration_point_data_step_np1_.at(block_index)
        .emplace_back(new Field<FieldType::HostSymTensorIntPt>(h_view_step_np1));

    if (set_initial_value) {
      int num_elem    = h_view_step_n.extent(0);
      int num_int_pt  = h_view_step_n.extent(1);
      int num_entries = 6;
      for (int i_elem = 0; i_elem < num_elem; ++i_elem) {
        for (int i_int_pt = 0; i_int_pt < num_int_pt; ++i_int_pt) {
          for (int i_entry = 0; i_entry < num_entries; ++i_entry) {
            h_view_step_n(i_elem, i_int_pt, i_entry)   = initial_value.at(i_entry);
            h_view_step_np1(i_elem, i_int_pt, i_entry) = initial_value.at(i_entry);
          }
        }
      }
      Kokkos::deep_copy(d_view_step_n, h_view_step_n);
      Kokkos::deep_copy(d_view_step_np1, h_view_step_np1);
    }
  } else if (d_field_step_n->type() == FieldType::DeviceFullTensorIntPt) {
    auto field_step_n = dynamic_cast<Field<FieldType::DeviceFullTensorIntPt>*>(d_field_step_n);
    Field<FieldType::DeviceFullTensorIntPt>::View d_view_step_n = field_step_n->data();
    auto                                          h_view_step_n = Kokkos::create_mirror_view(d_view_step_n);
    host_integration_point_data_step_n_.at(block_index)
        .emplace_back(new Field<FieldType::HostFullTensorIntPt>(h_view_step_n));

    auto field_step_np1 = dynamic_cast<Field<FieldType::DeviceFullTensorIntPt>*>(d_field_step_np1);
    Field<FieldType::DeviceFullTensorIntPt>::View d_view_step_np1 = field_step_np1->data();
    auto                                          h_view_step_np1 = Kokkos::create_mirror_view(d_view_step_np1);
    host_integration_point_data_step_np1_.at(block_index)
        .emplace_back(new Field<FieldType::HostFullTensorIntPt>(h_view_step_np1));

    if (set_initial_value) {
      int num_elem    = h_view_step_n.extent(0);
      int num_int_pt  = h_view_step_n.extent(1);
      int num_entries = 9;
      for (int i_elem = 0; i_elem < num_elem; ++i_elem) {
        for (int i_int_pt = 0; i_int_pt < num_int_pt; ++i_int_pt) {
          for (int i_entry = 0; i_entry < num_entries; ++i_entry) {
            h_view_step_n(i_elem, i_int_pt, i_entry)   = initial_value.at(i_entry);
            h_view_step_np1(i_elem, i_int_pt, i_entry) = initial_value.at(i_entry);
          }
        }
      }
      Kokkos::deep_copy(d_view_step_n, h_view_step_n);
      Kokkos::deep_copy(d_view_step_np1, h_view_step_np1);
    }
  } else {
    throw std::invalid_argument(
        "\nError:  Invalid host data length in "
        "nimble_kokkos::ModelData::AllocateElementData().\n");
  }

  field_id_to_host_integration_point_data_index_.at(block_index)[field_id] =
      host_integration_point_data_step_n_.at(block_index).size() - 1;

  return field_id;
}

std::vector<int>
ModelData::GetBlockIds() const
{
  std::vector<int> block_ids;
  for (auto const& entry : block_id_to_integration_point_data_index_) { block_ids.push_back(entry.first); }
  return block_ids;
}

void
ModelData::InitializeBlocks(
    nimble::DataManager&                                data_manager,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base)
{
  bool store_unrotated_stress(true);

  const auto& mesh_      = data_manager.GetMesh();
  const auto& parser_    = data_manager.GetParser();
  auto&       field_ids_ = data_manager.GetFieldIDs();

  auto       material_factory_ptr = dynamic_cast<nimble_kokkos::MaterialFactory*>(material_factory_base.get());
  const auto num_blocks           = static_cast<int>(mesh_.GetNumBlocks());

  //
  // Blocks
  //
  std::vector<int> block_ids = mesh_.GetBlockIds();
  for (int i = 0; i < num_blocks; i++) {
    int                block_id                  = block_ids.at(i);
    std::string const& model_material_parameters = parser_.GetModelMaterialParameters(block_id);
    int                num_elements_in_block     = mesh_.GetNumElementsInBlock(block_id);
    blocks_[block_id]                            = nimble_kokkos::Block();
    blocks_.at(block_id).Initialize(model_material_parameters, num_elements_in_block, *material_factory_ptr);
    //
    // MPI version use model_data.DeclareElementData(block_id,
    // data_labels_and_lengths);
    //
    std::vector<double> initial_value(9, 0.0);
    initial_value[0] = initial_value[1] = initial_value[2] = 1.0;
    field_ids_.deformation_gradient                        = AllocateIntegrationPointData(
        block_id, nimble::FULL_TENSOR, "deformation_gradient", num_elements_in_block, initial_value);
    // volume-averaged quantities for I/O are stored as element data
    AllocateElementData(block_id, nimble::FULL_TENSOR, "deformation_gradient", num_elements_in_block);

    field_ids_.stress =
        AllocateIntegrationPointData(block_id, nimble::SYMMETRIC_TENSOR, "stress", num_elements_in_block);
    if (store_unrotated_stress) {
      field_ids_.unrotated_stress =
          AllocateIntegrationPointData(block_id, nimble::SYMMETRIC_TENSOR, "stress", num_elements_in_block);
    }

    // volume-averaged quantities for I/O are stored as element data
    AllocateElementData(block_id, nimble::SYMMETRIC_TENSOR, "stress", num_elements_in_block);

    if (parser_.GetOutputFieldString().find("volume") != std::string::npos) {
      AllocateElementData(block_id, nimble::SCALAR, "volume", num_elements_in_block);
    }

    if (blocks_.at(block_id).GetMaterialPointer()->NumStateVariables() > 0) {
      //
      // Verify which state variables are actually needed
      //
      field_ids_.state_sym_tensor =
          AllocateIntegrationPointData(block_id, nimble::SYMMETRIC_TENSOR, "state_sym_tensor", num_elements_in_block);
      field_ids_.state_full_tensor =
          AllocateIntegrationPointData(block_id, nimble::FULL_TENSOR, "state_full_tensor", num_elements_in_block);
      field_ids_.state_scalar =
          AllocateIntegrationPointData(block_id, nimble::SCALAR, "state_scalar", num_elements_in_block);
      field_ids_.state_vec3D =
          AllocateIntegrationPointData(block_id, nimble::VECTOR, "state_vec3D", num_elements_in_block);
    }
  }

  // Initialize gathered containers when using explicit scheme
  if (parser_.TimeIntegrationScheme() == "explicit") InitializeGatheredVectors(mesh_);

  InitializeBlockData(data_manager);
}

void
ModelData::UpdateStates(const nimble::DataManager& data_manager)
{
  const auto& field_ids_ = data_manager.GetFieldIDs();

  // Copy STEP_NP1 data to STEP_N
  int block_index = 0;
  for (auto& block_it : blocks_) {
    int  block_id = block_it.first;
    auto deformation_gradient_step_n_d =
        GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.deformation_gradient, nimble::STEP_N);
    auto unrotated_stress_step_n_d =
        GetDeviceSymTensorIntegrationPointData(block_id, field_ids_.unrotated_stress, nimble::STEP_N);
    auto stress_step_n_d = GetDeviceSymTensorIntegrationPointData(block_id, field_ids_.stress, nimble::STEP_N);
    auto deformation_gradient_step_np1_d =
        GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.deformation_gradient, nimble::STEP_NP1);
    auto unrotated_stress_step_np1_d =
        GetDeviceSymTensorIntegrationPointData(block_id, field_ids_.unrotated_stress, nimble::STEP_NP1);
    auto stress_step_np1_d = GetDeviceSymTensorIntegrationPointData(block_id, field_ids_.stress, nimble::STEP_NP1);
    Kokkos::deep_copy(deformation_gradient_step_n_d, deformation_gradient_step_np1_d);
    Kokkos::deep_copy(unrotated_stress_step_n_d, unrotated_stress_step_np1_d);
    Kokkos::deep_copy(stress_step_n_d, stress_step_np1_d);
    //
    if (blocks_.at(block_id).GetMaterialPointer()->NumStateVariables() > 0) {
      if (field_ids_.state_sym_tensor >= 0) {
        auto state_sym_tensor_n =
            GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.state_sym_tensor, nimble::STEP_N);
        auto state_sym_tensor_np1 =
            GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.state_sym_tensor, nimble::STEP_NP1);
        Kokkos::deep_copy(state_sym_tensor_n, state_sym_tensor_np1);
      }
      //
      if (field_ids_.state_full_tensor >= 0) {
        auto state_full_tensor_n =
            GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.state_full_tensor, nimble::STEP_N);
        auto state_full_tensor_np1 =
            GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.state_full_tensor, nimble::STEP_NP1);
        Kokkos::deep_copy(state_full_tensor_n, state_full_tensor_np1);
      }
      //
      if (field_ids_.state_scalar >= 0) {
        auto state_scalar_n =
            GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.state_scalar, nimble::STEP_N);
        auto state_scalar_np1 =
            GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.state_scalar, nimble::STEP_NP1);
        Kokkos::deep_copy(state_scalar_n, state_scalar_np1);
      }
      //
      if (field_ids_.state_vec3D >= 0) {
        auto state_vector_n = GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.state_vec3D, nimble::STEP_N);
        auto state_vector_np1 =
            GetDeviceFullTensorIntegrationPointData(block_id, field_ids_.state_vec3D, nimble::STEP_NP1);
        Kokkos::deep_copy(state_vector_n, state_vector_np1);
      }
    }
    //
    block_index += 1;
  }
}

nimble::Viewify<1>
ModelData::GetScalarNodeData(int field_id)
{
  auto field_view   = GetHostScalarNodeData(field_id);
  auto field_size   = static_cast<int>(field_view.extent(0));
  auto field_stride = static_cast<int>(field_view.stride_0());
  return {field_view.data(), {field_size}, {field_stride}};
}

nimble::Viewify<2>
ModelData::GetVectorNodeData(int field_id)
{
  auto field_view = GetHostVectorNodeData(field_id);
  auto size0      = static_cast<int>(field_view.extent(0));
  auto size1      = static_cast<int>(field_view.extent(1));
  auto stride0    = static_cast<int>(field_view.stride_0());
  auto stride1    = static_cast<int>(field_view.stride_1());
  return {field_view.data(), {size0, size1}, {stride0, stride1}};
}

void
ModelData::InitializeGatheredVectors(const nimble::GenesisMesh& mesh_)
{
  int num_blocks = static_cast<int>(mesh_.GetNumBlocks());

  gathered_reference_coordinate_d.resize(
      num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_reference_coordinates", 1));
  gathered_displacement_d.resize(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_displacement", 1));
  gathered_internal_force_d.resize(
      num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_internal_force", 1));
  gathered_contact_force_d.resize(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_contact_force", 1));

  int block_index = 0;
  for (const auto& block_it : blocks_) {
    int block_id          = block_it.first;
    int num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    Kokkos::resize(gathered_reference_coordinate_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_displacement_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_internal_force_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_contact_force_d.at(block_index), num_elem_in_block);
    block_index += 1;
  }
}

void
ModelData::ComputeLumpedMass(nimble::DataManager& data_manager)
{
  const auto& mesh_               = data_manager.GetMesh();
  const auto& parser_             = data_manager.GetParser();
  auto&       field_ids_          = data_manager.GetFieldIDs();
  auto        vector_communicator = data_manager.GetVectorCommunicator();

  int num_nodes  = static_cast<int>(mesh_.GetNumNodes());
  int num_blocks = static_cast<int>(mesh_.GetNumBlocks());

  std::vector<nimble_kokkos::DeviceScalarNodeGatheredView> gathered_lumped_mass_d(
      num_blocks, nimble_kokkos::DeviceScalarNodeGatheredView("gathered_lumped_mass", 1));
  int block_index = 0;
  for (const auto& block_it : blocks_) {
    int block_id          = block_it.first;
    int num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    Kokkos::resize(gathered_lumped_mass_d.at(block_index), num_elem_in_block);
    block_index += 1;
  }

  auto lumped_mass_h = GetHostScalarNodeData(field_ids_.lumped_mass);
  Kokkos::deep_copy(lumped_mass_h, (double)(0.0));

  auto lumped_mass_d = GetDeviceScalarNodeData(field_ids_.lumped_mass);
  Kokkos::deep_copy(lumped_mass_d, (double)(0.0));

  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto displacement         = GetVectorNodeData("displacement");

  critical_time_step_ = std::numeric_limits<double>::max();

  block_index = 0;
  for (auto& block_it : blocks_) {
    int                   block_id           = block_it.first;
    nimble_kokkos::Block& block              = block_it.second;
    nimble::Element*      element_d          = block.GetDeviceElement();
    double                density            = block.GetDensity();
    int                   num_elem_in_block  = mesh_.GetNumElementsInBlock(block_id);
    int                   num_nodes_per_elem = mesh_.GetNumNodesPerElement(block_id);
    int                   elem_conn_length   = num_elem_in_block * num_nodes_per_elem;
    int const*            elem_conn          = mesh_.GetConnectivity(block_id);

    nimble_kokkos::HostElementConnectivityView elem_conn_h("element_connectivity_h", elem_conn_length);
    for (int i = 0; i < elem_conn_length; i++) { elem_conn_h(i) = elem_conn[i]; }
    auto&& elem_conn_d = block.GetDeviceElementConnectivityView();
    Kokkos::resize(elem_conn_d, elem_conn_length);
    Kokkos::deep_copy(elem_conn_d, elem_conn_h);

    auto gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
    auto gathered_lumped_mass_block_d          = gathered_lumped_mass_d.at(block_index);

    GatherVectorNodeData(
        field_ids_.reference_coordinates,
        num_elem_in_block,
        num_nodes_per_elem,
        elem_conn_d,
        gathered_reference_coordinate_block_d);

    // COMPUTE LUMPED MASS
    Kokkos::parallel_for(
        "Lumped Mass", num_elem_in_block, KOKKOS_LAMBDA(const int i_elem) {
          auto element_reference_coordinate_d =
              Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          auto element_lumped_mass_d = Kokkos::subview(gathered_lumped_mass_block_d, i_elem, Kokkos::ALL);
          element_d->ComputeLumpedMass(density, element_reference_coordinate_d, element_lumped_mass_d);
        });

    // SCATTER TO NODE DATA
    ScatterScalarNodeData(
        field_ids_.lumped_mass, num_elem_in_block, num_nodes_per_elem, elem_conn_d, gathered_lumped_mass_block_d);

    double block_critical_time_step =
        block.ComputeCriticalTimeStep(reference_coordinate, displacement, num_elem_in_block, elem_conn);
    if (block_critical_time_step < critical_time_step_) { critical_time_step_ = block_critical_time_step; }

    block_index += 1;
  }
  Kokkos::deep_copy(lumped_mass_h, lumped_mass_d);

  // MPI vector reduction on lumped mass
  std::vector<double> mpi_scalar_buffer(num_nodes);
  for (unsigned int i = 0; i < num_nodes; i++) { mpi_scalar_buffer[i] = lumped_mass_h(i); }
  vector_communicator->VectorReduction(1, mpi_scalar_buffer.data());
  for (int i = 0; i < num_nodes; i++) { lumped_mass_h(i) = mpi_scalar_buffer[i]; }
}

void
ModelData::ComputeInternalForce(
    nimble::DataManager&      data_manager,
    double                    time_previous,
    double                    time_current,
    bool                      is_output_step,
    const nimble::Viewify<2>& displacement,
    nimble::Viewify<2>&       force)
{
  //
  // This version does not use the Viewify objects displacement and force
  // It may need to be updated for quasi-static simulations
  //

  const auto& mesh      = data_manager.GetMesh();
  const auto& field_ids = data_manager.GetFieldIDs();

  auto block_material_interface_factory = data_manager.GetBlockMaterialInterfaceFactory();

  nimble_kokkos::DeviceVectorNodeView internal_force_h = GetHostVectorNodeData(field_ids.internal_force);
  nimble_kokkos::DeviceVectorNodeView internal_force_d = GetDeviceVectorNodeData(field_ids.internal_force);
  Kokkos::deep_copy(internal_force_d, (double)(0.0));

  // Compute element-level kinematics
  constexpr int mpi_vector_dim = 3;

  int block_index = 0;
  for (auto& block_it : blocks_) {
    //
    const int block_id           = block_it.first;
    const int num_elem_in_block  = mesh.GetNumElementsInBlock(block_id);
    const int num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);

    nimble_kokkos::Block& block     = block_it.second;
    nimble::Element*      element_d = block.GetDeviceElement();

    auto elem_conn_d                           = block.GetDeviceElementConnectivityView();
    auto gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
    auto gathered_displacement_block_d         = gathered_displacement_d.at(block_index);
    auto gathered_internal_force_block_d       = gathered_internal_force_d.at(block_index);

    GatherVectorNodeData(
        field_ids.reference_coordinates,
        num_elem_in_block,
        num_nodes_per_elem,
        elem_conn_d,
        gathered_reference_coordinate_block_d);

    GatherVectorNodeData(
        field_ids.displacement, num_elem_in_block, num_nodes_per_elem, elem_conn_d, gathered_displacement_block_d);

    auto deformation_gradient_step_np1_d =
        GetDeviceFullTensorIntegrationPointData(block_id, field_ids.deformation_gradient, nimble::STEP_NP1);

    //
    // COMPUTE DEFORMATION GRADIENTS
    //
    Kokkos::parallel_for(
        "Deformation Gradient", num_elem_in_block, KOKKOS_LAMBDA(const int i_elem) {
          nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d =
              Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL(), Kokkos::ALL());
          nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d =
              Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL(), Kokkos::ALL());
          nimble_kokkos::DeviceFullTensorIntPtSubView element_deformation_gradient_step_np1_d =
              Kokkos::subview(deformation_gradient_step_np1_d, i_elem, Kokkos::ALL(), Kokkos::ALL());
          element_d->ComputeDeformationGradients(
              element_reference_coordinate_d, element_displacement_d, element_deformation_gradient_step_np1_d);
        });

    //
    // Insert the update for the state variables
    //

    block_index += 1;
  }

  if (block_material_interface_factory) {
    auto block_material_interface =
        block_material_interface_factory->create(time_previous, time_current, field_ids, block_data_, this);
    block_material_interface->ComputeStress();
  }

  //
  // Stress divergence
  //
  block_index = 0;
  for (auto& block_it : blocks_) {
    const int block_id           = block_it.first;
    const int num_elem_in_block  = mesh.GetNumElementsInBlock(block_id);
    const int num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);

    nimble_kokkos::Block& block     = block_it.second;
    nimble::Element*      element_d = block.GetDeviceElement();

    nimble_kokkos::DeviceElementConnectivityView elem_conn_d = block.GetDeviceElementConnectivityView();
    nimble_kokkos::DeviceVectorNodeGatheredView  gathered_reference_coordinate_block_d =
        gathered_reference_coordinate_d.at(block_index);
    nimble_kokkos::DeviceVectorNodeGatheredView gathered_displacement_block_d = gathered_displacement_d.at(block_index);
    nimble_kokkos::DeviceVectorNodeGatheredView gathered_internal_force_block_d =
        gathered_internal_force_d.at(block_index);

    nimble_kokkos::DeviceSymTensorIntPtView stress_step_np1_d =
        GetDeviceSymTensorIntegrationPointData(block_id, field_ids.stress, nimble::STEP_NP1);

    // COMPUTE NODAL FORCES
    Kokkos::parallel_for(
        "Force", num_elem_in_block, KOKKOS_LAMBDA(const int i_elem) {
          nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d =
              Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d =
              Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceSymTensorIntPtSubView element_stress_step_np1_d =
              Kokkos::subview(stress_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          nimble_kokkos::DeviceVectorNodeGatheredSubView element_internal_force_d =
              Kokkos::subview(gathered_internal_force_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
          element_d->ComputeNodalForces(
              element_reference_coordinate_d,
              element_displacement_d,
              element_stress_step_np1_d,
              element_internal_force_d);
        });

    ScatterVectorNodeData(
        field_ids.internal_force, num_elem_in_block, num_nodes_per_elem, elem_conn_d, gathered_internal_force_block_d);

    block_index += 1;
  }  // loop over blocks

  Kokkos::deep_copy(internal_force_h, internal_force_d);

  auto myVectorCommunicator = data_manager.GetVectorCommunicator();
  myVectorCommunicator->VectorReduction(mpi_vector_dim, internal_force_h);
}

void
ModelData::InitializeBlockData(nimble::DataManager& data_manager)
{
  //
  // Build up block data for stress computation
  //
  const auto& mesh = data_manager.GetMesh();
  for (auto&& block_it : blocks_) {
    int                   block_id                           = block_it.first;
    nimble_kokkos::Block& block                              = block_it.second;
    nimble::Material*     material_d                         = block.GetDeviceMaterialModel();
    int                   num_elem_in_block                  = mesh.GetNumElementsInBlock(block_id);
    int                   num_integration_points_per_element = block.GetHostElement()->NumIntegrationPointsPerElement();
    block_data_.emplace_back(&block, material_d, block_id, num_elem_in_block, num_integration_points_per_element);
  }
}

void
ModelData::InitializeExodusOutput(nimble::DataManager& data_manager)
{
  const auto& mesh_   = data_manager.GetMesh();
  const auto& parser_ = data_manager.GetParser();

  // Initialize the exodus-output-manager
  exodus_output_manager_->SpecifyOutputFields(*this, parser_.GetOutputFieldString());

  output_node_component_labels_    = std::move(exodus_output_manager_->GetNodeDataLabelsForOutput());
  output_element_component_labels_ = std::move(exodus_output_manager_->GetElementDataLabelsForOutput());

  derived_output_element_data_labels_.clear();
  std::vector<int> block_ids = mesh_.GetBlockIds();
  for (auto block_id : block_ids) {
    derived_output_element_data_labels_[block_id] = std::vector<std::string>();  // TODO eliminate this
  }

  auto& field_ids = data_manager.GetFieldIDs();
  displacement_h_ = GetHostVectorNodeData(field_ids.displacement);
  displacement_d_ = GetDeviceVectorNodeData(field_ids.displacement);

  velocity_h_ = GetHostVectorNodeData(field_ids.velocity);
  velocity_d_ = GetDeviceVectorNodeData(field_ids.velocity);
}

void
ModelData::WriteExodusOutput(nimble::DataManager& data_manager, double time_current)
{
  auto        mesh_         = data_manager.GetMesh();
  const auto& parser_       = data_manager.GetParser();
  auto        exodus_output = data_manager.GetExodusOutput();

  Kokkos::deep_copy(displacement_d_, displacement_h_);
  Kokkos::deep_copy(velocity_d_, velocity_h_);

  exodus_output_manager_->ComputeElementData(
      mesh_, (*this), blocks_, gathered_reference_coordinate_d, gathered_displacement_d);

  auto const& node_data_output = exodus_output_manager_->GetNodeDataForOutput(*this);
  auto const& elem_data_output = exodus_output_manager_->GetElementDataForOutput(*this);

  std::vector<double>                             glbl_data;
  std::map<int, std::vector<std::vector<double>>> drvd_elem_data;

  exodus_output->WriteStep(
      time_current,
      glbl_data,
      node_data_output,
      output_element_component_labels_,
      elem_data_output,
      derived_output_element_data_labels_,
      drvd_elem_data);
}

std::vector<std::string>
ModelData::GetScalarNodeDataLabels() const
{
  std::vector<std::string> node_data_labels;
  for (auto const& entry : field_label_to_field_id_map_) {
    std::string const& field_label = entry.first;
    int                field_id    = entry.second;
    for (auto const& node_entry : field_id_to_host_node_data_index_) {
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

std::vector<std::string>
ModelData::GetVectorNodeDataLabels() const
{
  std::vector<std::string> node_data_labels;
  for (auto const& entry : field_label_to_field_id_map_) {
    std::string const& field_label = entry.first;
    int                field_id    = entry.second;
    for (auto const& node_entry : field_id_to_host_node_data_index_) {
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

std::vector<std::string>
ModelData::GetSymmetricTensorIntegrationPointDataLabels(int block_id) const
{
  int                      block_index = block_id_to_integration_point_data_index_.at(block_id);
  int                      num_blocks  = static_cast<int>(block_id_to_integration_point_data_index_.size());
  std::vector<std::string> ipt_data_labels;
  for (auto const& entry : field_label_to_field_id_map_) {
    std::string const& field_label = entry.first;
    int                field_id    = entry.second;
    for (auto const& ipt_entry : field_id_to_device_integration_point_data_index_.at(block_index)) {
      int ipt_data_field_id = ipt_entry.first;
      if (field_id == ipt_data_field_id) {
        int ipt_data_index = ipt_entry.second;
        if (device_integration_point_data_step_np1_.at(block_index).at(ipt_data_index)->type() ==
            FieldType::DeviceSymTensorIntPt) {
          if (std::find(ipt_data_labels.begin(), ipt_data_labels.end(), field_label) == ipt_data_labels.end()) {
            ipt_data_labels.push_back(field_label);
          }
        }
      }
    }
  }
  return ipt_data_labels;
}

std::vector<std::string>
ModelData::GetFullTensorIntegrationPointDataLabels(int block_id) const
{
  int                      block_index = block_id_to_integration_point_data_index_.at(block_id);
  std::vector<std::string> ipt_data_labels;
  for (auto const& entry : field_label_to_field_id_map_) {
    std::string const& field_label = entry.first;
    int                field_id    = entry.second;
    for (auto const& ipt_entry : field_id_to_device_integration_point_data_index_.at(block_index)) {
      int ipt_data_field_id = ipt_entry.first;
      if (field_id == ipt_data_field_id) {
        int ipt_data_index = ipt_entry.second;
        if (device_integration_point_data_step_np1_.at(block_index).at(ipt_data_index)->type() ==
            FieldType::DeviceFullTensorIntPt) {
          if (std::find(ipt_data_labels.begin(), ipt_data_labels.end(), field_label) == ipt_data_labels.end()) {
            ipt_data_labels.push_back(field_label);
          }
        }
      }
    }
  }
  return ipt_data_labels;
}

HostScalarNodeView
ModelData::GetHostScalarNodeData(int field_id)
{
  int        index             = field_id_to_host_node_data_index_.at(field_id);
  FieldBase* base_field_ptr    = host_node_data_.at(index).get();
  auto       derived_field_ptr = dynamic_cast<Field<FieldType::HostScalarNode>*>(base_field_ptr);
  return derived_field_ptr->data();
}

HostVectorNodeView
ModelData::GetHostVectorNodeData(int field_id)
{
  int        index             = field_id_to_host_node_data_index_.at(field_id);
  FieldBase* base_field_ptr    = host_node_data_.at(index).get();
  auto       derived_field_ptr = dynamic_cast<Field<FieldType::HostVectorNode>*>(base_field_ptr);
  return derived_field_ptr->data();
}

HostScalarElemView
ModelData::GetHostScalarElementData(int block_id, int field_id)
{
  int        block_index    = block_id_to_element_data_index_.at(block_id);
  int        data_index     = field_id_to_host_element_data_index_.at(block_index).at(field_id);
  FieldBase* base_field_ptr = nullptr;
  base_field_ptr            = host_element_data_.at(block_index).at(data_index).get();
  auto derived_field_ptr    = dynamic_cast<Field<FieldType::HostScalarElem>*>(base_field_ptr);
  return derived_field_ptr->data();
}

HostSymTensorIntPtView
ModelData::GetHostSymTensorIntegrationPointData(int block_id, int field_id, nimble::Step step)
{
  int        block_index    = block_id_to_integration_point_data_index_.at(block_id);
  int        data_index     = field_id_to_host_integration_point_data_index_.at(block_index).at(field_id);
  FieldBase* base_field_ptr = nullptr;
  if (step == nimble::STEP_N) {
    base_field_ptr = host_integration_point_data_step_n_.at(block_index).at(data_index).get();
  } else if (step == nimble::STEP_NP1) {
    base_field_ptr = host_integration_point_data_step_np1_.at(block_index).at(data_index).get();
  }
  auto derived_field_ptr = dynamic_cast<Field<FieldType::HostSymTensorIntPt>*>(base_field_ptr);
  return derived_field_ptr->data();
}

HostFullTensorIntPtView
ModelData::GetHostFullTensorIntegrationPointData(int block_id, int field_id, nimble::Step step)
{
  int        block_index    = block_id_to_integration_point_data_index_.at(block_id);
  int        data_index     = field_id_to_host_integration_point_data_index_.at(block_index).at(field_id);
  FieldBase* base_field_ptr = nullptr;
  if (step == nimble::STEP_N) {
    base_field_ptr = host_integration_point_data_step_n_.at(block_index).at(data_index).get();
  } else if (step == nimble::STEP_NP1) {
    base_field_ptr = host_integration_point_data_step_np1_.at(block_index).at(data_index).get();
  }
  auto derived_field_ptr = dynamic_cast<Field<FieldType::HostFullTensorIntPt>*>(base_field_ptr);
  return derived_field_ptr->data();
}

HostSymTensorElemView
ModelData::GetHostSymTensorElementData(int block_id, int field_id)
{
  int        block_index    = block_id_to_element_data_index_.at(block_id);
  int        data_index     = field_id_to_host_element_data_index_.at(block_index).at(field_id);
  FieldBase* base_field_ptr = nullptr;
  base_field_ptr            = host_element_data_.at(block_index).at(data_index).get();
  auto derived_field_ptr    = dynamic_cast<Field<FieldType::HostSymTensorElem>*>(base_field_ptr);
  return derived_field_ptr->data();
}

HostFullTensorElemView
ModelData::GetHostFullTensorElementData(int block_id, int field_id)
{
  int        block_index    = block_id_to_element_data_index_.at(block_id);
  int        data_index     = field_id_to_host_element_data_index_.at(block_index).at(field_id);
  FieldBase* base_field_ptr = nullptr;
  base_field_ptr            = host_element_data_.at(block_index).at(data_index).get();
  auto derived_field_ptr    = dynamic_cast<Field<FieldType::HostFullTensorElem>*>(base_field_ptr);
  return derived_field_ptr->data();
}

template <FieldType ft>
typename Field<ft>::View
ModelData::GetDeviceElementData(int block_id, int field_id)
{
  int  block_index       = block_id_to_element_data_index_.at(block_id);
  int  data_index        = field_id_to_device_element_data_index_.at(block_index).at(field_id);
  auto base_field_ptr    = device_element_data_.at(block_index).at(data_index).get();
  auto derived_field_ptr = dynamic_cast<Field<ft>*>(base_field_ptr);
  return derived_field_ptr->data();
}

template <FieldType ft>
typename Field<ft>::View
ModelData::GetDeviceIntPointData(int block_id, int field_id, nimble::Step step)
{
  int        block_index    = block_id_to_integration_point_data_index_.at(block_id);
  int        data_index     = field_id_to_device_integration_point_data_index_.at(block_index).at(field_id);
  FieldBase* base_field_ptr = nullptr;
  if (step == nimble::STEP_N) {
    base_field_ptr = device_integration_point_data_step_n_.at(block_index).at(data_index).get();
  } else if (step == nimble::STEP_NP1) {
    base_field_ptr = device_integration_point_data_step_np1_.at(block_index).at(data_index).get();
  }
  auto derived_field = dynamic_cast<Field<ft>&>(*base_field_ptr);
  return derived_field.data();
}

DeviceScalarNodeView
ModelData::GetDeviceScalarNodeData(int field_id)
{
  int        index             = field_id_to_device_node_data_index_.at(field_id);
  FieldBase* base_field_ptr    = device_node_data_.at(index).get();
  auto       derived_field_ptr = dynamic_cast<Field<FieldType::DeviceScalarNode>*>(base_field_ptr);
  return derived_field_ptr->data();
}

DeviceVectorNodeView
ModelData::GetDeviceVectorNodeData(int field_id)
{
  int        index             = field_id_to_device_node_data_index_.at(field_id);
  FieldBase* base_field_ptr    = device_node_data_.at(index).get();
  auto       derived_field_ptr = dynamic_cast<Field<FieldType::DeviceVectorNode>*>(base_field_ptr);
  return derived_field_ptr->data();
}

DeviceSymTensorIntPtView
ModelData::GetDeviceSymTensorIntegrationPointData(int block_id, int field_id, nimble::Step step)
{
  return GetDeviceIntPointData<FieldType::DeviceSymTensorIntPt>(block_id, field_id, step);
}

DeviceFullTensorIntPtView
ModelData::GetDeviceFullTensorIntegrationPointData(int block_id, int field_id, nimble::Step step)
{
  return GetDeviceIntPointData<FieldType::DeviceFullTensorIntPt>(block_id, field_id, step);
}

DeviceScalarIntPtView
ModelData::GetDeviceScalarIntegrationPointData(int block_id, int field_id, nimble::Step step)
{
  return GetDeviceIntPointData<FieldType::DeviceScalarIntPt>(block_id, field_id, step);
}

DeviceVectorIntPtView
ModelData::GetDeviceVectorIntegrationPointData(int block_id, int field_id, nimble::Step step)
{
  return GetDeviceIntPointData<FieldType::DeviceVectorIntPt>(block_id, field_id, step);
}

DeviceScalarElemView
ModelData::GetDeviceScalarElementData(int block_id, int field_id)
{
  return GetDeviceElementData<FieldType::DeviceScalarElem>(block_id, field_id);
}

DeviceSymTensorElemView
ModelData::GetDeviceSymTensorElementData(int block_id, int field_id)
{
  return GetDeviceElementData<FieldType::DeviceSymTensorElem>(block_id, field_id);
}

DeviceFullTensorElemView
ModelData::GetDeviceFullTensorElementData(int block_id, int field_id)
{
  return GetDeviceElementData<FieldType::DeviceFullTensorElem>(block_id, field_id);
}

DeviceScalarNodeGatheredView
ModelData::GatherScalarNodeData(
    int                                  field_id,
    int                                  num_elements,
    int                                  num_nodes_per_element,
    const DeviceElementConnectivityView& elem_conn_d,
    DeviceScalarNodeGatheredView         gathered_view_d)
{
  int        index             = field_id_to_device_node_data_index_.at(field_id);
  FieldBase* base_field_ptr    = device_node_data_.at(index).get();
  auto       derived_field_ptr = dynamic_cast<Field<FieldType::DeviceScalarNode>*>(base_field_ptr);
  auto       data              = derived_field_ptr->data();
  Kokkos::parallel_for(
      "GatherScalarNodeData", num_elements, KOKKOS_LAMBDA(const int i_elem) {
        for (int i_node = 0; i_node < num_nodes_per_element; i_node++) {
          gathered_view_d(i_elem, i_node) = data(elem_conn_d(num_nodes_per_element * i_elem + i_node));
        }
      });
  return gathered_view_d;
}

DeviceVectorNodeGatheredView
ModelData::GatherVectorNodeData(
    int                                  field_id,
    int                                  num_elements,
    int                                  num_nodes_per_element,
    const DeviceElementConnectivityView& elem_conn_d,
    DeviceVectorNodeGatheredView         gathered_view_d)
{
  int        index             = field_id_to_device_node_data_index_.at(field_id);
  FieldBase* base_field_ptr    = device_node_data_.at(index).get();
  auto       derived_field_ptr = dynamic_cast<Field<FieldType::DeviceVectorNode>*>(base_field_ptr);
  auto       data              = derived_field_ptr->data();
  Kokkos::parallel_for(
      "GatherVectorNodeData", num_elements, KOKKOS_LAMBDA(const int i_elem) {
        for (int i_node = 0; i_node < num_nodes_per_element; i_node++) {
          int node_index = elem_conn_d(num_nodes_per_element * i_elem + i_node);
          for (int i_coord = 0; i_coord < 3; i_coord++) {
            gathered_view_d(i_elem, i_node, i_coord) = data(node_index, i_coord);
          }
        }
      });
  return gathered_view_d;
}

void
ModelData::ScatterScalarNodeData(
    int                                  field_id,
    int                                  num_elements,
    int                                  num_nodes_per_element,
    const DeviceElementConnectivityView& elem_conn_d,
    const DeviceScalarNodeGatheredView&  gathered_view_d)
{
  int        index             = field_id_to_device_node_data_index_.at(field_id);
  FieldBase* base_field_ptr    = device_node_data_.at(index).get();
  auto       derived_field_ptr = dynamic_cast<Field<FieldType::DeviceScalarNode>*>(base_field_ptr);
  Field<FieldType::DeviceScalarNode>::AtomicView data = derived_field_ptr->data();
  Kokkos::parallel_for(
      "ScatterScalarNodeData", num_elements, KOKKOS_LAMBDA(const int i_elem) {
        for (int i_node = 0; i_node < num_nodes_per_element; i_node++) {
          data(elem_conn_d(num_nodes_per_element * i_elem + i_node)) += gathered_view_d(i_elem, i_node);
        }
      });
}

void
ModelData::ScatterVectorNodeData(
    int                                  field_id,
    int                                  num_elements,
    int                                  num_nodes_per_element,
    const DeviceElementConnectivityView& elem_conn_d,
    const DeviceVectorNodeGatheredView&  gathered_view_d)
{
  int        index             = field_id_to_device_node_data_index_.at(field_id);
  FieldBase* base_field_ptr    = device_node_data_.at(index).get();
  auto       derived_field_ptr = dynamic_cast<Field<FieldType::DeviceVectorNode>*>(base_field_ptr);
  Field<FieldType::DeviceVectorNode>::AtomicView data = derived_field_ptr->data();
  Kokkos::parallel_for(
      "ScatterVectorNodeData", num_elements, KOKKOS_LAMBDA(const int i_elem) {
        for (int i_node = 0; i_node < num_nodes_per_element; i_node++) {
          int node_index = elem_conn_d(num_nodes_per_element * i_elem + i_node);
          for (int i_coord = 0; i_coord < 3; i_coord++) {
            data(node_index, i_coord) += gathered_view_d(i_elem, i_node, i_coord);
          }
        }
      });
}

#ifndef KOKKOS_ENABLE_QTHREADS
void
ModelData::ScatterScalarNodeDataUsingKokkosScatterView(
    int                                  field_id,
    int                                  num_elements,
    int                                  num_nodes_per_element,
    const DeviceElementConnectivityView& elem_conn_d,
    const DeviceScalarNodeGatheredView&  gathered_view_d)
{
  int        index             = field_id_to_device_node_data_index_.at(field_id);
  FieldBase* base_field_ptr    = device_node_data_.at(index).get();
  auto       derived_field_ptr = dynamic_cast<Field<FieldType::DeviceScalarNode>*>(base_field_ptr);
  auto       data              = derived_field_ptr->data();
  auto       scatter_view =
      Kokkos::Experimental::create_scatter_view(data);  // DJL it is a terrible idea to allocate this here
  scatter_view.reset();
  Kokkos::parallel_for(
      "GatherVectorNodeData", num_elements, KOKKOS_LAMBDA(const int i_elem) {
        auto scattered_access = scatter_view.access();
        for (int i_node = 0; i_node < num_nodes_per_element; i_node++) {
          scattered_access(elem_conn_d(num_nodes_per_element * i_elem + i_node)) += gathered_view_d(i_elem, i_node);
        }
      });
  Kokkos::Experimental::contribute(data, scatter_view);
}
#endif

void
ModelData::UpdateWithNewVelocity(nimble::DataManager& data_manager, double dt)
{
  Kokkos::deep_copy(velocity_d_, velocity_h_);
}

void
ModelData::UpdateWithNewDisplacement(nimble::DataManager& data_manager, double dt)
{
  Kokkos::deep_copy(displacement_d_, displacement_h_);
}

}  // namespace nimble_kokkos
