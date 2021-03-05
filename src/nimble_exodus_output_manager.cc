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

#include <sstream>

#include "nimble_data_manager.h"
#include "nimble_exodus_output_manager.h"
#include "nimble_kokkos_model_data.h"
#include "nimble_utils.h"

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
          node_data_types_.push_back(FieldType::HostScalarNode);
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
    for (auto const & block_id : block_ids) {
      elem_data_labels_[block_id] = std::vector<std::string>();
      elem_data_iptdata_field_ids_[block_id] = std::vector<int>();
      elem_data_edata_field_ids_[block_id] = std::vector<int>();
      elem_data_types_[block_id] = std::vector<FieldType>();
      elem_data_integration_point_index_[block_id] = std::vector<int>();
      elem_data_components_[block_id] = std::vector<int>();
      elem_data_[block_id] = std::vector< std::vector<double> >();
      sym_tensor_field_ids_requiring_volume_average_[block_id] = std::vector<int>();
      full_tensor_field_ids_requiring_volume_average_[block_id] = std::vector<int>();
    }

    for (unsigned int i_block=0 ; i_block<block_ids.size() ; ++i_block) {
      int block_id = block_ids[i_block];

      // volume is special, it is not tied to integration point data
      for (auto const & requested_label : requested_labels) {
        if (requested_label == "volume") {
          output_element_volume_ = true;
          volume_field_id_ = model_data.GetFieldId(requested_label);
          int num_elem = model_data.GetDeviceScalarElementData(block_id, volume_field_id_).extent(0);
          elem_data_labels_[block_id].push_back(requested_label);
          elem_data_iptdata_field_ids_[block_id].push_back(volume_field_id_);
          elem_data_edata_field_ids_[block_id].push_back(volume_field_id_);
          elem_data_types_[block_id].push_back(FieldType::HostScalarElem);
          elem_data_integration_point_index_[block_id].push_back(-1);
          elem_data_components_[block_id].push_back(0);
          elem_data_[block_id].push_back(std::vector<double>(num_elem, 0.0));
        }
      }

      std::vector<std::string> symmetric_tensor_integration_point_data_labels = model_data.GetSymmetricTensorIntegrationPointDataLabels(block_id);
      for (auto const & requested_label : requested_labels) {
        for (auto& ipt_label : symmetric_tensor_integration_point_data_labels) {
          std::string field_label = requested_label;
          bool single_integration_point = false;
          bool require_volume_average = true;
          int integration_point_index = -1;
          if(requested_label.substr(0,3) == "ipt") {
            field_label = requested_label.substr(6, requested_label.size() - 6);
            single_integration_point = true;
            require_volume_average = false;
            integration_point_index = atoi(requested_label.substr(4,2).c_str()) - 1;
          }
          if (field_label == ipt_label) {
            int iptdata_field_id = model_data.GetFieldId(ipt_label);
            int edata_field_id = iptdata_field_id;
            if (require_volume_average) {
              sym_tensor_field_ids_requiring_volume_average_[block_id].push_back(iptdata_field_id);
            }
            int num_elem = model_data.GetDeviceSymTensorIntegrationPointData(block_id, iptdata_field_id, nimble::STEP_NP1).extent(0);
            if (single_integration_point) {
              edata_field_id = model_data.AllocateElementData(block_id, nimble::SYMMETRIC_TENSOR, requested_label, num_elem);
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

      std::vector<std::string> full_tensor_integration_point_data_labels = model_data.GetFullTensorIntegrationPointDataLabels(block_id);
      for (auto const & requested_label : requested_labels) {
        // Handle case where user has requested output for a specific integration point
        std::string field_label = requested_label;
        bool single_integration_point = false;
        bool require_volume_average = true;
        int integration_point_index = -1;
        if(requested_label.substr(0,3) == "ipt") {
          field_label = requested_label.substr(6, requested_label.size() - 6);
          single_integration_point = true;
          require_volume_average = false;
          integration_point_index = atoi(requested_label.substr(4,2).c_str()) - 1;
        }
        for (auto& ipt_label : full_tensor_integration_point_data_labels) {
          if (field_label == ipt_label) {
            int iptdata_field_id = model_data.GetFieldId(ipt_label);
            int edata_field_id = iptdata_field_id;
            if (require_volume_average) {
              full_tensor_field_ids_requiring_volume_average_[block_id].push_back(iptdata_field_id);
            }
            int num_elem = model_data.GetDeviceFullTensorIntegrationPointData(block_id, iptdata_field_id, nimble::STEP_NP1).extent(0);
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

  void ExodusOutputManager::ComputeElementData(nimble::GenesisMesh& mesh,
                                               nimble_kokkos::ModelData& model_data,
                                               std::map<int, nimble_kokkos::Block>& blocks,
                                               std::vector<nimble_kokkos::DeviceVectorNodeGatheredView>& gathered_reference_coordinate_d,
                                               std::vector<nimble_kokkos::DeviceVectorNodeGatheredView>& gathered_displacement_d) {
    int block_index;
    std::map<int, nimble_kokkos::Block>::iterator block_it;

    // Element volume
    if (output_element_volume_) {
      for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
        int block_id = block_it->first;
        nimble_kokkos::Block& block = block_it->second;
        nimble::Element* element_d = block.GetDeviceElement();
        int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
        int num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);
        nimble_kokkos::DeviceElementConnectivityView elem_conn_d = block.GetDeviceElementConnectivityView();
        nimble_kokkos::DeviceVectorNodeGatheredView gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
        nimble_kokkos::DeviceVectorNodeGatheredView gathered_displacement_block_d = gathered_displacement_d.at(block_index);
        nimble_kokkos::DeviceScalarElemView volume_d = model_data.GetDeviceScalarElementData(block_id, volume_field_id_);
        Kokkos::parallel_for("Element Volume", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d = Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceScalarElemSingleEntryView element_volume_d = Kokkos::subview(volume_d, i_elem);
            element_d->ComputeVolume(element_reference_coordinate_d,
                                     element_displacement_d,
                                     element_volume_d);
          });
        nimble_kokkos::HostScalarElemView volume_h = model_data.GetHostScalarElementData(block_id, volume_field_id_);
        deep_copy(volume_h, volume_d);
      }
    }

    // Extract data for specific integration points
    for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
      int block_id = block_it->first;
      int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      for (unsigned int i_data=0 ; i_data < elem_data_labels_.at(block_id).size() ; ++i_data) {
        int integration_point_index = elem_data_integration_point_index_.at(block_id).at(i_data);
        if (integration_point_index != -1) {
          int iptdata_field_id = elem_data_iptdata_field_ids_.at(block_id).at(i_data);
          int edata_field_id = elem_data_edata_field_ids_.at(block_id).at(i_data);
          FieldType field_type = elem_data_types_.at(block_id).at(i_data);
          if (field_type == FieldType::HostSymTensorElem){
            nimble_kokkos::DeviceSymTensorIntPtView sym_tensor_data_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id, iptdata_field_id, nimble::STEP_NP1);
            nimble_kokkos::DeviceSymTensorElemView sym_tensor_data_single_int_pt_d = model_data.GetDeviceSymTensorElementData(block_id, edata_field_id);
            Kokkos::parallel_for("Extract Symmetric Tensor Integration Point Data for Output", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
                nimble_kokkos::DeviceSymTensorIntPtSubView element_sym_tensor_step_np1_d = Kokkos::subview(sym_tensor_data_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
                nimble_kokkos::DeviceSymTensorElemSingleEntryView element_sym_tensor_single_int_pt_d = Kokkos::subview(sym_tensor_data_single_int_pt_d, i_elem, Kokkos::ALL);
                for (int i=0 ; i<element_sym_tensor_single_int_pt_d.extent(0) ; i++) {
                  element_sym_tensor_single_int_pt_d(i) = element_sym_tensor_step_np1_d(integration_point_index,i);
                }
              });
            nimble_kokkos::HostSymTensorElemView sym_tensor_data_single_int_pt_h = model_data.GetHostSymTensorElementData(block_id, edata_field_id);
            deep_copy(sym_tensor_data_single_int_pt_h, sym_tensor_data_single_int_pt_d);

          }
          else if (field_type == FieldType::HostFullTensorElem) {
            nimble_kokkos::DeviceFullTensorIntPtView full_tensor_data_step_np1_d = model_data.GetDeviceFullTensorIntegrationPointData(block_id, iptdata_field_id, nimble::STEP_NP1);
            nimble_kokkos::DeviceFullTensorElemView full_tensor_data_single_int_pt_d = model_data.GetDeviceFullTensorElementData(block_id, edata_field_id);
            Kokkos::parallel_for("Extract Full Tensor Integration Point Data for Output", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
                nimble_kokkos::DeviceFullTensorIntPtSubView element_full_tensor_step_np1_d = Kokkos::subview(full_tensor_data_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
                nimble_kokkos::DeviceFullTensorElemSingleEntryView element_full_tensor_single_int_pt_d = Kokkos::subview(full_tensor_data_single_int_pt_d, i_elem, Kokkos::ALL);
                for (int i=0 ; i<element_full_tensor_single_int_pt_d.extent(0) ; i++) {
                  element_full_tensor_single_int_pt_d(i) = element_full_tensor_step_np1_d(integration_point_index,i);
                }
              });
            nimble_kokkos::HostFullTensorElemView full_tensor_data_single_int_pt_h = model_data.GetHostFullTensorElementData(block_id, edata_field_id);
            deep_copy(full_tensor_data_single_int_pt_h, full_tensor_data_single_int_pt_d);
          }
        }
      }
    }


    // Volume averaging of symmetric and full tensors stored at integration points
    for (block_index=0, block_it=blocks.begin(); block_it!=blocks.end() ; block_index++, block_it++) {
      int block_id = block_it->first;
      nimble_kokkos::Block& block = block_it->second;
      nimble::Element* element_d = block.GetDeviceElement();
      int num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
      int num_nodes_per_elem = mesh.GetNumNodesPerElement(block_id);
      nimble_kokkos::DeviceElementConnectivityView elem_conn_d = block.GetDeviceElementConnectivityView();
      nimble_kokkos::DeviceVectorNodeGatheredView gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
      nimble_kokkos::DeviceVectorNodeGatheredView gathered_displacement_block_d = gathered_displacement_d.at(block_index);
      for(auto & field_id : sym_tensor_field_ids_requiring_volume_average_.at(block_id)) {
        nimble_kokkos::DeviceSymTensorIntPtView sym_tensor_data_step_np1_d = model_data.GetDeviceSymTensorIntegrationPointData(block_id, field_id, nimble::STEP_NP1);
        nimble_kokkos::DeviceSymTensorElemView sym_tensor_data_vol_ave_d = model_data.GetDeviceSymTensorElementData(block_id, field_id);
        Kokkos::parallel_for("Volume Averaging Sym Tensor", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d = Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceSymTensorIntPtSubView element_sym_tensor_step_np1_d = Kokkos::subview(sym_tensor_data_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceSymTensorElemSingleEntryView element_sym_tensor_vol_ave_d = Kokkos::subview(sym_tensor_data_vol_ave_d, i_elem, Kokkos::ALL);
            element_d->ComputeVolumeAverageSymTensor(element_reference_coordinate_d,
                                                     element_displacement_d,
                                                     element_sym_tensor_step_np1_d,
                                                     element_sym_tensor_vol_ave_d);
          });
        nimble_kokkos::HostSymTensorElemView sym_tensor_data_vol_ave_h = model_data.GetHostSymTensorElementData(block_id, field_id);
        deep_copy(sym_tensor_data_vol_ave_h, sym_tensor_data_vol_ave_d);
      }
      for(auto & field_id : full_tensor_field_ids_requiring_volume_average_.at(block_id)) {
        nimble_kokkos::DeviceFullTensorIntPtView full_tensor_data_step_np1_d = model_data.GetDeviceFullTensorIntegrationPointData(block_id, field_id, nimble::STEP_NP1);
        nimble_kokkos::DeviceFullTensorElemView full_tensor_data_vol_ave_d = model_data.GetDeviceFullTensorElementData(block_id, field_id);
        Kokkos::parallel_for("Volume Averaging Full Tensor", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceVectorNodeGatheredSubView element_displacement_d = Kokkos::subview(gathered_displacement_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceFullTensorIntPtSubView element_full_tensor_step_np1_d = Kokkos::subview(full_tensor_data_step_np1_d, i_elem, Kokkos::ALL, Kokkos::ALL);
            nimble_kokkos::DeviceFullTensorElemSingleEntryView element_full_tensor_vol_ave_d = Kokkos::subview(full_tensor_data_vol_ave_d, i_elem, Kokkos::ALL);
            element_d->ComputeVolumeAverageFullTensor(element_reference_coordinate_d,
                                                      element_displacement_d,
                                                      element_full_tensor_step_np1_d,
                                                      element_full_tensor_vol_ave_d);
          });
        nimble_kokkos::HostFullTensorElemView full_tensor_data_vol_ave_h = model_data.GetHostFullTensorElementData(block_id, field_id);
        deep_copy(full_tensor_data_vol_ave_h, full_tensor_data_vol_ave_d);
      }
    }
  }

  std::vector< std::vector<double> > ExodusOutputManager::GetNodeDataForOutput(nimble_kokkos::ModelData& model_data) {
    for (unsigned int i_data=0 ; i_data < node_data_labels_.size() ; ++i_data) {
      int field_id = node_data_field_ids_.at(i_data);
      FieldType field_type = node_data_types_.at(i_data);
      if (field_type == FieldType::HostScalarNode) {
        HostScalarNodeView data = model_data.GetHostScalarNodeData(field_id);
        for (unsigned int i=0 ; i<node_data_[i_data].size() ; i++) {
          node_data_[i_data][i] = data(i);
        }
      }
      else if (field_type == FieldType::HostVectorNode) {
        HostVectorNodeView data = model_data.GetHostVectorNodeData(field_id);
        int component = node_data_components_[i_data];
        for (unsigned int i=0 ; i<node_data_[i_data].size() ; i++) {
          node_data_[i_data][i] = data(i, component);
        }
      }
    }
    return node_data_;
  }

  std::map<int, std::vector< std::vector<double> > > ExodusOutputManager::GetElementDataForOutput(nimble_kokkos::ModelData& model_data) {
    std::vector<int> block_ids = model_data.GetBlockIds();
    for (auto const & block_id : block_ids) {
      for (unsigned int i_data=0 ; i_data < elem_data_labels_.at(block_id).size() ; ++i_data) {
        int edata_field_id = elem_data_edata_field_ids_.at(block_id).at(i_data);
        FieldType field_type = elem_data_types_.at(block_id).at(i_data);
        if (field_type == FieldType::HostScalarElem) {
          HostScalarElemView data = model_data.GetHostScalarElementData(block_id, edata_field_id);
          for (unsigned int i=0 ; i<elem_data_.at(block_id)[i_data].size() ; i++) {
            elem_data_.at(block_id)[i_data][i] = data(i);
          }
        }
        else if (field_type == FieldType::HostSymTensorElem) {
          HostSymTensorElemView data = model_data.GetHostSymTensorElementData(block_id, edata_field_id);
          int component = elem_data_components_.at(block_id).at(i_data);
          for (unsigned int i=0 ; i<elem_data_.at(block_id)[i_data].size() ; i++) {
            elem_data_.at(block_id)[i_data][i] = data(i, component);
          }
        }
        else if (field_type == FieldType::HostFullTensorElem) {
          HostFullTensorElemView data = model_data.GetHostFullTensorElementData(block_id, edata_field_id);
          int component = elem_data_components_.at(block_id).at(i_data);
          for (unsigned int i=0 ; i<elem_data_.at(block_id)[i_data].size() ; i++) {
            elem_data_.at(block_id)[i_data][i] = data(i, component);
          }
        }
      }
    }
    return elem_data_;
  }

}
