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

#include <stdexcept>

#include "nimble_data_manager.h"
#include "nimble_kokkos_model_data.h"
#include "nimble_kokkos_material_factory.h"
#include "nimble_vector_communicator.h"

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


void ModelData::InitializeBlocks(
    nimble::DataManager &data_manager,
    const std::shared_ptr<MaterialFactoryType> &material_factory_base)
{

  bool store_unrotated_stress(true);

  const auto& mesh_ = data_manager.GetMesh();
  const auto& parser_ = data_manager.GetParser();
  auto &field_ids_ = data_manager.GetFieldIDs();

  auto material_factory_ptr = dynamic_cast< nimble_kokkos::MaterialFactory* >(material_factory_base.get());
  const auto num_blocks = static_cast<int>(mesh_.GetNumBlocks());

  //
  // Blocks
  //
  std::vector<int> block_ids = mesh_.GetBlockIds();
  for (int i=0 ; i<num_blocks ; i++){
    int block_id = block_ids.at(i);
    std::string const & macro_material_parameters = parser_.GetMacroscaleMaterialParameters(block_id);
    std::map<int, std::string> const & rve_material_parameters = parser_.GetMicroscaleMaterialParameters();
    std::string rve_bc_strategy = parser_.GetMicroscaleBoundaryConditionStrategy();
    int num_elements_in_block = mesh_.GetNumElementsInBlock(block_id);
    blocks_[block_id] = nimble_kokkos::Block();
    blocks_.at(block_id).Initialize(macro_material_parameters, num_elements_in_block, *material_factory_ptr);
    //
    // MPI version use model_data.DeclareElementData(block_id, data_labels_and_lengths);
    //
    std::vector<double> initial_value(9, 0.0);
    initial_value[0] = initial_value[1] = initial_value[2] = 1.0;
    field_ids_.deformation_gradient = AllocateIntegrationPointData(block_id, nimble::FULL_TENSOR,
                                                                   "deformation_gradient",
                                                                   num_elements_in_block, initial_value);
    // volume-averaged quantities for I/O are stored as element data
    AllocateElementData(block_id, nimble::FULL_TENSOR, "deformation_gradient", num_elements_in_block);

    field_ids_.stress = AllocateIntegrationPointData(block_id, nimble::SYMMETRIC_TENSOR, "stress",
                                                     num_elements_in_block);
    if (store_unrotated_stress) {
      field_ids_.unrotated_stress = AllocateIntegrationPointData(block_id, nimble::SYMMETRIC_TENSOR, "stress",
                                                                 num_elements_in_block);
    }

    // volume-averaged quantities for I/O are stored as element data
    AllocateElementData(block_id, nimble::SYMMETRIC_TENSOR, "stress", num_elements_in_block);

    if (parser_.GetOutputFieldString().find("volume") != std::string::npos) {
      AllocateElementData(block_id, nimble::SCALAR, "volume", num_elements_in_block);
    }
  }

  // Initialize the exodus-output-manager
  auto &exodus_output_manager = GetExodusOutputManager();
  exodus_output_manager.SpecifyOutputFields(*this, parser_.GetOutputFieldString());

  output_node_component_labels_ = std::move(exodus_output_manager.GetNodeDataLabelsForOutput());
  output_element_component_labels_ = std::move(exodus_output_manager.GetElementDataLabelsForOutput());

  derived_output_element_data_labels_.clear();
  for (auto block_id : block_ids) {
    derived_output_element_data_labels_[block_id] = std::vector<std::string>(); // TODO eliminate this
  }

  // Initialize gathered containers when using explicit scheme
  if (parser_.TimeIntegrationScheme() == "explicit")
    InitializeGatheredVectors(mesh_);

}


void ModelData::SwapStates(const nimble::DataManager &data_manager)
{
  const auto &field_ids_ = data_manager.GetFieldIDs();

  // Copy STEP_NP1 data to STEP_N
  int block_index = 0;
  for (auto &block_it : blocks_) {
    int block_id = block_it.first;
    auto deformation_gradient_step_n_d = GetDeviceFullTensorIntegrationPointData(
        block_id, field_ids_.deformation_gradient, nimble::STEP_N);
    auto unrotated_stress_step_n_d = GetDeviceSymTensorIntegrationPointData(block_id,
                                                                            field_ids_.unrotated_stress,
                                                                            nimble::STEP_N);
    auto stress_step_n_d = GetDeviceSymTensorIntegrationPointData(block_id, field_ids_.stress,
                                                                  nimble::STEP_N);
    auto deformation_gradient_step_np1_d = GetDeviceFullTensorIntegrationPointData(
        block_id, field_ids_.deformation_gradient, nimble::STEP_NP1);
    auto unrotated_stress_step_np1_d = GetDeviceSymTensorIntegrationPointData(block_id,
                                                                              field_ids_.unrotated_stress,
                                                                              nimble::STEP_NP1);
    auto stress_step_np1_d = GetDeviceSymTensorIntegrationPointData(block_id, field_ids_.stress,
                                                                    nimble::STEP_NP1);
    Kokkos::deep_copy(deformation_gradient_step_n_d, deformation_gradient_step_np1_d);
    Kokkos::deep_copy(unrotated_stress_step_n_d, unrotated_stress_step_np1_d);
    Kokkos::deep_copy(stress_step_n_d, stress_step_np1_d);
    block_index += 1;
  }
}


Viewify ModelData::GetScalarNodeData(const std::string& label)
{
  auto field_id = GetFieldId(label);
  if (field_id < 0) {
    std::string code = " Field " + label + " Not Allocated ";
    throw std::runtime_error(code);
  }
  return {GetHostScalarNodeData(field_id).data(), 1};
}


void ModelData::InitializeGatheredVectors(const nimble::GenesisMesh &mesh_)
{
  int num_blocks = static_cast<int>(mesh_.GetNumBlocks());

  gathered_reference_coordinate_d.resize(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_reference_coordinates", 1));
  gathered_displacement_d.resize(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_displacement", 1));
  gathered_internal_force_d.resize(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_internal_force", 1));
  gathered_contact_force_d.resize(num_blocks, nimble_kokkos::DeviceVectorNodeGatheredView("gathered_contact_force", 1));


  int block_index = 0;
  for (const auto &block_it : blocks_) {
    int block_id = block_it.first;
    int num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    Kokkos::resize(gathered_reference_coordinate_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_displacement_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_internal_force_d.at(block_index), num_elem_in_block);
    Kokkos::resize(gathered_contact_force_d.at(block_index), num_elem_in_block);
    block_index += 1;
  }

}


void ModelData::ComputeLumpedMass(nimble::DataManager &data_manager)
{

  const auto& mesh_ = data_manager.GetMesh();
  const auto& parser_ = data_manager.GetParser();
  auto &field_ids_ = data_manager.GetFieldIDs();
  auto vector_communicator = data_manager.GetVectorCommunicator();

  int num_nodes = static_cast<int>(mesh_.GetNumNodes());
  int num_blocks = static_cast<int>(mesh_.GetNumBlocks());

  std::vector<nimble_kokkos::DeviceScalarNodeGatheredView> gathered_lumped_mass_d(num_blocks, nimble_kokkos::DeviceScalarNodeGatheredView("gathered_lumped_mass", 1));
  int block_index = 0;
  for (const auto &block_it : blocks_) {
    int block_id = block_it.first;
    int num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    Kokkos::resize(gathered_lumped_mass_d.at(block_index), num_elem_in_block);
    block_index += 1;
  }

  auto lumped_mass_h = GetHostScalarNodeData(field_ids_.lumped_mass);
  Kokkos::deep_copy(lumped_mass_h, (double)(0.0));

  auto lumped_mass_d = GetDeviceScalarNodeData(field_ids_.lumped_mass);
  Kokkos::deep_copy(lumped_mass_d, (double)(0.0));

  block_index = 0;
  for (auto &block_it : blocks_) {
    int block_id = block_it.first;
    nimble_kokkos::Block& block = block_it.second;
    nimble::Element* element_d = block.GetDeviceElement();
    double density = block.GetDensity();
    int num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    int num_nodes_per_elem = mesh_.GetNumNodesPerElement(block_id);
    int elem_conn_length = num_elem_in_block * num_nodes_per_elem;
    int const * elem_conn = mesh_.GetConnectivity(block_id);

    nimble_kokkos::HostElementConnectivityView elem_conn_h("element_connectivity_h", elem_conn_length);
    for (int i=0 ; i<elem_conn_length ; i++) {
      elem_conn_h(i) = elem_conn[i];
    }
    auto &&elem_conn_d = block.GetDeviceElementConnectivityView();
    Kokkos::resize(elem_conn_d, elem_conn_length);
    Kokkos::deep_copy(elem_conn_d, elem_conn_h);

    auto gathered_reference_coordinate_block_d = gathered_reference_coordinate_d.at(block_index);
    auto gathered_lumped_mass_block_d = gathered_lumped_mass_d.at(block_index);

    GatherVectorNodeData(field_ids_.reference_coordinates,
                         num_elem_in_block,
                         num_nodes_per_elem,
                         elem_conn_d,
                         gathered_reference_coordinate_block_d);

    // COMPUTE LUMPED MASS
    Kokkos::parallel_for("Lumped Mass", num_elem_in_block, KOKKOS_LAMBDA (const int i_elem) {
      auto element_reference_coordinate_d = Kokkos::subview(gathered_reference_coordinate_block_d, i_elem, Kokkos::ALL, Kokkos::ALL);
      auto element_lumped_mass_d = Kokkos::subview(gathered_lumped_mass_block_d, i_elem, Kokkos::ALL);
      element_d->ComputeLumpedMass(density, element_reference_coordinate_d, element_lumped_mass_d);
    });

    // SCATTER TO NODE DATA
    ScatterScalarNodeData(field_ids_.lumped_mass, num_elem_in_block,
                          num_nodes_per_elem, elem_conn_d,
                          gathered_lumped_mass_block_d);

    block_index += 1;
  }
  Kokkos::deep_copy(lumped_mass_h, lumped_mass_d);

  // MPI vector reduction on lumped mass
  std::vector<double> mpi_scalar_buffer(num_nodes);
  for (unsigned int i=0 ; i<num_nodes ; i++) {
    mpi_scalar_buffer[i] = lumped_mass_h(i);
  }
  vector_communicator->VectorReduction(1, mpi_scalar_buffer.data());
  for (int i=0 ; i<num_nodes ; i++) {
    lumped_mass_h(i) = mpi_scalar_buffer[i];
  }

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
  Field< FieldType::HostVectorNode >* derived_field_ptr = dynamic_cast< Field< FieldType::HostVectorNode >* >(base_field_ptr);
  return derived_field_ptr->data();
}

HostScalarElemView ModelData::GetHostScalarElementData(int block_id,
                                                       int field_id) {
  int block_index = block_id_to_element_data_index_.at(block_id);
  int data_index = field_id_to_host_element_data_index_.at(block_index).at(field_id);
  FieldBase* base_field_ptr(0);
  base_field_ptr = host_element_data_.at(block_index).at(data_index).get();
  Field< FieldType::HostScalarElem >* derived_field_ptr = dynamic_cast< Field< FieldType::HostScalarElem >* >(base_field_ptr);
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
  Field< FieldType::HostSymTensorIntPt >* derived_field_ptr = dynamic_cast< Field< FieldType::HostSymTensorIntPt >* >(base_field_ptr);
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
  FieldBase *base_field_ptr = device_node_data_.at(index).get();
  auto derived_field_ptr =
      dynamic_cast<Field<FieldType::DeviceVectorNode> *>(base_field_ptr);
  auto data = derived_field_ptr->data();
  Kokkos::parallel_for("GatherVectorNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
    for (int i_node = 0; i_node < num_nodes_per_element; i_node++) {
      int node_index = elem_conn_d(num_nodes_per_element * i_elem + i_node);
      for (int i_coord = 0; i_coord < 3; i_coord++) {
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
  Kokkos::parallel_for("ScatterScalarNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
    for (int i_node=0 ; i_node < num_nodes_per_element ; i_node++) {
      data(elem_conn_d(num_nodes_per_element*i_elem + i_node)) += gathered_view_d(i_elem, i_node);
    }
  });
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
  Kokkos::parallel_for("ScatterVectorNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
    for (int i_node=0 ; i_node < num_nodes_per_element ; i_node++) {
      int node_index = elem_conn_d(num_nodes_per_element*i_elem + i_node);
      for (int i_coord=0 ; i_coord < 3 ; i_coord++) {
        data(node_index, i_coord) += gathered_view_d(i_elem, i_node, i_coord);
      }
    }
  });
}

void ModelData::SetReferenceCoordinates(const nimble::GenesisMesh &mesh)
{
  const double * const ref_coord_x = mesh.GetCoordinatesX();
  const double * const ref_coord_y = mesh.GetCoordinatesY();
  const double * const ref_coord_z = mesh.GetCoordinatesZ();
  auto field_id = GetFieldId("reference_coordinate");
  nimble_kokkos::HostVectorNodeView reference_coordinate_h = GetHostVectorNodeData(field_id);
  nimble_kokkos::DeviceVectorNodeView reference_coordinate_d = GetDeviceVectorNodeData(field_id);
  int num_nodes = static_cast<int>(mesh.GetNumNodes());
  for (int i=0 ; i<num_nodes ; i++) {
    reference_coordinate_h(i, 0) = ref_coord_x[i];
    reference_coordinate_h(i, 1) = ref_coord_y[i];
    reference_coordinate_h(i, 2) = ref_coord_z[i];
  }
  Kokkos::deep_copy(reference_coordinate_d, reference_coordinate_h);
}

#ifndef KOKKOS_ENABLE_QTHREADS
void ModelData::ScatterScalarNodeDataUsingKokkosScatterView(int field_id,
                                                            int num_elements,
                                                            int num_nodes_per_element,
                                                            DeviceElementConnectivityView elem_conn_d,
                                                            DeviceScalarNodeGatheredView gathered_view_d) {
  int index = field_id_to_device_node_data_index_.at(field_id);
  FieldBase* base_field_ptr = device_node_data_.at(index).get();
  Field< FieldType::DeviceScalarNode >* derived_field_ptr = static_cast< Field< FieldType::DeviceScalarNode >* >(base_field_ptr);
  auto data = derived_field_ptr->data();
  auto scatter_view = Kokkos::Experimental::create_scatter_view(data); // DJL it is a terrible idea to allocate this here
  scatter_view.reset();
  Kokkos::parallel_for("GatherVectorNodeData", num_elements, KOKKOS_LAMBDA (const int i_elem) {
    auto scattered_access = scatter_view.access();
    for (int i_node=0 ; i_node < num_nodes_per_element ; i_node++) {
      scattered_access(elem_conn_d(num_nodes_per_element*i_elem + i_node)) += gathered_view_d(i_elem, i_node);
    }
  });
  Kokkos::Experimental::contribute(data, scatter_view);
}
#endif


} // namespace nimble
