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

#include "nimble_model_data.h"

#include <algorithm>
#include <cctype>
#include <sstream>

#include "nimble_data_manager.h"
#include "nimble_genesis_mesh.h"
#include "nimble_macros.h"
#include "nimble_material_factory.h"
#include "nimble_parser.h"
#include "nimble_vector_communicator.h"
#include "nimble_view.h"

#ifdef NIMBLE_HAVE_UQ
#include "uq/nimble_uq_block.h"
#endif

namespace nimble {

int
ModelData::GetFieldId(const std::string& label) const
{
  for (auto const& id_field_pair : data_fields_) {
    if (id_field_pair.second.label_ == label) { return id_field_pair.first; }
  }
  return -1;
}

Field
ModelData::GetField(int field_id)
{
  return data_fields_.at(field_id);
}

int
ModelData::AllocateNodeData(Length length, std::string label, int num_objects)
{
  Field field;
  field.label_    = label;
  field.relation_ = NODE;
  field.length_   = length;
  AssignFieldId(field);
  int array_size = num_objects;
  array_size *= LengthToInt(length, dim_);
  node_data_[field.id_] = std::vector<double>(array_size);
  return field.id_;
}

double*
ModelData::GetNodeData(int field_id)
{
  return node_data_.at(field_id).data();
}

void
ModelData::GetNodeDataForOutput(std::vector<std::vector<double>>& single_component_arrays)
{
  unsigned int number_of_output_fields = output_node_component_labels_.size();
  if (single_component_arrays.size() != number_of_output_fields) {
    single_component_arrays.resize(number_of_output_fields);
  }

  for (unsigned int i_output_label = 0; i_output_label < output_node_component_labels_.size(); i_output_label++) {
    std::string output_label = output_node_component_labels_[i_output_label];
    for (std::map<int, std::vector<double>>::const_iterator it = node_data_.begin(); it != node_data_.end(); it++) {
      Field                    field            = data_fields_.at(it->first);
      std::vector<std::string> component_labels = GetComponentLabels(field.label_, field.length_, dim_);
      for (unsigned int i_component = 0; i_component < component_labels.size(); i_component++) {
        if (component_labels[i_component] == output_label) {
          int array_size = it->second.size() / LengthToInt(field.length_, dim_);
          if (single_component_arrays[i_output_label].size() != array_size) {
            single_component_arrays[i_output_label].resize(array_size);
          }
          GetNodeDataComponent(field.id_, i_component, &single_component_arrays[i_output_label][0]);
        }
      }
    }
  }
}

void
ModelData::DeclareElementData(int block_id, std::vector<std::pair<std::string, Length>> const& data_labels_and_lengths)
{
  if (element_data_fields_.find(block_id) == element_data_fields_.end()) {
    element_data_fields_[block_id] = std::vector<int>();
  }

  for (std::pair<std::string, Length> const& entry : data_labels_and_lengths) {
    Field field;
    field.label_    = entry.first;
    field.relation_ = ELEMENT;
    field.length_   = entry.second;
    AssignFieldId(field);
    element_data_fields_[block_id].push_back(field.id_);
  }
}

void
ModelData::AllocateElementData(std::map<int, int> const& num_integration_points_in_each_block)
{
  for (std::map<int, std::vector<int>>::const_iterator it = element_data_fields_.begin();
       it != element_data_fields_.end();
       it++) {
    int                      block_id                        = it->first;
    std::vector<int> const&  field_ids                       = it->second;
    int                      array_length_for_each_int_point = 0;
    std::vector<std::string> component_labels_for_block;

    for (const int& field_id : field_ids) {
      Field field = GetField(field_id);
      array_length_for_each_int_point += LengthToInt(field.length_, dim_);
      std::vector<std::string> component_labels = GetComponentLabels(field.label_, field.length_, dim_);
      component_labels_for_block.insert(
          component_labels_for_block.end(), component_labels.begin(), component_labels.end());
    }

    block_ids_.push_back(block_id);
    element_component_labels_[block_id]           = component_labels_for_block;
    output_element_component_labels_[block_id]    = std::vector<std::string>();
    derived_output_element_data_labels_[block_id] = std::vector<std::string>();
    int array_length            = array_length_for_each_int_point * num_integration_points_in_each_block.at(block_id);
    element_data_n_[block_id]   = std::vector<double>(array_length);
    element_data_np1_[block_id] = std::vector<double>(array_length);
  }
}

void
ModelData::GetElementDataForOutput(std::map<int, std::vector<std::vector<double>>>& single_component_arrays)
{
  // Allocate space in single_component_arrays, if necessary
  for (auto& output_element_component_label : output_element_component_labels_) {
    int block_id = output_element_component_label.first;
    if (single_component_arrays.find(block_id) == single_component_arrays.end()) {
      single_component_arrays[block_id] = std::vector<std::vector<double>>();
    }
    unsigned int num_output_fields_for_block = output_element_component_label.second.size();
    if (single_component_arrays[block_id].size() != num_output_fields_for_block) {
      single_component_arrays[block_id].resize(num_output_fields_for_block);
    }
  }

  // Copy the data, in component form, to single_component_arrays
  for (auto& output_element_component_label : output_element_component_labels_) {
    int                             block_id = output_element_component_label.first;
    std::vector<std::string> const& output_element_component_labels_for_this_block =
        output_element_component_label.second;
    unsigned int num_element_variables = element_component_labels_.at(block_id).size();
    unsigned int num_output_variables  = output_element_component_labels_for_this_block.size();
    unsigned int array_size            = element_data_np1_.at(block_id).size() / num_element_variables;
    for (unsigned int output_array_index = 0; output_array_index < num_output_variables; output_array_index++) {
      std::string          output_label = output_element_component_labels_for_this_block[output_array_index];
      std::vector<double>& output_array = single_component_arrays[block_id][output_array_index];
      if (output_array.size() != array_size) { output_array.resize(array_size); }
      // determine the offset into the element data array
      int offset = -1;
      for (unsigned int i = 0; i < element_component_labels_[block_id].size(); i++) {
        if (output_label == element_component_labels_[block_id][i]) { offset = i; }
      }
      if (offset == -1) {
        NIMBLE_ABORT(
            "\n**** Error in ModelData::GetElementDataForOutput(), output "
            "label not found.\n");
      }
      for (int i = 0; i < array_size; i++) {
        output_array[i] = element_data_np1_[block_id][i * num_element_variables + offset];
      }
    }
  }
}

void
ModelData::SpecifyOutputFields(const std::string& output_field_string)
{
  // Parse the string into individual field labels
  std::vector<std::string> requested_output_labels;
  std::stringstream        ss(output_field_string);
  std::string              entry;
  while (std::getline(ss, entry, ' ')) { requested_output_labels.push_back(entry); }

  std::vector<std::string> node_noncomponent_labels;
  std::vector<std::string> node_component_labels;
  std::vector<std::string> element_noncomponent_labels;
  std::vector<std::string> element_component_labels;
  std::vector<std::string> element_int_pt_noncomponent_labels;
  std::vector<std::string> element_int_pt_component_labels;
  for (std::map<int, Field>::const_iterator field_it = data_fields_.begin(); field_it != data_fields_.end();
       field_it++) {
    std::string label    = field_it->second.label_;
    Length      length   = field_it->second.length_;
    Relation    relation = field_it->second.relation_;
    if (relation == NODE) {
      node_noncomponent_labels.push_back(label);
      std::vector<std::string> component_labels = GetComponentLabels(label, length, dim_);
      for (const auto& component_label : component_labels) { node_component_labels.push_back(component_label); }
    } else if (relation == ELEMENT) {
      if (!HasIntegrationPointPrefix(label)) {
        NIMBLE_ABORT(
            "\n**** Error, ModelData::SpecifyOutputFields() expected "
            "integration point prefix on data label \"" +
            label + "\".\n");
      }
      element_int_pt_noncomponent_labels.push_back(label);
      std::vector<std::string> int_pt_component_labels = GetComponentLabels(label, length, dim_);
      for (const auto& int_pt_component_label : int_pt_component_labels) {
        element_int_pt_component_labels.push_back(int_pt_component_label);
      }
      std::string label_without_prefix = RemoveIntegrationPointPrefix(label);
      if (std::find(element_noncomponent_labels.begin(), element_noncomponent_labels.end(), label_without_prefix) ==
          element_noncomponent_labels.end()) {
        element_noncomponent_labels.push_back(label_without_prefix);
        std::vector<std::string> component_labels = GetComponentLabels(label_without_prefix, length, dim_);
        for (const auto& component_label : component_labels) { element_component_labels.push_back(component_label); }
      }
    }
  }

  // Record the field labels for output
  for (auto const& requested_label : requested_output_labels) {
    // Requested label is node component data
    if (std::find(node_component_labels.begin(), node_component_labels.end(), requested_label) !=
        node_component_labels.end()) {
      output_node_component_labels_.push_back(requested_label);
    }

    // Requested label is a vector or tensor at a node (output all the
    // vector/tensor components)
    else if (
        std::find(node_noncomponent_labels.begin(), node_noncomponent_labels.end(), requested_label) !=
        node_noncomponent_labels.end()) {
      Length                   length           = LabelToLength(requested_label, data_fields_, dim_);
      std::vector<std::string> component_labels = GetComponentLabels(requested_label, length, dim_);
      for (auto const& label : component_labels) { output_node_component_labels_.push_back(label); }
    }

    // Requested label is a single component of a vector or tensor at the
    // integration points that will be volume averaged over the element
    else if (
        std::find(element_component_labels.begin(), element_component_labels.end(), requested_label) !=
        element_component_labels.end()) {
      for (auto const& block_id_label_pair : element_component_labels_) {
        int  block_id              = block_id_label_pair.first;
        bool label_exists_on_block = false;
        for (auto const& label_on_block : block_id_label_pair.second) {
          if (requested_label == RemoveIntegrationPointPrefix(label_on_block)) { label_exists_on_block = true; }
        }
        if (label_exists_on_block) { derived_output_element_data_labels_.at(block_id).push_back(requested_label); }
      }
    }

    // Requested label matches a component label at a specific integration point
    else if (
        std::find(element_int_pt_component_labels.begin(), element_int_pt_component_labels.end(), requested_label) !=
        element_int_pt_component_labels.end()) {
      for (auto const& block_id_label_pair : element_component_labels_) {
        int  block_id              = block_id_label_pair.first;
        bool label_exists_on_block = false;
        for (auto const& label_on_block : block_id_label_pair.second) {
          if (requested_label == label_on_block) { label_exists_on_block = true; }
        }
        if (label_exists_on_block) { output_element_component_labels_.at(block_id).push_back(requested_label); }
      }
    }

    // Requested label is for a scalar, vector, or tensor at the integration
    // points that will be volume averaged over the element (output all the
    // vector/tensor components)
    else if (
        std::find(element_noncomponent_labels.begin(), element_noncomponent_labels.end(), requested_label) !=
        element_noncomponent_labels.end()) {
      Length                   length           = LabelToLength(requested_label, data_fields_, dim_);
      std::vector<std::string> component_labels = GetComponentLabels(requested_label, length, dim_);
      for (auto const& label : component_labels) {
        for (auto const& it : element_component_labels_) {
          int                             block_id                                     = it.first;
          std::vector<std::string> const& available_element_component_labels_for_block = it.second;
          for (auto const& available_label : available_element_component_labels_for_block) {
            std::string label_without_int_pt_prefix = RemoveIntegrationPointPrefix(available_label);
            if (label == label_without_int_pt_prefix) {
              if (std::find(
                      derived_output_element_data_labels_.at(block_id).begin(),
                      derived_output_element_data_labels_.at(block_id).end(),
                      label) == derived_output_element_data_labels_.at(block_id).end()) {
                derived_output_element_data_labels_.at(block_id).push_back(label);
              }
            }
          }
        }
      }
    }

    // Requested label is for a scalar, vector, or tensor at a specific
    // integration point (output all the vector/tensor components)
    else if (
        std::find(
            element_int_pt_noncomponent_labels.begin(), element_int_pt_noncomponent_labels.end(), requested_label) !=
        element_int_pt_noncomponent_labels.end()) {
      Length                   length           = LabelToLength(requested_label, data_fields_, dim_);
      std::vector<std::string> component_labels = GetComponentLabels(requested_label, length, dim_);
      for (auto const& label : component_labels) {
        for (auto const& it : element_component_labels_) {
          int                             block_id                                     = it.first;
          std::vector<std::string> const& available_element_component_labels_for_block = it.second;
          for (auto const& available_label : available_element_component_labels_for_block) {
            if (label == available_label) { output_element_component_labels_.at(block_id).push_back(label); }
          }
        }
      }
    }

    // Special cases
    else if (requested_label == "volume") {
      for (auto const& block_id : block_ids_) {
        derived_output_element_data_labels_.at(block_id).push_back(requested_label);
      }
    }

    else {
      NIMBLE_ABORT(
          "\nError:  ModelData::SpecifyOutputFields(), unable to process "
          "requested output \"" +
          requested_label + "\".\n");
    }
  }
}

void
ModelData::AssignFieldId(Field& field)
{
  int         field_id = -1;
  std::string name(field.label_);
  std::transform(field.label_.begin(), field.label_.end(), name.begin(), ::tolower);
  for (std::map<int, Field>::const_iterator it = data_fields_.begin(); it != data_fields_.end(); it++) {
    std::string temp_name(it->second.label_);
    std::transform(it->second.label_.begin(), it->second.label_.end(), temp_name.begin(), ::tolower);
    if (name == temp_name) {
      if (field.relation_ == it->second.relation_ && field.length_ == it->second.length_) {
        field_id = it->second.id_;
        break;
      } else {
        NIMBLE_ABORT(
            "\nError:  ModelData::AssignFieldId(), Inconsistent fields cannot "
            "be assigned the same label.  Label = " +
            field.label_ + "\n");
      }
    }
  }
  if (field_id == -1) {
    field_id               = data_fields_.size();
    field.id_              = field_id;
    data_fields_[field_id] = field;
  } else {
    field.id_ = field_id;
  }
}

void
ModelData::GetNodeDataComponent(int field_id, int component, double* const component_data)
{
  Field&               field          = data_fields_.at(field_id);
  int                  num_components = LengthToInt(field.length_, dim_);
  std::vector<double>& data           = node_data_.at(field_id);

  if (component >= num_components) {
    NIMBLE_ABORT("\nError:  Invalid component in ModelData::GetNodeDataComponent\n");
  }

  for (unsigned int i = 0; i < data.size() / num_components; i++) {
    component_data[i] = data[i * num_components + component];
  }
}

template <typename TBlock>
void
ModelData::EmplaceBlocks(
    nimble::DataManager&                                data_manager,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base)
{
  const auto& mesh_   = data_manager.GetMesh();
  const auto& parser_ = data_manager.GetParser();

  auto material_factory_ptr = dynamic_cast<nimble::MaterialFactory*>(material_factory_base.get());

  const std::vector<int>& block_ids = mesh_.GetBlockIds();
  for (int block_id : block_ids) {
    std::string const& model_material_parameters = parser_.GetModelMaterialParameters(block_id);
    auto               block_ptr                 = std::shared_ptr<nimble::Block>(new TBlock());
    block_ptr->Initialize(model_material_parameters, *material_factory_ptr);
    std::vector<std::pair<std::string, nimble::Length>> data_labels_and_lengths;
    block_ptr->GetDataLabelsAndLengths(data_labels_and_lengths);
    DeclareElementData(block_id, data_labels_and_lengths);
    blocks_[block_id] = block_ptr;
  }
}

template void
ModelData::EmplaceBlocks<nimble::Block>(
    nimble::DataManager&                                data_manager,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base);

#ifdef NIMBLE_HAVE_UQ
template void
ModelData::EmplaceBlocks<::nimble_uq::Block>(
    nimble::DataManager&                                data_manager,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base);
#endif

void
ModelData::AllocateInitializeElementData(
    nimble::DataManager&                                data_manager,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base)
{
  const auto& mesh_   = data_manager.GetMesh();
  const auto& parser_ = data_manager.GetParser();

  auto material_factory_ptr = dynamic_cast<nimble::MaterialFactory*>(material_factory_base.get());

  std::map<int, int> num_elem_in_each_block = mesh_.GetNumElementsInBlock();
  AllocateElementData(num_elem_in_each_block);
  SpecifyOutputFields(parser_.GetOutputFieldString());
  std::map<int, std::vector<std::string>> const& elem_data_labels         = GetElementDataLabels();
  std::map<int, std::vector<std::string>> const& derived_elem_data_labels = GetDerivedElementDataLabelsForOutput();

  // Initialize the element data
  for (auto& block_it : blocks_) {
    int                     block_id          = block_it.first;
    int                     num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    std::vector<int> const& elem_global_ids   = mesh_.GetElementGlobalIdsInBlock(block_id);
    auto&                   block             = block_it.second;
    std::vector<double>&    elem_data_n       = GetElementDataOld(block_id);
    std::vector<double>&    elem_data_np1     = GetElementDataNew(block_id);
    block->InitializeElementData(
        num_elem_in_block,
        elem_global_ids,
        elem_data_labels.at(block_id),
        derived_elem_data_labels.at(block_id),
        elem_data_n,
        elem_data_np1,
        *material_factory_ptr,
        data_manager);
  }
}
void
ModelData::InitializeBlocks(
    nimble::DataManager&                                data_manager,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base)
{
  EmplaceBlocks<nimble::Block>(data_manager, material_factory_base);
  AllocateInitializeElementData(data_manager, material_factory_base);
}

void
ModelData::ComputeLumpedMass(nimble::DataManager& data_manager)
{
  const auto& mesh_   = data_manager.GetMesh();
  const auto& parser_ = data_manager.GetParser();

  auto lumped_mass_field_id = GetFieldId("lumped_mass");

  auto lumped_mass          = GetNodeData(lumped_mass_field_id);
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto displacement         = GetVectorNodeData("displacement");

  //
  // Computed the lumped mass matrix (diagonal matrix) and the critical time
  // step
  //
  critical_time_step_ = std::numeric_limits<double>::max();
  for (auto& block_it : blocks_) {
    int        block_id          = block_it.first;
    int        num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    int const* elem_conn         = mesh_.GetConnectivity(block_id);
    auto&      block             = block_it.second;
    block->ComputeLumpedMassMatrix(reference_coordinate.data(), num_elem_in_block, elem_conn, lumped_mass);
    double block_critical_time_step =
        block->ComputeCriticalTimeStep(reference_coordinate, displacement, num_elem_in_block, elem_conn);
    if (block_critical_time_step < critical_time_step_) { critical_time_step_ = block_critical_time_step; }
  }

  //
  // Perform a vector reduction on lumped mass.  This is a scalar nodal
  // quantity.
  //
  auto          vector_comm      = data_manager.GetVectorCommunicator();
  constexpr int scalar_dimension = 1;
  vector_comm->VectorReduction(scalar_dimension, lumped_mass);
}

nimble::Viewify<1>
ModelData::GetScalarNodeData(int field_id)
{
  auto& field_vector = node_data_.at(field_id);
  auto  isize        = static_cast<int>(field_vector.size());
  return {field_vector.data(), {isize}, {1}};
}

nimble::Viewify<2>
ModelData::GetVectorNodeData(int field_id)
{
  auto& field_vector = node_data_.at(field_id);
  auto  isize        = static_cast<int>(field_vector.size()) / 3;
  return {field_vector.data(), {isize, 3}, {3, 1}};
}

void
ModelData::InitializeExodusOutput(nimble::DataManager& data_manager)
{
  for (auto& block_it : blocks_) {
    int block_id                 = block_it.first;
    derived_elem_data_[block_id] = std::vector<std::vector<double>>();
  }
}

void
ModelData::WriteExodusOutput(nimble::DataManager& data_manager, double time_current)
{
  const auto& mesh_         = data_manager.GetMesh();
  auto        exodus_output = data_manager.GetExodusOutput();

  std::vector<double> global_data;

  auto reference_coord_ = GetNodeData("reference_coordinate");
  auto displacement     = GetNodeData("displacement");

  for (auto& block_it : blocks_) {
    int         block_id          = block_it.first;
    auto        block             = block_it.second;
    int         num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    int const*  elem_conn         = mesh_.GetConnectivity(block_id);
    const auto& elem_data_np1     = GetElementDataNew(block_id);
    block->ComputeDerivedElementData(
        reference_coord_,
        displacement,
        num_elem_in_block,
        elem_conn,
        element_component_labels_.at(block_id).size(),
        elem_data_np1,
        derived_output_element_data_labels_.at(block_id).size(),
        derived_elem_data_.at(block_id));
  }

  GetNodeDataForOutput(node_data_for_output_);
  GetElementDataForOutput(elem_data_for_output_);

  exodus_output->WriteStep(
      time_current,
      global_data,
      node_data_for_output_,
      output_element_component_labels_,
      elem_data_for_output_,
      derived_output_element_data_labels_,
      derived_elem_data_);
}

double*
ModelData::GetNodeData(const std::string& label)
{
  const auto field_id = GetFieldId(label);
  if (field_id < 0) {
    std::string code = " Field " + label + " Not Allocated ";
    NIMBLE_ABORT(code);
  }
  return node_data_.at(field_id).data();
}

void
ModelData::ComputeExternalForce(
    nimble::DataManager& data_manager,
    double               time_previous,
    double               time_current,
    bool                 is_output_step)
{
  auto force = GetVectorNodeData("external_force");
  force.zero();
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
  const auto& mesh = data_manager.GetMesh();

  force.zero();

  auto reference_coord = GetNodeData("reference_coordinate");
  auto velocity        = GetNodeData("velocity");

  for (auto& block_it : blocks_) {
    int                        block_id          = block_it.first;
    int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    int const*                 elem_conn         = mesh.GetConnectivity(block_id);
    std::vector<int> const&    elem_global_ids   = mesh.GetElementGlobalIdsInBlock(block_id);
    auto&                      block             = block_it.second;
    std::vector<double> const& elem_data_n       = GetElementDataOld(block_id);
    std::vector<double>&       elem_data_np1     = GetElementDataNew(block_id);
    block->ComputeInternalForce(
        reference_coord,
        displacement.data(),
        velocity,
        force.data(),
        time_previous,
        time_current,
        num_elem_in_block,
        elem_conn,
        elem_global_ids.data(),
        element_component_labels_.at(block_id),
        elem_data_n,
        elem_data_np1,
        data_manager,
        is_output_step);
  }

  // DJL
  // Perform a vector reduction on internal force.  This is a vector nodal
  // quantity.
  auto          vector_comm      = data_manager.GetVectorCommunicator();
  constexpr int vector_dimension = 3;
  vector_comm->VectorReduction(vector_dimension, force.data());
}

}  // namespace nimble
