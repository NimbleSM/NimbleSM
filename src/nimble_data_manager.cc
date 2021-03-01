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

#include "nimble_data_manager.h"

#ifdef NIMBLE_HAVE_KOKKOS
#include "nimble_kokkos_model_data.h"
#include "nimble_kokkos_material_factory.h"
#endif

#include "nimble_material_factory.h"
#include "nimble_model_data.h"

#include <algorithm>
#include <stdexcept>

namespace nimble {


DataManager::DataManager(const nimble::Parser &parser,
                         const nimble::GenesisMesh &mesh,
                         const nimble::GenesisMesh &rve_mesh)
    : parser_(parser), mesh_(mesh), rve_mesh_(rve_mesh),
      macroscale_data_(), rve_data_()
{
#ifdef NIMBLE_HAVE_KOKKOS
  if (parser_.UseKokkos()) {
    macroscale_data_ = std::make_shared<nimble_kokkos::ModelData>();
    return;
  }
#endif
  macroscale_data_ = std::make_shared<nimble::ModelData>();
}


void DataManager::Initialize(const std::shared_ptr<MaterialFactoryType> &material_factory_base)
{

  const auto dim = static_cast<int>(mesh_.GetDim());
  const auto num_nodes = static_cast<int>(mesh_.GetNumNodes());

  macroscale_data_->SetDimension(dim);

  field_ids_.lumped_mass = macroscale_data_->AllocateNodeData(nimble::SCALAR, "lumped_mass", num_nodes);
  field_ids_.reference_coordinates = macroscale_data_->AllocateNodeData(nimble::VECTOR, "reference_coordinate", num_nodes);
  field_ids_.displacement = macroscale_data_->AllocateNodeData(nimble::VECTOR, "displacement", num_nodes);
  field_ids_.velocity = macroscale_data_->AllocateNodeData(nimble::VECTOR, "velocity", num_nodes);
  field_ids_.acceleration = macroscale_data_->AllocateNodeData(nimble::VECTOR, "acceleration", num_nodes);
  field_ids_.internal_force = macroscale_data_->AllocateNodeData(nimble::VECTOR, "internal_force", num_nodes);
  field_ids_.contact_force = macroscale_data_->AllocateNodeData(nimble::VECTOR, "contact_force", num_nodes);

  if (!parser_.UseKokkos()) {
    macroscale_data_->AllocateNodeData(nimble::VECTOR, "trial_displacement", num_nodes);
    macroscale_data_->AllocateNodeData(nimble::VECTOR, "displacement_fluctuation", num_nodes);
    macroscale_data_->AllocateNodeData(nimble::VECTOR, "trial_internal_force", num_nodes);
    macroscale_data_->AllocateNodeData(nimble::VECTOR, "external_force", num_nodes);
    macroscale_data_->AllocateNodeData(nimble::SCALAR, "skin_node", num_nodes);
  }

  macroscale_data_->SetReferenceCoordinates(mesh_);
  //
  std::shared_ptr<nimble::MaterialFactoryBase> material_factory;
  if (parser_.UseKokkos()) {
#ifdef NIMBLE_HAVE_KOKKOS
    material_factory = std::shared_ptr<nimble::MaterialFactoryBase>(new nimble_kokkos::MaterialFactory);
#else
    throw std::runtime_error(" Inconsistent Environment for Data Type");
#endif
  }
  else {
    material_factory = std::shared_ptr<nimble::MaterialFactoryBase>(new nimble::MaterialFactory);
  }

  //
  // Blocks
  //
  if (parser_.UseKokkos()) {
    Initialize_Blocks_Kokkos(material_factory);
  }
  else {
    Initialize_Blocks(material_factory);
  }

}


void DataManager::Initialize_Blocks(const std::shared_ptr<MaterialFactoryType>& material_factory)
{
  const auto num_blocks = static_cast<int>(mesh_.GetNumBlocks());

  auto macro_data = dynamic_cast< nimble::ModelData* >(macroscale_data_.get());
  auto material_factory_ptr = dynamic_cast< nimble::MaterialFactory* >(material_factory.get());

  std::map<int, nimble::Block>& blocks = macro_data->GetBlocks();
  std::map<int, nimble::Block>::iterator block_it;
  std::vector<int> block_ids = mesh_.GetBlockIds();
  for (int i=0 ; i< num_blocks ; i++){
    int block_id = block_ids[i];
    std::string const & macro_material_parameters = parser_.GetMacroscaleMaterialParameters(block_id);
    std::map<int, std::string> const & rve_material_parameters = parser_.GetMicroscaleMaterialParameters();
    std::string rve_bc_strategy = parser_.GetMicroscaleBoundaryConditionStrategy();
    blocks[block_id] = nimble::Block();
    blocks[block_id].Initialize(macro_material_parameters, rve_material_parameters, rve_mesh_, rve_bc_strategy, *material_factory_ptr);
    std::vector< std::pair<std::string, nimble::Length> > data_labels_and_lengths;
    blocks[block_id].GetDataLabelsAndLengths(data_labels_and_lengths);
    macro_data->DeclareElementData(block_id, data_labels_and_lengths);
  }


#ifdef NIMBLE_HAVE_UQ
  // configure & allocate
  if (parser_.HasUq())
  {
    uq_model_ = std::shared_ptr< nimble::UqModel >(new nimble::UqModel(dim,num_nodes,num_blocks));
    uq_model.ParseConfiguration(parser->UqModelString());
    std::map<std::string, std::string> lines = parser_.UqParamsStrings();
    for(std::map<std::string, std::string>::iterator it = lines.begin(); it != lines.end(); it++){
      std::string material_key = it->first;
      int block_id = parser->GetBlockIdFromMaterial( material_key );
      std::string uq_params_this_material = it->second;
      uq_model.ParseBlockInput( uq_params_this_material, block_id, blocks[block_id] );
    }
    // initialize
    uq_model.Initialize(*mesh_,macro_data);
  }
#endif

  std::map<int, int> num_elem_in_each_block = mesh_.GetNumElementsInBlock();
  macro_data->AllocateElementData(num_elem_in_each_block);
  macro_data->SpecifyOutputFields(parser_.GetOutputFieldString());
  std::map<int, std::vector<std::string> > const & elem_data_labels = macroscale_data_->GetElementDataLabels();
  std::map<int, std::vector<std::string> > const & derived_elem_data_labels = macroscale_data_->GetDerivedElementDataLabelsForOutput();

  // Initialize the element data
  std::vector<int> rve_output_elem_ids = parser_.MicroscaleOutputElementIds();
  for (block_it=blocks.begin(); block_it!=blocks.end() ; block_it++) {
    int block_id = block_it->first;
    int num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    std::vector<int> const & elem_global_ids = mesh_.GetElementGlobalIdsInBlock(block_id);
    nimble::Block& block = block_it->second;
    std::vector<double> & elem_data_n = macro_data->GetElementDataOld(block_id);
    std::vector<double> & elem_data_np1 = macro_data->GetElementDataNew(block_id);
    block.InitializeElementData(num_elem_in_block,
                                elem_global_ids,
                                rve_output_elem_ids,
                                elem_data_labels.at(block_id),
                                derived_elem_data_labels.at(block_id),
                                elem_data_n,
                                elem_data_np1,
                                *material_factory_ptr,
                                *this);
  }

}


void DataManager::Initialize_Blocks_Kokkos(
    const std::shared_ptr<MaterialFactoryType>& material_factory
)
{
#ifdef NIMBLE_HAVE_KOKKOS
  bool store_unrotated_stress(true);

  const auto num_blocks = static_cast<int>(mesh_.GetNumBlocks());

  auto model_data = dynamic_cast< nimble_kokkos::ModelData* >(macroscale_data_.get());
  auto material_factory_ptr = dynamic_cast< nimble_kokkos::MaterialFactory* >(material_factory.get());

  //
  // Blocks
  //
  auto &blocks = model_data->GetBlocks();
  std::vector<int> block_ids = mesh_.GetBlockIds();
  for (int i=0 ; i<num_blocks ; i++){
    int block_id = block_ids.at(i);
    std::string const & macro_material_parameters = parser_.GetMacroscaleMaterialParameters(block_id);
    std::map<int, std::string> const & rve_material_parameters = parser_.GetMicroscaleMaterialParameters();
    std::string rve_bc_strategy = parser_.GetMicroscaleBoundaryConditionStrategy();
    int num_elements_in_block = mesh_.GetNumElementsInBlock(block_id);
    blocks[block_id] = nimble_kokkos::Block();
    blocks.at(block_id).Initialize(macro_material_parameters, num_elements_in_block, *material_factory_ptr);
    //
    // MPI version use model_data.DeclareElementData(block_id, data_labels_and_lengths);
    //
    std::vector<double> initial_value(9, 0.0);
    initial_value[0] = initial_value[1] = initial_value[2] = 1.0;
    field_ids_.deformation_gradient = model_data->AllocateIntegrationPointData(block_id, nimble::FULL_TENSOR,
                                                                             "deformation_gradient",
                                                                             num_elements_in_block, initial_value);
    // volume-averaged quantities for I/O are stored as element data
    model_data->AllocateElementData(block_id, nimble::FULL_TENSOR, "deformation_gradient", num_elements_in_block);

    field_ids_.stress = model_data->AllocateIntegrationPointData(block_id, nimble::SYMMETRIC_TENSOR, "stress",
                                                               num_elements_in_block);
    if (store_unrotated_stress) {
      field_ids_.unrotated_stress = model_data->AllocateIntegrationPointData(block_id, nimble::SYMMETRIC_TENSOR, "stress",
                                                                           num_elements_in_block);
    }

    // volume-averaged quantities for I/O are stored as element data
    model_data->AllocateElementData(block_id, nimble::SYMMETRIC_TENSOR, "stress", num_elements_in_block);

    if (parser_.GetOutputFieldString().find("volume") != std::string::npos) {
      model_data->AllocateElementData(block_id, nimble::SCALAR, "volume", num_elements_in_block);
    }
  }

  // Initialize the exodus-output-manager
  auto &exodus_output_manager = model_data->GetExodusOutputManager();
  exodus_output_manager.SpecifyOutputFields(*model_data, parser_.GetOutputFieldString());

  model_data->SetNodeDataLabelsForOutput(exodus_output_manager.GetNodeDataLabelsForOutput());
  model_data->SetElementDataLabelsForOutput(exodus_output_manager.GetElementDataLabelsForOutput());

  std::map<int, std::vector<std::string> > derived_elem_data_labels;
  for (auto block_id : block_ids) {
    derived_elem_data_labels[block_id] = std::vector<std::string>(); // TODO eliminate this
  }

  model_data->SetDerivedElementDataLabelsForOutput(std::move(derived_elem_data_labels));

#endif
}


RVEData& DataManager::AllocateRVEData(int global_element_id,
                                      int integration_point_id) {
  std::pair<int, int> id_pair(global_element_id, integration_point_id);
  if (rve_data_.find(id_pair) != rve_data_.end()) {
    throw std::logic_error("\n****Error in DataManager::AllocateRVEData, requested global_element_id and integration_point_id have already been allocated.\n");
  }
  rve_data_[id_pair] = RVEData();
  RVEData& rve_data = rve_data_[id_pair];
  return rve_data;
}


RVEData& DataManager::GetRVEData(int global_element_id,
                                   int integration_point_id) {
  std::pair<int, int> id_pair(global_element_id, integration_point_id);
  RVEData& rve_data = rve_data_.at(id_pair);
  return rve_data;
}

} // namespace nimble
