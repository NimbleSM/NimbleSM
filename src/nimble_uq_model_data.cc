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

#ifdef NIMBLE_HAVE_UQ

#include<map>

#include "nimble_boundary_condition_manager.h"
#include "nimble_data_manager.h"
#include "nimble_material_factory.h"
#include "nimble_genesis_mesh.h"
#include "nimble_model_data.h"
#include "nimble_parser.h"
#include "nimble_uq.h"
#include "nimble_uq_block.h"
#include "nimble_uq_model_data.h"
#include "nimble_vector_communicator.h"
#include "nimble_view.h"

namespace nimble_uq {

void
ModelData::InitializeBlocks(
    nimble::DataManager&                                data_manager,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base)
{
  //
  const auto& mesh_     = data_manager.GetMesh();
  const auto& parser_   = data_manager.GetParser();
  const auto& rve_mesh_ = data_manager.GetRVEMesh();

  auto material_factory_ptr = dynamic_cast<nimble::MaterialFactory*>(material_factory_base.get());

  nimble::ModelData::EmplaceBlocks< nimble_uq::Block >(data_manager, material_factory_base);

  std::map<int, int> num_elem_in_each_block = mesh_.GetNumElementsInBlock();
  AllocateElementData(num_elem_in_each_block);
  SpecifyOutputFields(parser_.GetOutputFieldString());
  std::map<int, std::vector<std::string>> const& elem_data_labels         = GetElementDataLabels();
  std::map<int, std::vector<std::string>> const& derived_elem_data_labels = GetDerivedElementDataLabelsForOutput();

  // Initialize the element data
  std::vector<int> rve_output_elem_ids = parser_.MicroscaleOutputElementIds();
  for (auto block_it : blocks_) {
    int                     block_id          = block_it.first;
    int                     num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    std::vector<int> const& elem_global_ids   = mesh_.GetElementGlobalIdsInBlock(block_id);
    nimble::Block&          block             = block_it.second;
    std::vector<double>&    elem_data_n       = GetElementDataOld(block_id);
    std::vector<double>&    elem_data_np1     = GetElementDataNew(block_id);
    block.InitializeElementData(
        num_elem_in_block,
        elem_global_ids,
        rve_output_elem_ids,
        elem_data_labels.at(block_id),
        derived_elem_data_labels.at(block_id),
        elem_data_n,
        elem_data_np1,
        *material_factory_ptr,
        data_manager);
  }
  //
  //
  //
  //
  // Make a copy of cast blocks
  //
  for (auto block_it : blocks_) {
    int id = block_it.first;
    nimble_uq::Block* cast_block = dynamic_cast<nimble_uq::Block*>(&block_it.second);
    // Add test if nullptr
    uq_blocks_[id] = cast_block;
  }
  int num_nodes = mesh_.GetNumNodes();
  int num_blocks = mesh_.GetNumBlocks();
  std::vector<int> block_ids = mesh_.GetBlockIds();
  uq_model_ = std::shared_ptr<nimble::UqModel>(new nimble::UqModel(dim_, &mesh_, this));
  uq_model_->ParseConfiguration(parser_.UqModelString());
  std::map<std::string, std::string> lines = parser_.UqParamsStrings();
  for (auto & line : lines) {
    std::string material_key            = line.first;
    int         block_id                = parser_.GetBlockIdFromMaterial(material_key);
    std::string uq_params_this_material = line.second;
    std::string const & nominal_params_string = parser_.GetMacroscaleMaterialParameters(block_id);
    bool block_id_present = std::find(block_ids.begin(), block_ids.end(), block_id) != block_ids.end() ;
    uq_model_->ParseBlockInput( uq_params_this_material, block_id, nominal_params_string,
                               material_factory_base, block_id_present, uq_blocks_ );
  }
  //
  // initialize
  //
  uq_model_->Initialize();
  //
  uq_model_->Setup();
  /// FOR CALL TO BCS ===================
  int num_samples = uq_model_->GetNumSamples();
  for (int nuq = 0; nuq < num_samples; nuq++) {
    double* v = uq_model_->Velocities()[nuq];
    //
    // We should replace {0, 3} with {length, 3}
    //
    bc_offnom_velocity_views_.push_back(nimble::Viewify<2>(v, {0, 3}, {3, 1}));
  }
}

void
ModelData::WriteExodusOutput(nimble::DataManager& data_manager, double time_current)
{
  nimble::ModelData::WriteExodusOutput(data_manager, time_current);

  if (!uq_model_->Initialized()) return;

  const auto& mesh_     = data_manager.GetMesh();
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  for (auto block_it : blocks_) {
    int            block_id          = block_it.first;
    int            num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    int const*     elem_conn         = mesh_.GetConnectivity(block_id);
    nimble::Block& block             = block_it.second;
    uq_model_->PerformAnalyses(reference_coordinate.data(), num_elem_in_block, elem_conn, block_id, block);
  }
  // UH -- DEBUG TO FIX
  //uq_model_->Write(step);
}

void
ModelData::ApplyInitialConditions(nimble::DataManager& data_manager)
{
  auto bc                   = data_manager.GetBoundaryConditionManager();
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto velocity             = GetVectorNodeData("velocity");
  bc->ApplyInitialConditions(reference_coordinate, velocity, bc_offnom_velocity_views_);
}

void
ModelData::ApplyKinematicConditions(nimble::DataManager& data_manager, double time_current, double time_previous)
{
  auto bc                   = data_manager.GetBoundaryConditionManager();
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto displacement         = GetVectorNodeData("displacement");
  auto velocity             = GetVectorNodeData("velocity");
  bc->ApplyKinematicBC(
      time_current, time_previous, reference_coordinate, displacement, velocity, bc_offnom_velocity_views_);
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

  auto& rve_macroscale_deformation_gradient = data_manager.GetRVEDeformationGradient();

  for (auto& block_it : uq_blocks_) {
    int                        block_id          = block_it.first;
    int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    int const*                 elem_conn         = mesh.GetConnectivity(block_id);
    std::vector<int> const&    elem_global_ids   = mesh.GetElementGlobalIdsInBlock(block_id);
    auto                       block             =  block_it.second;
    std::vector<double> const& elem_data_n       = GetElementDataOld(block_id);
    std::vector<double>&       elem_data_np1     = GetElementDataNew(block_id);
    bool                       is_off_nominal    = false;  // HACK
    std::vector<double>        uq_params_this_sample;      // HACK
    block->ComputeInternalForce(
        reference_coord,
        displacement.data(),
        velocity,
        rve_macroscale_deformation_gradient.data(),
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
        is_output_step,
        is_off_nominal,         // UQ
        uq_params_this_sample  // UQ
        );
  }

  // Perform a vector reduction on internal force.  This is a vector nodal
  // quantity.
  auto          vector_comm      = data_manager.GetVectorCommunicator();
  constexpr int vector_dimension = 3;
  vector_comm->VectorReduction(vector_dimension, force.data());
}

void
ModelData::UpdateWithNewVelocity(nimble::DataManager& data_manager, double dt)
{
  uq_model_->UpdateVelocity(dt);
}

void
ModelData::UpdateWithNewDisplacement(nimble::DataManager& data_manager, double dt)
{
  // advance approximate trajectories
  uq_model_->UpdateDisplacement(dt);
  uq_model_->Prep();
}

}  // namespace nimble_uq

#endif
