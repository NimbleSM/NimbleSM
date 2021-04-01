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

#include "nimble_block.h"
#include "nimble_data_manager.h"
#include "nimble_model_data.h"
#include "nimble_uq.h"
#include "nimble_view.h"

namespace nimble_uq {

void
ModelData::InitializeBlocks(
    nimble::DataManager&                        data_manager,
    const std::shared_ptr<MaterialFactoryType>& material_factory_base)
{
  nimble::ModelData::InitializeBlocks(data_manager, material_factory_base);
  //
  uq_model_ = std::shared_ptr<nimble::UqModel>(new nimble::UqModel(dim, num_nodes, num_blocks));
  uq_model_->ParseConfiguration(parser->UqModelString());
  std::map<std::string, std::string> lines = parser_.UqParamsStrings();
  for (std::map<std::string, std::string>::iterator it = lines.begin(); it != lines.end(); it++) {
    std::string material_key            = it->first;
    int         block_id                = parser_.GetBlockIdFromMaterial(material_key);
    std::string uq_params_this_material = it->second;
    uq_model_->ParseBlockInput(uq_params_this_material, block_id, blocks[block_id]);
  }
  //
  // initialize
  //
  uq_model_->Initialize(mesh_, this);
  //
  uq_model_->Setup();
  /// FOR CALL TO BCS ===================
  int num_samples = uq_model.GetNumSamples();
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

  for (block_it = blocks.begin(); block_it != blocks.end(); block_it++) {
    int            block_id          = block_it->first;
    int            num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    int const*     elem_conn         = mesh.GetConnectivity(block_id);
    nimble::Block& block             = block_it->second;
    uq_model_->PerformAnalyses(reference_coordinate, num_elem_in_block, elem_conn, block_id, block);
  }
  uq_model_->Write(step);
}

void
ModelData::ApplyInitialConditions(DataManager& data_manager)
{
  auto bc                   = data_manager.GetBoundaryConditionManager();
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto velocity             = GetVectorNodeData("velocity");
  bc->ApplyInitialConditions(reference_coordinate, velocity, bc_offnom_velocity_views_);
}

void
ModelData::ApplyKinematicConditions(DataManager& data_manager, double time_current, double time_previous)
{
  auto bc                   = data_manager.GetBoundaryConditionManager();
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto displacement         = GetVectorNodeData("displacement");
  auto velocity             = GetVectorNodeData("velocity");
  bc->ApplyKinematicBC(
      time_current, time_previous, reference_coordinate, displacement, velocity, bc_offnom_velocity_views_);
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
