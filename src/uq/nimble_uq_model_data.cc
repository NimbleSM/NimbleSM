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

#include "nimble_uq_model_data.h"

#include <map>

#include "nimble_boundary_condition_manager.h"
#include "nimble_data_manager.h"
#include "nimble_genesis_mesh.h"
#include "nimble_material_factory.h"
#include "nimble_model_data.h"
#include "nimble_parser.h"
#include "nimble_uq.h"
#include "nimble_uq_block.h"
#include "nimble_vector_communicator.h"
#include "nimble_view.h"

namespace nimble_uq {

void
ModelData::InitializeBlocks(
    nimble::DataManager&                                data_manager,
    const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base)
{
  //
  const auto& mesh_   = data_manager.GetMesh();
  const auto& parser_ = data_manager.GetParser();

  nimble::ModelData::EmplaceBlocks<nimble_uq::Block>(data_manager, material_factory_base);

  std::vector<int> block_ids = mesh_.GetBlockIds();
  uq_model_                  = std::shared_ptr<nimble::UqModel>(new nimble::UqModel(dim_, &mesh_, this));
  uq_model_->ParseConfiguration(parser_.UqModelString());
  std::map<std::string, std::string> lines = parser_.UqParamsStrings();
  for (auto& line : lines) {
    std::string        material_key            = line.first;
    int                block_id                = parser_.GetBlockIdFromMaterial(material_key);
    std::string        uq_params_this_material = line.second;
    std::string const& nominal_params_string   = parser_.GetModelMaterialParameters(block_id);
    bool               block_id_present = std::find(block_ids.begin(), block_ids.end(), block_id) != block_ids.end();
    uq_model_->ParseBlockInput(uq_params_this_material, block_id, nominal_params_string, material_factory_base, block_id_present, blocks_);
  }
  //
  uq_model_->Initialize();
  //
  std::cout << "START with AllocateInitializeElementData\n" << std::flush;
  AllocateInitializeElementData(data_manager, material_factory_base);
  std::cout << "DONE with AllocateInitializeElementData\n" << std::flush;
  //
  uq_model_->Setup();
  int num_samples = uq_model_->GetNumSamples();
  for (int i = 0; i < num_samples; i++) {
    double* u = uq_model_->Displacements()[i];
    double* v = uq_model_->Velocities()[i];
    double* f = uq_model_->Forces()[i];
    //
    // We should replace {0, 3} with {length, 3} ASK ULRICH
    //
    offnom_displacement_views_.push_back(nimble::Viewify<2>(u, {0, 3}, {3, 1}));
    offnom_velocity_views_.    push_back(nimble::Viewify<2>(v, {0, 3}, {3, 1}));
    offnom_force_views_.       push_back(nimble::Viewify<2>(f, {0, 3}, {3, 1}));
  }
}

void
ModelData::WriteExodusOutput(nimble::DataManager& data_manager, double time_current)
{
  nimble::ModelData::WriteExodusOutput(data_manager, time_current);

  if (!uq_model_->Initialized()) return;

  const auto& mesh_                = data_manager.GetMesh();
  auto        reference_coordinate = GetVectorNodeData("reference_coordinate");
  for (auto block_it : blocks_) {
    int        block_id          = block_it.first;
    int        num_elem_in_block = mesh_.GetNumElementsInBlock(block_id);
    int const* elem_conn         = mesh_.GetConnectivity(block_id);
    auto&      block             = block_it.second;
    uq_model_->PerformAnalyses(reference_coordinate.data(), num_elem_in_block, elem_conn, block_id, block);
  }

  uq_model_->Write(time_current);
}

void
ModelData::ApplyInitialConditions(nimble::DataManager& data_manager)
{
  auto bc                   = data_manager.GetBoundaryConditionManager();
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto velocity             = GetVectorNodeData("velocity");
  bc->ApplyInitialConditions(reference_coordinate, velocity, offnom_velocity_views_); // applied to all
}

void
ModelData::ApplyKinematicConditions(nimble::DataManager& data_manager, double time_current, double time_previous)
{
  auto bc                   = data_manager.GetBoundaryConditionManager();
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto displacement         = GetVectorNodeData("displacement");
  auto velocity             = GetVectorNodeData("velocity");
  bc->ApplyKinematicBC(
      time_current, time_previous, reference_coordinate, displacement, velocity, offnom_velocity_views_); // applied to all
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
  int num_samples = uq_model_->GetNumSamples();
  int num_exact_samples = uq_model_->GetNumExactSamples();

  force.zero();
  for(int i=0; i < num_samples; i++){ offnom_force_views_[i].zero(); }

  auto reference_coord = GetNodeData("reference_coordinate");
  auto velocity        = GetNodeData("velocity"); // ASK double * not View?

  for (auto& block_it : blocks_) {
    int                        block_id          = block_it.first;
    int                        num_elem_in_block = mesh.GetNumElementsInBlock(block_id);
    int const*                 elem_conn         = mesh.GetConnectivity(block_id);
    std::vector<int> const&    elem_global_ids   = mesh.GetElementGlobalIdsInBlock(block_id);
    std::vector<double> const& elem_data_n       = GetElementDataOld(block_id);
    std::vector<double>&       elem_data_np1     = GetElementDataNew(block_id);
    auto                       block             = dynamic_cast<nimble_uq::Block*>(block_it.second.get());
//  std::vector<double> const parameters();
    for(int i=0; i < num_exact_samples; i++){
      int ii = i-1;
      bool is_off_nominal = (i > 0);
      is_off_nominal = false; // HACK <<<<<<
      auto u = displacement;
      auto v = velocity;
      auto f = force;
      if (is_off_nominal) {
        u  = offnom_displacement_views_[ii];
        v  = uq_model_->Velocities()[ii];
        f  = offnom_force_views_[ii];
//      parameters  = uq_model_->GetParameters(ii); 
      }
      ii = (ii > -1) ? ii : 0; // HACK
      std::vector<double> const parameters  = uq_model_->GetParameters(ii);  // HACK ii = -1 for nominal
      block->ComputeInternalForce(
        reference_coord,
        u.data(),
        v,
        f.data(),
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
        is_off_nominal, parameters  // UQ
      );
    }
  }

  // Perform a vector reduction on internal force.  This is a vector nodal
  // quantity.
  auto          vector_comm      = data_manager.GetVectorCommunicator();
  constexpr int vector_dimension = 3;
  vector_comm->VectorReduction(vector_dimension, force.data());
  int num_exact_trajectories = uq_model_->GetNumExactSamples();
  for(int i=0; i <= num_exact_trajectories; i++){
    vector_comm->VectorReduction(vector_dimension,offnom_force_views_[i].data());
  }

  // Now apply closure to estimate approximate forces from the exact samples
  uq_model_->ApplyClosure();  
}

// NOTE need data manager?
// advance approximate trajectories
void
ModelData::UpdateWithNewVelocity(nimble::DataManager& data_manager, double dt)
{
//if (!uq_model_->initialized()) { return; }
  auto mass  = GetNodeData("lumped_mass"); // ASK double * not View?
  int n = uq_model_->GetNumSamples();
  for (int s = 0; s < n; ++s) {
    double* f = uq_model_->Forces()[s];
    double* v = uq_model_->Velocities()[s];
    for (int i = 0; i < 3; ++i) {
      double m = mass[i / 3];
      double a = (1.0 / m) * f[i]; // NOTE store inv_mass
      v[i] += dt * a;
    }
  }

}

// advance approximate trajectories
void
ModelData::UpdateWithNewDisplacement(nimble::DataManager& data_manager, double dt)
{
//if (!uq_model_->initialized()) { return; }
  int n = uq_model_->GetNumSamples();
  for (int s = 0; s < n; ++s) {
    double* u = uq_model_->Displacements()[s];
    double* v = uq_model_->Velocities()[s];
    for (int i = 0; i < 3; ++i) { u[i] += dt * v[i]; }
  }
}

}  // namespace nimble_uq

#endif
