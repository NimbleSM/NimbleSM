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

#ifdef NIMBLE_HAVE_UQ
#include "nimble_uq.h"
#endif

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
