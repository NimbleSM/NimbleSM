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

#ifdef NIMBLE_HAVE_UQ
#include "nimble_uq.h"
#include "nimble_uq_model_data.h"
#endif

#ifdef NIMBLE_HAVE_KOKKOS
#include "nimble_kokkos_material_factory.h"
#include "nimble_kokkos_model_data.h"
#endif

#include <algorithm>
#include <stdexcept>

#include "nimble_block_material_interface_factory_base.h"
#include "nimble_boundary_condition_manager.h"
#include "nimble_material_factory.h"
#include "nimble_model_data.h"
#include "nimble_model_data_base.h"
#include "nimble_parser.h"
#include "nimble_vector_communicator.h"

namespace nimble {

DataManager::DataManager(
    const nimble::Parser&      parser,
    const nimble::GenesisMesh& mesh,
    const nimble::GenesisMesh& rve_mesh)
    : parser_(parser),
      mesh_(mesh),
      rve_mesh_(rve_mesh),
      macroscale_data_(),
      rve_data_(),
      field_ids_(),
      vector_communicator_(nullptr),
      boundary_condition_(new nimble::BoundaryConditionManager())
{
  Initialize();
}

void
DataManager::Initialize()
{
  const auto dim       = static_cast<int>(mesh_.GetDim());
  const auto num_nodes = static_cast<int>(mesh_.GetNumNodes());

  //--- Create VectorCommunicator
#ifdef NIMBLE_HAVE_TRILINOS
  auto comm = (parser_.UseTpetra()) ? Tpetra::getDefaultComm() : Teuchos::RCP<const Teuchos::Comm<int>>();
#else
  int comm = 0;
#endif
  vector_communicator_ = std::make_shared<nimble::VectorCommunicator>(dim, num_nodes, comm);

  std::vector<int> global_node_ids(num_nodes);
  int const* const global_node_ids_ptr = mesh_.GetNodeGlobalIds();
  for (int n = 0; n < num_nodes; ++n) { global_node_ids[n] = global_node_ids_ptr[n]; }

  // DJL
  // Here is where the initialization occurs for MPI operations
  // In this call, each rank determines which nodes are shared with which other
  // ranks This information is stored so that the vector reductions will work
  // later
  vector_communicator_->Initialize(global_node_ids);

  //--- Create ModelData
  if (parser_.UseUQ()) {
#ifdef NIMBLE_HAVE_UQ
    macroscale_data_ = std::make_shared<nimble_uq::ModelData>();
#else
    throw std::runtime_error(" Wrong environment !\n");
#endif
  } else if (parser_.UseKokkos()) {
#ifdef NIMBLE_HAVE_KOKKOS
    macroscale_data_ = std::make_shared<nimble_kokkos::ModelData>();
#else
    throw std::runtime_error(" Wrong environment !\n");
#endif
  } else {
    macroscale_data_ = std::make_shared<nimble::ModelData>();
  }

  macroscale_data_->SetDimension(dim);

  //
  // Initialize the boundary condition manager
  //

  std::map<int, std::string> const&      node_set_names          = mesh_.GetNodeSetNames();
  std::map<int, std::vector<int>> const& node_sets               = mesh_.GetNodeSets();
  std::vector<std::string> const&        bc_strings              = parser_.GetBoundaryConditionStrings();
  std::string                            time_integration_scheme = parser_.TimeIntegrationScheme();
  boundary_condition_->Initialize(node_set_names, node_sets, bc_strings, dim, time_integration_scheme);

  //
  // Initialize vectors for storing fields
  //

  if (time_integration_scheme == "explicit")
    field_ids_.lumped_mass = macroscale_data_->AllocateNodeData(nimble::SCALAR, "lumped_mass", num_nodes);

  field_ids_.reference_coordinates =
      macroscale_data_->AllocateNodeData(nimble::VECTOR, "reference_coordinate", num_nodes);
  field_ids_.displacement = macroscale_data_->AllocateNodeData(nimble::VECTOR, "displacement", num_nodes);
  field_ids_.velocity     = macroscale_data_->AllocateNodeData(nimble::VECTOR, "velocity", num_nodes);
  field_ids_.acceleration = macroscale_data_->AllocateNodeData(nimble::VECTOR, "acceleration", num_nodes);

  field_ids_.internal_force = macroscale_data_->AllocateNodeData(nimble::VECTOR, "internal_force", num_nodes);
  field_ids_.external_force = macroscale_data_->AllocateNodeData(nimble::VECTOR, "external_force", num_nodes);

  field_ids_.contact_force = macroscale_data_->AllocateNodeData(nimble::VECTOR, "contact_force", num_nodes);

  if (time_integration_scheme == "quasistatic") {
    //
    // These variables are used in the "quasi-static" simulations
    //
    macroscale_data_->AllocateNodeData(nimble::VECTOR, "trial_displacement", num_nodes);
    macroscale_data_->AllocateNodeData(nimble::VECTOR, "displacement_fluctuation", num_nodes);
    macroscale_data_->AllocateNodeData(nimble::VECTOR, "trial_internal_force", num_nodes);
    macroscale_data_->AllocateNodeData(nimble::SCALAR, "skin_node", num_nodes);
  }

  macroscale_data_->SetReferenceCoordinates(mesh_);

  //
  // Initialization of RVE data
  //
  rve_macroscale_deformation_gradient_.resize(dim * dim, 0.0);
  if (parser_.TimeIntegrationScheme() == "explicit") {
    for (int i = 0; i < dim; i++) rve_macroscale_deformation_gradient_[i] = 1.0;
  }
}

RVEData&
DataManager::AllocateRVEData(int global_element_id, int integration_point_id)
{
  std::pair<int, int> id_pair(global_element_id, integration_point_id);
  if (rve_data_.find(id_pair) != rve_data_.end()) {
    throw std::logic_error(
        "\n****Error in DataManager::AllocateRVEData, requested "
        "global_element_id and integration_point_id have already been "
        "allocated.\n");
  }
  rve_data_[id_pair] = RVEData();
  RVEData& rve_data  = rve_data_[id_pair];
  return rve_data;
}

RVEData&
DataManager::GetRVEData(int global_element_id, int integration_point_id)
{
  std::pair<int, int> id_pair(global_element_id, integration_point_id);
  RVEData&            rve_data = rve_data_.at(id_pair);
  return rve_data;
}

void
DataManager::InitializeOutput(const std::string& filename)
{
  std::vector<std::string> global_data_labels;

  exodus_output_ = std::shared_ptr<nimble::ExodusOutput>(new nimble::ExodusOutput);
  exodus_output_->Initialize(filename, mesh_);

  auto& node_data_labels_for_output = macroscale_data_->GetNodeDataLabelsForOutput();
  auto& elem_data_labels_for_output = macroscale_data_->GetElementDataLabelsForOutput();
  auto& derived_elem_data_labels    = macroscale_data_->GetDerivedElementDataLabelsForOutput();

  macroscale_data_->InitializeExodusOutput(*this);

  exodus_output_->InitializeDatabase(
      mesh_, global_data_labels, node_data_labels_for_output, elem_data_labels_for_output, derived_elem_data_labels);
}

void
DataManager::WriteOutput(double time_current)
{
  macroscale_data_->WriteExodusOutput(*this, time_current);
}

void
DataManager::SetBlockMaterialInterfaceFactory(
    const std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase>& block_material_factory)
{
  block_material_factory_ = block_material_factory;
}

const std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase>&
DataManager::GetBlockMaterialInterfaceFactory() const
{
  return block_material_factory_;
}

}  // namespace nimble
