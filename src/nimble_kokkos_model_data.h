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

#ifndef NIMBLE_KOKKOS_MODEL_DATA_H
#define NIMBLE_KOKKOS_MODEL_DATA_H

#include <Kokkos_Core.hpp>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "nimble_block_material_interface_base.h"
#include "nimble_data_utils.h"
#include "nimble_defs.h"
#include "nimble_kokkos_block.h"
#include "nimble_model_data_base.h"

namespace nimble {

class GenesisMesh;
class DataManager;
class MaterialFactoryBase;

}  // namespace nimble

namespace nimble_kokkos {

class ExodusOutputManager;

class ModelData : public nimble::ModelDataBase
{
 public:
  ModelData();

  ~ModelData() override;

  //--- Common interface from nimble::ModelDataBase

  /// \brief Allocate data storage for a node-based quantity
  ///
  /// \param length
  /// \param label
  /// \param num_objects
  /// \return Field ID for the data allocated
  int
  AllocateNodeData(nimble::Length length, std::string label, int num_objects) override;

  /// \brief Returns the field ID for a specific label
  ///
  /// \param field_label Label for a stored quantity
  /// \return Field ID to identify the data storage
  int
  GetFieldId(const std::string& field_label) const override
  {
    return field_label_to_field_id_map_.at(field_label);
  }

  /// \brief Initialize the different blocks in the mesh
  ///
  /// \param data_manager Reference to the data manager
  /// \param material_factory_base Shared pointer to the material factory
  void
  InitializeBlocks(
      nimble::DataManager&                                data_manager,
      const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base) override;

  /// \brief Copy time state (n+1) into time state (n)
  ///
  /// \param data_manager Reference to the data manager
  void
  UpdateStates(const nimble::DataManager& data_manager) override;

  using ModelDataBase::GetScalarNodeData;
  using ModelDataBase::GetVectorNodeData;

  /// \brief Get view of scalar quantity defined on nodes
  ///
  /// \param field_id the field id (see DataManager::GetFieldIDs())
  /// \return Viewify<1> object for scalar quantity
  nimble::Viewify<1>
  GetScalarNodeData(int field_id) override;

  /// \brief Get view of vector quantity defined on nodes
  ///
  /// \param field_id the field id (see DataManager::GetFieldIDs())
  /// \return Viewify<2> object for vector quantity
  nimble::Viewify<2>
  GetVectorNodeData(int field_id) override;

  void
  ComputeLumpedMass(nimble::DataManager& data_manager) override;

  void
  InitializeExodusOutput(nimble::DataManager& data_manager) override;

  void
  WriteExodusOutput(nimble::DataManager& data_manager, double time_current) override;

  /// \brief Compute the internal force
  ///
  /// \param[in] data_manager
  /// \param[in] time_previous
  /// \param[in] time_current
  /// \param[in] is_output_step
  /// \param[in] displacement
  /// \param[out] internal_force  Output for internal force
  void
  ComputeInternalForce(
      nimble::DataManager&      data_manager,
      double                    time_previous,
      double                    time_current,
      bool                      is_output_step,
      const nimble::Viewify<2>& displacement,
      nimble::Viewify<2>&       force) override;

  /// \brief Update model with new velocity
  ///
  /// \param[in] data_manager Reference to the data manager
  /// \param[in] dt Current time step
  ///
  /// \note This routine wll synchronize the host and device velocities.
  ///
  void
  UpdateWithNewVelocity(nimble::DataManager& data_manager, double dt) override;

  /// \brief Update model with new displacement
  ///
  /// \param[in] data_manager Reference to the data manager
  /// \param[in] dt Current time step
  ///
  /// \note This routine wll synchronize the host and device displacements.
  ///
  void
  UpdateWithNewDisplacement(nimble::DataManager& data_manager, double dt) override;

  //--- Specific routines

  int
  AllocateElementData(int block_id, nimble::Length length, std::string label, int num_objects);

  int
  AllocateIntegrationPointData(
      int                 block_id,
      nimble::Length      length,
      std::string         label,
      int                 num_objects,
      std::vector<double> initial_value = std::vector<double>());

  std::vector<int>
  GetBlockIds() const;

  std::vector<std::string>
  GetScalarNodeDataLabels() const;

  std::vector<std::string>
  GetVectorNodeDataLabels() const;

  std::vector<std::string>
  GetSymmetricTensorIntegrationPointDataLabels(int block_id) const;

  std::vector<std::string>
  GetFullTensorIntegrationPointDataLabels(int block_id) const;

  HostScalarNodeView
  GetHostScalarNodeData(int field_id);

  HostVectorNodeView
  GetHostVectorNodeData(int field_id);

  HostSymTensorIntPtView
  GetHostSymTensorIntegrationPointData(int block_id, int field_id, nimble::Step step);

  HostFullTensorIntPtView
  GetHostFullTensorIntegrationPointData(int block_id, int field_id, nimble::Step step);

  HostScalarElemView
  GetHostScalarElementData(int block_id, int field_id);

  HostSymTensorElemView
  GetHostSymTensorElementData(int block_id, int field_id);

  HostFullTensorElemView
  GetHostFullTensorElementData(int block_id, int field_id);

  DeviceScalarNodeView
  GetDeviceScalarNodeData(int field_id);

  DeviceVectorNodeView
  GetDeviceVectorNodeData(int field_id);

  DeviceSymTensorIntPtView
  GetDeviceSymTensorIntegrationPointData(int block_id, int field_id, nimble::Step step);

  DeviceFullTensorIntPtView
  GetDeviceFullTensorIntegrationPointData(int block_id, int field_id, nimble::Step step);

  DeviceScalarIntPtView
  GetDeviceScalarIntegrationPointData(int block_id, int field_id, nimble::Step step);

  DeviceVectorIntPtView
  GetDeviceVectorIntegrationPointData(int block_id, int field_id, nimble::Step step);

  DeviceScalarElemView
  GetDeviceScalarElementData(int block_id, int field_id);

  DeviceSymTensorElemView
  GetDeviceSymTensorElementData(int block_id, int field_id);

  DeviceFullTensorElemView
  GetDeviceFullTensorElementData(int block_id, int field_id);

  DeviceScalarNodeGatheredView
  GatherScalarNodeData(
      int                                  field_id,
      int                                  num_elements,
      int                                  num_nodes_per_element,
      const DeviceElementConnectivityView& elem_conn_d,
      DeviceScalarNodeGatheredView         gathered_view_d);

  DeviceVectorNodeGatheredView
  GatherVectorNodeData(
      int                                  field_id,
      int                                  num_elements,
      int                                  num_nodes_per_element,
      const DeviceElementConnectivityView& elem_conn_d,
      DeviceVectorNodeGatheredView         gathered_view_d);

  void
  ScatterScalarNodeData(
      int                                  field_id,
      int                                  num_elements,
      int                                  num_nodes_per_element,
      const DeviceElementConnectivityView& elem_conn_d,
      const DeviceScalarNodeGatheredView&  gathered_view_d);

  void
  ScatterVectorNodeData(
      int                                  field_id,
      int                                  num_elements,
      int                                  num_nodes_per_element,
      const DeviceElementConnectivityView& elem_conn_d,
      const DeviceVectorNodeGatheredView&  gathered_view_d);

#ifndef KOKKOS_ENABLE_QTHREADS
  void
  ScatterScalarNodeDataUsingKokkosScatterView(
      int                                  field_id,
      int                                  num_elements,
      int                                  num_nodes_per_element,
      const DeviceElementConnectivityView& elem_conn_d,
      const DeviceScalarNodeGatheredView&  gathered_view_d);
#endif

 protected:
  void
  InitializeGatheredVectors(const nimble::GenesisMesh& mesh_);

  /// \brief Initialize block data for material information
  ///
  /// \param data_manager
  void
  InitializeBlockData(nimble::DataManager& data_manager);

  template <FieldType ft>
  typename Field<ft>::View
  GetDeviceElementData(int block_id, int field_id);

  template <FieldType ft>
  typename Field<ft>::View
  GetDeviceIntPointData(int block_id, int field_id, nimble::Step step);

 protected:
  using Data = std::unique_ptr<FieldBase>;

  std::map<std::string, int> field_label_to_field_id_map_;

  //! Blocks
  std::map<int, nimble_kokkos::Block> blocks_;

  std::vector<Data>  host_node_data_;
  std::vector<Data>  device_node_data_;
  std::map<int, int> field_id_to_host_node_data_index_;
  std::map<int, int> field_id_to_device_node_data_index_;

  std::map<int, int>              block_id_to_element_data_index_;
  std::vector<std::vector<Data>>  host_element_data_;
  std::vector<std::vector<Data>>  device_element_data_;
  std::vector<std::map<int, int>> field_id_to_host_element_data_index_;
  std::vector<std::map<int, int>> field_id_to_device_element_data_index_;

  std::map<int, int>              block_id_to_integration_point_data_index_;
  std::vector<std::vector<Data>>  host_integration_point_data_step_n_;
  std::vector<std::vector<Data>>  host_integration_point_data_step_np1_;
  std::vector<std::vector<Data>>  device_integration_point_data_step_n_;
  std::vector<std::vector<Data>>  device_integration_point_data_step_np1_;
  std::vector<std::map<int, int>> field_id_to_host_integration_point_data_index_;
  std::vector<std::map<int, int>> field_id_to_device_integration_point_data_index_;

  std::unique_ptr<nimble_kokkos::ExodusOutputManager> exodus_output_manager_;

  //--- Work arrays
  std::vector<nimble_kokkos::DeviceVectorNodeGatheredView> gathered_reference_coordinate_d;
  std::vector<nimble_kokkos::DeviceVectorNodeGatheredView> gathered_displacement_d;
  std::vector<nimble_kokkos::DeviceVectorNodeGatheredView> gathered_internal_force_d;
  std::vector<nimble_kokkos::DeviceVectorNodeGatheredView> gathered_contact_force_d;

  //--- Data for Exodus output
  nimble_kokkos::HostVectorNodeView   displacement_h_;
  nimble_kokkos::DeviceVectorNodeView displacement_d_;

  nimble_kokkos::HostVectorNodeView   velocity_h_;
  nimble_kokkos::DeviceVectorNodeView velocity_d_;

  //--- Block data for materials
  std::vector<nimble::BlockData> block_data_;
};

}  // namespace nimble_kokkos

#endif
