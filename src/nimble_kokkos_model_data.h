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
#include <string>
#include <iostream>
#include <vector>
#include <memory>

#include "nimble_defs.h"
#include "nimble_exodus_output_manager.h"
#include "nimble_data_utils.h"
#include "nimble_model_data_base.h"


namespace nimble {

class GenesisMesh;
class DataManager;

}

namespace nimble_kokkos {

class ModelData : public nimble::ModelDataBase
{

 public:

  ModelData() = default;

  ~ModelData() override = default;

  //--- Common interface from nimble::ModelDataBase

  /// \brief Allocate data storage for a node-based quantity
  ///
  /// \param length
  /// \param label
  /// \param num_objects
  /// \return Field ID for the data allocated
  int AllocateNodeData(nimble::Length length,
                       std::string label,
                       int num_objects) override;

  /// \brief Returns the field ID for a specific label
  ///
  /// \param field_label Label for a stored quantity
  /// \return Field ID to identify the data storage
  int GetFieldId(const std::string& field_label) const override
  { return field_label_to_field_id_map_.at(field_label); }

  /// \brief Set the reference coordinates
  ///
  /// \param mesh Reference to the global mesh
  void SetReferenceCoordinates(const nimble::GenesisMesh &mesh) override;

  /// \brief Initialize the different blocks in the mesh
  ///
  /// \param data_manager Reference to the data manager
  /// \param material_factory_base Shared pointer to the material factory
  void InitializeBlocks(nimble::DataManager &data_manager,
                        const std::shared_ptr<MaterialFactoryType> &material_factory_base) override;

  /// \brief Swap states between time n and time (n+1)
  ///
  /// \param data_manager Reference to the data manager
  virtual void SwapStates(const nimble::DataManager &data_manager) override;

  //--- Specific routines

  int AllocateElementData(int block_id,
                          nimble::Length length,
                          std::string label,
                          int num_objects);

  int AllocateIntegrationPointData(int block_id,
                                   nimble::Length length,
                                   std::string label,
                                   int num_objects,
                                   std::vector<double> initial_value = std::vector<double>());

  std::vector<int> GetBlockIds() const ;

  std::vector<std::string> GetScalarNodeDataLabels() const ;

  std::vector<std::string> GetVectorNodeDataLabels() const ;

  std::vector<std::string> GetSymmetricTensorIntegrationPointDataLabels(int block_id) const ;

  std::vector<std::string> GetFullTensorIntegrationPointDataLabels(int block_id) const ;

  HostScalarNodeView GetHostScalarNodeData(int field_id);

  HostVectorNodeView GetHostVectorNodeData(int field_id);

  HostSymTensorIntPtView GetHostSymTensorIntegrationPointData(int block_id,
                                                              int field_id,
                                                              nimble::Step step);

  HostFullTensorIntPtView GetHostFullTensorIntegrationPointData(int block_id,
                                                                int field_id,
                                                                nimble::Step step);

  HostScalarElemView GetHostScalarElementData(int block_id,
                                              int field_id);

  HostSymTensorElemView GetHostSymTensorElementData(int block_id,
                                                    int field_id);

  HostFullTensorElemView GetHostFullTensorElementData(int block_id,
                                                      int field_id);

  DeviceScalarNodeView GetDeviceScalarNodeData(int field_id);

  DeviceVectorNodeView GetDeviceVectorNodeData(int field_id);

  DeviceSymTensorIntPtView GetDeviceSymTensorIntegrationPointData(int block_id,
                                                                  int field_id,
                                                                  nimble::Step step);

  DeviceFullTensorIntPtView GetDeviceFullTensorIntegrationPointData(int block_id,
                                                                    int field_id,
                                                                    nimble::Step step);

  DeviceScalarElemView GetDeviceScalarElementData(int block_id,
                                                  int field_id);

  DeviceSymTensorElemView GetDeviceSymTensorElementData(int block_id,
                                                        int field_id);

  DeviceFullTensorElemView GetDeviceFullTensorElementData(int block_id,
                                                          int field_id);

  DeviceScalarNodeGatheredView GatherScalarNodeData(int field_id,
                                                    int num_elements,
                                                    int num_nodes_per_element,
                                                    DeviceElementConnectivityView elem_conn_d,
                                                    DeviceScalarNodeGatheredView gathered_view_d);

  DeviceVectorNodeGatheredView GatherVectorNodeData(int field_id,
                                                    int num_elements,
                                                    int num_nodes_per_element,
                                                    DeviceElementConnectivityView elem_conn_d,
                                                    DeviceVectorNodeGatheredView gathered_view_d);

  void ScatterScalarNodeData(int field_id,
                             int num_elements,
                             int num_nodes_per_element,
                             DeviceElementConnectivityView elem_conn_d,
                             DeviceScalarNodeGatheredView gathered_view_d);

  void ScatterVectorNodeData(int field_id,
                             int num_elements,
                             int num_nodes_per_element,
                             DeviceElementConnectivityView elem_conn_d,
                             DeviceVectorNodeGatheredView gathered_view_d);

#ifndef KOKKOS_ENABLE_QTHREADS
  void ScatterScalarNodeDataUsingKokkosScatterView(int field_id,
                                                   int num_elements,
                                                   int num_nodes_per_element,
                                                   DeviceElementConnectivityView elem_conn_d,
                                                   DeviceScalarNodeGatheredView gathered_view_d);
#endif

  std::map<int, nimble_kokkos::Block>& GetBlocks() { return blocks_; }

  nimble_kokkos::ExodusOutputManager& GetExodusOutputManager()
  { return exodus_output_manager_; }

 protected:

  using Data = std::unique_ptr< FieldBase >;

  std::map<std::string, int> field_label_to_field_id_map_;

  //! Blocks
  std::map<int, nimble_kokkos::Block> blocks_;

  std::vector< Data > host_node_data_;
  std::vector< Data > device_node_data_;
  std::map<int, int> field_id_to_host_node_data_index_;
  std::map<int, int> field_id_to_device_node_data_index_;

  std::map<int, int> block_id_to_element_data_index_;
  std::vector< std::vector< Data > > host_element_data_;
  std::vector< std::vector< Data > > device_element_data_;
  std::vector< std::map<int, int> > field_id_to_host_element_data_index_;
  std::vector< std::map<int, int> > field_id_to_device_element_data_index_;

  std::map<int, int> block_id_to_integration_point_data_index_;
  std::vector< std::vector< Data > > host_integration_point_data_step_n_;
  std::vector< std::vector< Data > > host_integration_point_data_step_np1_;
  std::vector< std::vector< Data > > device_integration_point_data_step_n_;
  std::vector< std::vector< Data > > device_integration_point_data_step_np1_;
  std::vector< std::map<int, int> > field_id_to_host_integration_point_data_index_;
  std::vector< std::map<int, int> > field_id_to_device_integration_point_data_index_;

  nimble_kokkos::ExodusOutputManager exodus_output_manager_;

};

} // namespace

#endif
