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

#ifndef NIMBLE_MODEL_DATA_H
#define NIMBLE_MODEL_DATA_H

#include "nimble_block.h"
#include "nimble_defs.h"
#include "nimble_model_data_base.h"
#include "nimble_view.h"

#ifdef NIMBLE_HAVE_DARMA
#include "darma.h"
#else
#include <map>
#include <string>
#include <vector>
#endif

namespace nimble {

class DataManager;
class MaterialFactoryBase;
class VectorCommunicator;

class ModelData : public ModelDataBase
{
 public:
  ModelData() = default;

  ~ModelData() override = default;

  //--- Common interface from base class

  /// \brief Allocate data storage for a node-based quantity
  ///
  /// \param length
  /// \param label
  /// \param num_objects
  /// \return Field ID for the data allocated
  int
  AllocateNodeData(Length length, std::string label, int num_objects) override;

  /// \brief Returns the field ID for a specific label
  ///
  /// \param field_label Label for a stored quantity
  /// \return Field ID to identify the data storage
  int
  GetFieldId(const std::string& label) const override;

  /// \brief Initialize the different blocks in the mesh
  ///
  /// \param data_manager Reference to the data manager
  /// \param material_factory_base Shared pointer to the material factory
  void
  InitializeBlocks(nimble::DataManager& data_manager, const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base)
      override;

  /// \brief Copy time state (n+1) into time state (n)
  ///
  /// \param data_manager Reference to the data manager
  void
  UpdateStates(const nimble::DataManager& data_manager) override
  {
    element_data_n_.swap(element_data_np1_);
  }

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

  /// \brief Compute the external force
  ///
  /// \param data_manager Reference to the DataManager
  /// \param time_previous
  /// \param time_current
  /// \param is_output_step
  void
  ComputeExternalForce(
      nimble::DataManager& data_manager,
      double               time_previous,
      double               time_current,
      bool                 is_output_step) override;

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

  //--- Specific routines

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    ar | dim_ | critical_time_step_ | block_ids_ | blocks_ | data_fields_ | node_data_ | output_node_component_labels_ |
        element_data_fields_ | element_component_labels_ | element_data_n_ | element_data_np1_ |
        output_element_component_labels_ | derived_output_element_data_labels_ | globally_shared_nodes_ |
        global_node_id_to_local_node_id_;
  }
#endif

  Field
  GetField(int field_id);

  double*
  GetNodeData(int field_id);

  void
  GetNodeDataForOutput(std::vector<std::vector<double>>& single_component_arrays);

  void
  DeclareElementData(int block_id, std::vector<std::pair<std::string, Length>> const& data_labels_and_lengths);

  void
  AllocateElementData(std::map<int, int> const& num_integration_points_in_each_block);

  std::vector<double>&
  GetElementDataOld(int block_id)
  {
    return element_data_n_.at(block_id);
  }

  std::vector<double>&
  GetElementDataNew(int block_id)
  {
    return element_data_np1_.at(block_id);
  }

  void
  GetElementDataForOutput(std::map<int, std::vector<std::vector<double>>>& single_component_arrays);

  void
  SpecifyOutputFields(const std::string &output_field_string);

  std::map<int, std::shared_ptr< nimble::Block> >&
  GetBlocks()
  {
    return blocks_;
  }

  std::map<int, std::shared_ptr< nimble::Block> > const&
  GetBlocks() const
  {
    return blocks_;
  }

  std::vector<int>&
  GetGloballySharedNodes()
  {
    return globally_shared_nodes_;
  }

  std::map<int, int>&
  GetGlobalNodeIdToLocalNodeIdMap()
  {
    return global_node_id_to_local_node_id_;
  }

 protected:
  void
  AssignFieldId(Field& field);

  void
  GetNodeDataComponent(int field_id, int component, double* const component_data);

  double*
  GetNodeData(const std::string& label);

  template< typename TBlock >
  void EmplaceBlocks(nimble::DataManager& data_manager,
                     const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base);

  void AllocateInitializeElementData(nimble::DataManager&                                data_manager,
                                     const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base);

 protected:
  //! Block ids
  std::vector<int> block_ids_;

  //! Blocks
  std::map<int, std::shared_ptr<nimble::Block> > blocks_;

  //! Map key is the field_id, value is the corresponding Field.
  std::map<int, Field> data_fields_;

  //! Map key is the field_id, vector contains the nested data array for the
  //! given nodal field.
  std::map<int, std::vector<double>> node_data_;

  //! Map key is the block_id, the vector contains the field ids for the fields
  //! on that block.
  std::map<int, std::vector<int>> element_data_fields_;

  //! Map key is the block_id, vector contains full data array for that block at
  //! step N.
  std::map<int, std::vector<double>> element_data_n_;

  //! Map key is the block_id, vector contains full data array for that block at
  //! step N+1.
  std::map<int, std::vector<double>> element_data_np1_;

  //! List of node ids that are shared across multiple ranks
  std::vector<int> globally_shared_nodes_;

  //! Map from global node it to local node id
  std::map<int, int> global_node_id_to_local_node_id_;

  //! Information for Exodus output about node data
  std::vector<std::vector<double>> node_data_for_output_;

  //! Information for Exodus output about element data
  std::map<int, std::vector<std::vector<double>>> elem_data_for_output_;

  //! Information for Exodus output about element data
  std::map<int, std::vector<std::vector<double>>> derived_elem_data_;
};

}  // namespace nimble

#endif
