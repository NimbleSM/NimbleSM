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
#include "nimble_data_utils.h"
#include "nimble_defs.h"
#include "nimble_exodus_output.h"
#include "nimble_genesis_mesh.h"
#include "nimble_linear_solver.h"
#include "nimble_model_data_base.h"

#ifdef NIMBLE_HAVE_DARMA
#include "darma.h"
#else
#include <vector>
#include <string>
#include <map>
#endif

namespace nimble {

class DataManager;

class ModelData : public ModelDataBase {

public:

  ModelData() = default;

  ~ModelData() override = default;

  //--- Common interface from nimble::ModelDataBase

  int GetFieldId(const std::string& label) const override;

  int AllocateNodeData(Length length,
                       std::string label,
                       int num_objects) override;

  void SetReferenceCoordinates(const nimble::GenesisMesh &mesh) override;

  void InitializeBlocks(nimble::DataManager &data_manager,
                        const std::shared_ptr<MaterialFactoryType> &material_factory_base) override;

  //--- Specific routines

#ifdef NIMBLE_HAVE_DARMA
  template<typename ArchiveType>
    void serialize(ArchiveType& ar) {
      ar | dim_ | critical_time_step_ | block_ids_ | blocks_ | data_fields_ | node_data_ | output_node_component_labels_ | element_data_fields_ | element_component_labels_ | element_data_n_ | element_data_np1_ | output_element_component_labels_ | derived_output_element_data_labels_ | globally_shared_nodes_ | global_node_id_to_local_node_id_;
    }
#endif

  Field GetField(int field_id);

  double* GetNodeData(int field_id);

  void GetNodeDataForOutput(std::vector< std::vector<double> >& single_component_arrays);

  void DeclareElementData(int block_id,
                          std::vector< std::pair<std::string, Length> > const & data_labels_and_lengths);

  void AllocateElementData(std::map<int, int> const & num_integration_points_in_each_block);

  std::vector<double> & GetElementDataOld(int block_id) {
    return element_data_n_.at(block_id);
  }

  std::vector<double> & GetElementDataNew(int block_id) {
    return element_data_np1_.at(block_id);
  }

  void GetElementDataForOutput(std::map<int, std::vector< std::vector<double> > >& single_component_arrays);

  void SpecifyOutputFields(std::string output_field_string);

  std::map<int, Block>& GetBlocks() { return blocks_; }

  std::map<int, Block> const & GetBlocks() const { return blocks_; }

  std::vector<int>& GetGloballySharedNodes() { return globally_shared_nodes_; }

  std::map<int, int>& GetGlobalNodeIdToLocalNodeIdMap() { return global_node_id_to_local_node_id_; }

  void SwapStates() {
    element_data_n_.swap(element_data_np1_);
  }

protected:

  void AssignFieldId(Field& field);

  void GetNodeDataComponent(int field_id,
                            int component,
                            double* const component_data);

protected:

  //! Block ids
  std::vector<int> block_ids_;

  //! Blocks
  std::map<int, Block> blocks_;

  //! Map key is the field_id, value is the corresponding Field.
  std::map<int, Field> data_fields_;

  //! Map key is the field_id, vector contains the nested data array for the given nodal field.
  std::map<int, std::vector<double>> node_data_;

  //! Map key is the block_id, the vector contains the field ids for the fields on that block.
  std::map<int, std::vector<int>> element_data_fields_;

  //! Map key is the block_id, vector contains full data array for that block at step N.
  std::map<int, std::vector<double> > element_data_n_;

  //! Map key is the block_id, vector contains full data array for that block at step N+1.
  std::map<int, std::vector<double> > element_data_np1_;

  //! List of node ids that are shared across multiple ranks
  std::vector<int> globally_shared_nodes_;

  //! Map from global node it to local node id
  std::map<int, int> global_node_id_to_local_node_id_;

};

} // namespace nimble

#endif
