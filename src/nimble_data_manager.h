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

#ifndef NIMBLE_DATA_MANAGER_H
#define NIMBLE_DATA_MANAGER_H

#include "nimble_data_utils.h"
#include "nimble_block.h"
#include "nimble_exodus_output.h"
#include "nimble_linear_solver.h"

#ifdef NIMBLE_HAVE_DARMA
  #include "darma.h"
#else
  #include <vector>
  #include <string>
  #include <map>
#endif

#ifdef NIMBLE_HAVE_UQ
  #include "nimble_uq.h"
#endif

namespace nimble {

  class ModelData {

  public:

    ModelData() : dim_(3), critical_time_step_(0.0) 
#ifdef NIMBLE_HAVE_UQ
                , uq_model_()
#endif
                {}

    virtual ~ModelData() {}

#ifdef NIMBLE_HAVE_DARMA
    template<typename ArchiveType>
    void serialize(ArchiveType& ar) {
      ar | dim_ | critical_time_step_ | block_ids_ | blocks_ | data_fields_ | node_data_ | output_node_component_labels_ | element_data_fields_ | element_component_labels_ | element_data_n_ | element_data_np1_ | output_element_component_labels_ | derived_output_element_data_labels_ | globally_shared_nodes_ | global_node_id_to_local_node_id_;
    }
#endif

    void SetDimension(int dim);

    int GetDimension() const {return dim_; }

    void SetCriticalTimeStep(double time_step) { critical_time_step_ = time_step; }

    double GetCriticalTimeStep() const { return critical_time_step_; }

    int GetFieldId(std::string label);

    Field GetField(int field_id);

    int AllocateNodeData(Length length,
                         std::string label,
                         int num_objects);

    double* GetNodeData(int field_id);

    std::vector<std::string> const & GetNodeDataLabelsForOutput() const {
      return output_node_component_labels_;
    }

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

    std::map<int, std::vector<std::string> > GetElementDataLabels() {
      return element_component_labels_;
    }

    std::map<int, std::vector<std::string> > GetElementDataLabelsForOutput() {
      return output_element_component_labels_;
    }

    void GetElementDataForOutput(std::map<int, std::vector< std::vector<double> > >& single_component_arrays);

    void SpecifyOutputFields(std::string output_field_string);

    std::map<int, std::vector<std::string> > GetDerivedElementDataLabelsForOutput() {
      return derived_output_element_data_labels_;
    }

    std::map<int, Block>& GetBlocks() { return blocks_; }

    std::map<int, Block> const & GetBlocks() const { return blocks_; }

    std::vector<int>& GetGloballySharedNodes() { return globally_shared_nodes_; }

    std::map<int, int>& GetGlobalNodeIdToLocalNodeIdMap() { return global_node_id_to_local_node_id_; }

    void SwapStates() {
      element_data_n_.swap(element_data_np1_);
    }

#ifdef NIMBLE_HAVE_UQ
    void SetUqModel(std::string uq_model_string) {
      UqModel model(uq_model_string);
      uq_model_ = model;
    }

    UqModel const & GetUqModel() const { return uq_model_;}
#endif

  protected:

    void AssignFieldId(Field& field);

    void GetNodeDataComponent(int field_id,
                              int component,
                              double* const component_data);

    //! Problem dimension, either 2 or 3.
    int dim_;

    //! Critical time step
    double critical_time_step_;

    //! Block ids
    std::vector<int> block_ids_;

    //! Blocks
    std::map<int, Block> blocks_;

    //! Map key is the field_id, value is the corresponding Field.
    std::map<int, Field> data_fields_;

    //! Map key is the field_id, vector contains the nested data array for the given nodal field.
    std::map<int, std::vector<double>> node_data_;

    //! Output labels for node data that will be written to disk
    std::vector<std::string> output_node_component_labels_;

    //! Map key is the block_id, the vector contains the field ids for the fields on that block.
    std::map<int, std::vector<int>> element_data_fields_;

    //! Map key is the block_id, vector contains component-wise label for each scalar entry in the data array.
    std::map<int, std::vector<std::string>> element_component_labels_;

    //! Map key is the block_id, vector contains full data array for that block at step N.
    std::map<int, std::vector<double> > element_data_n_;

    //! Map key is the block_id, vector contains full data array for that block at step N+1.
    std::map<int, std::vector<double> > element_data_np1_;

    //! Output labels for element data that will be written to disk.
    std::map<int, std::vector<std::string> > output_element_component_labels_;

    //! Output labels for derived element data that will be written to disk.
    std::map<int, std::vector<std::string> > derived_output_element_data_labels_;

    //! List of node ids that are shared across multiple ranks
    std::vector<int> globally_shared_nodes_;

    //! Map from global node it to local node id
    std::map<int, int> global_node_id_to_local_node_id_;

#ifdef NIMBLE_HAVE_UQ
    //struct setting up the macroscale UQ model
    UqModel uq_model_;
#endif

  };

  struct RVEData {
    ModelData model_data_;
    std::vector<double> residual_vector_;
    std::vector<double> linear_solver_solution_;
    CRSMatrixContainer tangent_stiffness_;
    bool write_exodus_output_;
    ExodusOutput exodus_output_;
    std::map<int, std::vector< std::vector<double> > > derived_elem_data_;

#ifdef NIMBLE_HAVE_DARMA
    template<typename ArchiveType>
    void serialize(ArchiveType& ar) {
      ar | model_data_ | residual_vector_ | linear_solver_solution_ | tangent_stiffness_ | write_exodus_output_ | exodus_output_ | derived_elem_data_;
    }
#endif
  };

 class DataManager {

  public:

    DataManager() {}

    virtual ~DataManager() {}

#ifdef NIMBLE_HAVE_DARMA
    template<typename ArchiveType>
    void serialize(ArchiveType& ar) {
      ar | macroscale_data_ | rve_data_;
    }
#endif

    ModelData& GetMacroScaleData();

    RVEData& AllocateRVEData(int global_element_id,
                             int integration_point_id);

    RVEData& GetRVEData(int global_element_id,
                        int integration_point_id);

 protected:

    ModelData macroscale_data_;
    std::map<std::pair<int,int>, RVEData> rve_data_;
 };

} // namespace nimble

#endif
