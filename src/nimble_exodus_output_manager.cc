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

#include "nimble_exodus_output_manager.h"
#include <sstream>

#ifndef NIMBLE_HAVE_DARMA
  #include <vector>
#endif

#include <iostream>

namespace nimble_kokkos {

  void ExodusOutputManager::SpecifyOutputFields(nimble_kokkos::ModelData& model_data,
                                                std::string const & output_command_string) {

    // Parse the string into individual field labels
    std::vector<std::string> requested_labels;
    std::stringstream ss(output_command_string);
    std::string entry;
    while (std::getline(ss, entry, ' ')) {
      requested_labels.push_back(entry);
    }

    std::vector<std::string> scalar_node_data_labels = model_data.GetScalarNodeDataLabels();
    for (auto const & requested_label : requested_labels) {
      for (auto& node_label : scalar_node_data_labels) {
        if (requested_label == node_label) {
          int field_id = model_data.FieldId(node_label);
          int dim = model_data.GetHostScalarNodeData(field_id).extent(0);
          node_data_labels_.push_back(node_label);
          node_data_field_ids_.push_back(field_id);
          node_data_types_.push_back(FieldType::HostScalar);
          node_data_components_.push_back(0);
          node_data_.push_back(std::vector<double>(dim, 0.0));
        }
      }
    }

    std::vector<std::string> vector_node_data_labels = model_data.GetVectorNodeDataLabels();
    for (auto const & requested_label : requested_labels) {
      for (auto& node_label : vector_node_data_labels) {
        if (requested_label == node_label) {
          int field_id = model_data.FieldId(node_label);
          int dim = model_data.GetHostVectorNodeData(field_id).extent(0);
          // x component
          node_data_labels_.push_back(node_label + "_x");
          node_data_field_ids_.push_back(field_id);
          node_data_types_.push_back(FieldType::HostVector);
          node_data_components_.push_back(0);
          node_data_.push_back(std::vector<double>(dim, 0.0));
          // y component
          node_data_labels_.push_back(node_label + "_y");
          node_data_field_ids_.push_back(field_id);
          node_data_types_.push_back(FieldType::HostVector);
          node_data_components_.push_back(1);
          node_data_.push_back(std::vector<double>(dim, 0.0));
          // z component
          node_data_labels_.push_back(node_label + "_z");
          node_data_field_ids_.push_back(field_id);
          node_data_types_.push_back(FieldType::HostVector);
          node_data_components_.push_back(2);
          node_data_.push_back(std::vector<double>(dim, 0.0));
        }
      }
    }

  }

  std::vector< std::vector<double> >& ExodusOutputManager::GetNodeDataForOutput(nimble_kokkos::ModelData& model_data) {
    for (unsigned int i_data ; i_data < node_data_labels_.size() ; ++i_data) {
      int field_id = node_data_field_ids_.at(i_data);
      FieldType field_type = node_data_types_.at(i_data);
      if (field_type == FieldType::HostScalar) {
        HostScalarView data = model_data.GetHostScalarNodeData(field_id);
        for (unsigned int i=0 ; i<node_data_[i_data].size() ; i++) {
          node_data_[i_data][i] = data(i);
        }
      }
      else if (field_type == FieldType::HostVector) {
        HostVectorView data = model_data.GetHostVectorNodeData(field_id);
        int component = node_data_components_[i_data];
        for (unsigned int i=0 ; i<node_data_[i_data].size() ; i++) {
          node_data_[i_data][i] = data(i, component);
        }
      }
    }
    return node_data_;
  }

}
