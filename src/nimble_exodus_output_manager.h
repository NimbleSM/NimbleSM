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

#ifndef NIMBLE_EXODUS_OUTPUT_MANAGER_H
#define NIMBLE_EXODUS_OUTPUT_MANAGER_H

#include "nimble_kokkos_data_manager.h"

#ifdef NIMBLE_HAVE_DARMA
  #include "darma.h"
#else
  #include <string>
#endif

namespace nimble_kokkos {

  class ExodusOutputManager {

  public:

    ExodusOutputManager() {};

    void SpecifyOutputFields(nimble_kokkos::ModelData& model_data,
                             std::string const & output_command_string);

    std::vector<std::string> GetNodeDataLabelsForOutput() { return node_data_labels_; }

    std::vector< std::vector<double> > GetNodeDataForOutput(nimble_kokkos::ModelData& model_data);

    std::map<int, std::vector<std::string> > GetElementDataLabelsForOutput() { return elem_data_labels_; }

    std::map<int, std::vector< std::vector<double> > > GetElementDataForOutput(nimble_kokkos::ModelData& model_data);

    std::map<int, std::vector<int> > GetSymTensorFieldIdsRequiringVolumeAverage() { return sym_tensor_field_ids_requiring_volume_average_; }

    std::map<int, std::vector<int> > GetFullTensorFieldIdsRequiringVolumeAverage() { return full_tensor_field_ids_requiring_volume_average_; }

  private:

    std::vector<std::string> node_data_labels_;
    std::vector<int> node_data_field_ids_;
    std::vector<FieldType> node_data_types_;
    std::vector<int> node_data_components_;
    std::vector< std::vector<double> > node_data_;

    std::map<int, std::vector<std::string> > elem_data_labels_;
    std::map<int, std::vector<int> > elem_data_field_ids_;
    std::map<int, std::vector<FieldType> > elem_data_types_;
    std::map<int, std::vector<int> > elem_data_components_;
    std::map<int, std::vector< std::vector<double> > > elem_data_;

    std::map<int, std::vector<int> > sym_tensor_field_ids_requiring_volume_average_;
    std::map<int, std::vector<int> > full_tensor_field_ids_requiring_volume_average_;
  };

} // namespace nimble

#endif // NIMBLE_OUTPUT_EXODUS_H
