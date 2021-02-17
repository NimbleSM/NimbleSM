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

#ifndef NIMBLE_MODEL_DATA_BASE_H
#define NIMBLE_MODEL_DATA_BASE_H

#include "nimble_block.h"
#include "nimble_data_utils.h"
#include "nimble_exodus_output.h"
#include "nimble_genesis_mesh.h"
#include "nimble_linear_solver.h"

#include <vector>
#include <string>
#include <map>

namespace nimble {


struct FieldIds {
  int deformation_gradient = -1;
  int stress = -1;
  int unrotated_stress = -1;

  int reference_coordinates = -1;
  int displacement = -1;
  int velocity = -1;
  int acceleration = -1;

  int lumped_mass = -1;
  int internal_force = -1;
  int contact_force = -1;
};


class BaseModelData {

public:

  BaseModelData() = default;

  virtual ~BaseModelData() = default;

  void SetDimension(int dim) {
    if(dim != 2 && dim != 3){
      throw std::logic_error("\nError:  Invalid dimension in ModelData\n");
    }
    dim_ = dim;
  }

  int GetDimension() const {return dim_; }

  void SetCriticalTimeStep(double time_step) { critical_time_step_ = time_step; }

  double GetCriticalTimeStep() const { return critical_time_step_; }

  const std::vector<std::string> & GetNodeDataLabelsForOutput() const {
    return output_node_component_labels_;
  }

  const std::map<int, std::vector<std::string> > & GetElementDataLabels() const {
    return element_component_labels_;
  }

  const std::map<int, std::vector<std::string> > & GetElementDataLabelsForOutput() const {
    return output_element_component_labels_;
  }

  const std::map<int, std::vector<std::string> > & GetDerivedElementDataLabelsForOutput() const {
    return derived_output_element_data_labels_;
  }

  void SetDerivedElementDataLabelsForOutput(std::map<int, std::vector<std::string> > &&ref)
  {
    derived_output_element_data_labels_ = ref;
  }

  virtual int AllocateNodeData(Length length,
                               std::string label,
                               int num_objects) = 0;

  virtual int GetFieldId(const std::string& field_label) const = 0;

  virtual void SetReferenceCoordinates(const nimble::GenesisMesh &mesh) = 0;

protected:

  //! Problem dimension, either 2 or 3.
  int dim_ = 3;

  //! Critical time step
  double critical_time_step_ = 0.0;

  //! Output labels for node data that will be written to disk
  std::vector<std::string> output_node_component_labels_;

  //! Map key is the block_id, vector contains component-wise label for each scalar entry in the data array.
  std::map<int, std::vector<std::string>> element_component_labels_;

  //! Output labels for element data that will be written to disk.
  std::map<int, std::vector<std::string> > output_element_component_labels_;

  //! Output labels for derived element data that will be written to disk.
  std::map<int, std::vector<std::string> > derived_output_element_data_labels_;

};

} // namespace nimble

#endif
