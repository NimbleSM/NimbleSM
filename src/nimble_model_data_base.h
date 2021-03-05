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

class DataManager;

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


class ModelDataBase {

public:

  /// \brief Constructor
  ModelDataBase() = default;

  /// \brief Destructor
  virtual ~ModelDataBase() = default;

  /// \brief Set the spatial dimension
  ///
  /// \param dim Dimension
  void SetDimension(int dim) {
    if(dim != 2 && dim != 3){
      throw std::logic_error("\nError:  Invalid dimension in ModelData\n");
    }
    dim_ = dim;
  }

  //--- Virtual functions

  /// \brief Allocate data storage for a node-based quantity
  ///
  /// \param length
  /// \param label
  /// \param num_objects
  /// \return Field ID for the data allocated
  virtual int AllocateNodeData(Length length,
                               std::string label,
                               int num_objects) = 0;

  /// \brief Returns the field ID for a specific label
  ///
  /// \param field_label Label for a stored quantity
  /// \return Field ID to identify the data storage
  virtual int GetFieldId(const std::string& field_label) const = 0;

  /// \brief Set the reference coordinates
  ///
  /// \param mesh Reference to the global mesh
  virtual void SetReferenceCoordinates(const nimble::GenesisMesh &mesh) = 0;

  /// \brief Initialize the different blocks in the mesh
  ///
  /// \param data_manager Reference to the data manager
  /// \param material_factory_base Shared pointer to the material factory
  virtual void InitializeBlocks(nimble::DataManager &data_manager,
                                const std::shared_ptr<MaterialFactoryType> &material_factory_base) = 0;

  /// \brief Swap states between time n and time (n+1)
  ///
  /// \param data_manager Reference to the data manager
  virtual void SwapStates(const nimble::DataManager &data_manager) = 0;

  //--- Common interface routines

  /// \brief Get the spatial dimension
  ///
  /// \return Spatial dimension
  int GetDimension() const {return dim_; }

  /// \brief Set the critical time step
  ///
  /// \param time_step Critical time step to use.
  void SetCriticalTimeStep(double time_step) { critical_time_step_ = time_step; }

  /// \brief Get the critical time step
  ///
  /// \return Time step
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
