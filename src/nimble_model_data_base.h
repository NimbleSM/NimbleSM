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

#include <map>
#include <string>
#include <vector>
#include <sstream>

#include "nimble_block.h"
#include "nimble_data_utils.h"
#include "nimble_exodus_output.h"
#include "nimble_genesis_mesh.h"
#include "nimble_linear_solver.h"
#include "nimble_view.h"

namespace nimble {

class DataManager;
class MaterialFactoryBase;

struct FieldIds
{
  int deformation_gradient = -1;
  int stress               = -1;
  int unrotated_stress     = -1;

  int reference_coordinates = -1;
  int displacement          = -1;
  int velocity              = -1;
  int acceleration          = -1;

  int lumped_mass    = -1;
  int internal_force = -1;
  int contact_force  = -1;
  int external_force = -1;
};

class ModelDataBase
{
 public:
  /// \brief Constructor
  ModelDataBase() = default;

  /// \brief Destructor
  virtual ~ModelDataBase() = default;

  //--- Virtual functions

  /// \brief Allocate data storage for a node-based quantity
  ///
  /// \param length
  /// \param label
  /// \param num_objects
  /// \return Field ID for the data allocated
  virtual int
  AllocateNodeData(Length length, std::string label, int num_objects) = 0;

  int
  GetFieldIdChecked(const std::string& field_label) const {
    int id = GetFieldId(field_label);
    if (id < 0) {
      std::ostringstream errc;
      errc << "Field \"" << field_label << "\" not allocated";
      throw std::runtime_error(errc.str());
    }
    return id;
  }

  /// \brief Returns the field ID for a specific label
  ///
  /// \param field_label Label for a stored quantity
  /// \return Field ID to identify the data storage
  virtual int
  GetFieldId(const std::string& field_label) const = 0;

  /// \brief Initialize the different blocks in the mesh
  ///
  /// \param data_manager Reference to the data manager
  /// \param material_factory_base Shared pointer to the material factory
  virtual void
  InitializeBlocks(
      nimble::DataManager&                        data_manager,
      const std::shared_ptr<nimble::MaterialFactoryBase>& material_factory_base) = 0;

  /// \brief Copy time state (n+1) into time state (n)
  ///
  /// \param data_manager Reference to the data manager
  virtual void
  UpdateStates(const nimble::DataManager& data_manager) = 0;

  /// \brief Get view of scalar quantity defined on nodes
  ///
  /// \param field_id
  /// \return Viewify<1> object for scalar quantity
  virtual nimble::Viewify<1>
  GetScalarNodeData(const std::string& label) = 0;

  /// \brief Get view of vector quantity defined on nodes
  ///
  /// \param field_id
  /// \return Viewify<2> object for vector quantity
  nimble::Viewify<2>
  GetVectorNodeData(const std::string& label) {
    return GetVectorNodeData(GetFieldIdChecked(label));
  }

  virtual nimble::Viewify<2>
  GetVectorNodeData(int field_id) = 0;

  /// \brief Compute the lumped mass
  ///
  /// \param data_manager Reference to the data manager
  virtual void
  ComputeLumpedMass(nimble::DataManager& data_manager) = 0;

  virtual void
  InitializeExodusOutput(nimble::DataManager& data_manager)
  {
    throw std::runtime_error(" Exodus Output Not Implemented \n");
  }

  /// \brief Write output of simulation in Exodus format
  ///
  /// \param[in] data_manager Reference to data manager
  /// \param[in] time_current Time value
  virtual void
  WriteExodusOutput(nimble::DataManager& data_manager, double time_current)
  {
    throw std::runtime_error(" Exodus Output Not Implemented \n");
  }

  /// \brief Compute the external force
  ///
  /// \param data_manager Reference to the DataManager object
  /// \param time_previous
  /// \param time_current
  /// \param is_output_step
  ///
  /// \note This routine is a placeholder.
  virtual void
  ComputeExternalForce(
      nimble::DataManager& data_manager,
      double               time_previous,
      double               time_current,
      bool                 is_output_step){};

  /// \brief Compute the internal force
  ///
  /// \param[in] data_manager
  /// \param[in] time_previous
  /// \param[in] time_current
  /// \param[in] is_output_step
  /// \param[in] displacement
  /// \param[out] internal_force  Output for internal force
  virtual void
  ComputeInternalForce(
      nimble::DataManager&      data_manager,
      double                    time_previous,
      double                    time_current,
      bool                      is_output_step,
      const nimble::Viewify<2>& displacement,
      nimble::Viewify<2>&       force){};

  /// \brief Apply initial conditions
  virtual void
  ApplyInitialConditions(nimble::DataManager& data_manager);

  /// \brief Apply kinematic conditions
  virtual void
  ApplyKinematicConditions(nimble::DataManager& data_manager, double time_current, double time_previous);

  /// \brief Update model with new velocity
  ///
  /// \param[in] data_manager Reference to the data manager
  /// \param[in] dt Current time step
  ///
  /// \note This routine is usually empty.
  ///       The UQ model data is one case using this routine.
  virtual void
  UpdateWithNewVelocity(nimble::DataManager& data_manager, double dt)
  {
  }

  /// \brief Update model with new displacement
  ///
  /// \param[in] data_manager Reference to the data manager
  /// \param[in] dt Current time step
  ///
  /// \note This routine is usually empty.
  ///       The UQ model data is one case using this routine.
  virtual void
  UpdateWithNewDisplacement(nimble::DataManager& data_manager, double dt)
  {
  }

  //--- Common interface routines

  /// \brief Get the spatial dimension
  ///
  /// \return Spatial dimension
  int
  GetDimension() const
  {
    return dim_;
  }

  /// \brief Set the critical time step
  ///
  /// \param time_step Critical time step to use.
  void
  SetCriticalTimeStep(double time_step)
  {
    critical_time_step_ = time_step;
  }

  /// \brief Set spatial dimension
  ///
  /// \param dim Spatial dimension
  void
  SetDimension(int dim);

  /// \brief Set reference coordinates
  ///
  /// \param mesh Mesh
  void
  SetReferenceCoordinates(const nimble::GenesisMesh& mesh);

  /// \brief Get the critical time step
  ///
  /// \return Time step
  double
  GetCriticalTimeStep() const
  {
    return critical_time_step_;
  }

  const std::vector<std::string>&
  GetNodeDataLabelsForOutput() const
  {
    return output_node_component_labels_;
  }

  const std::map<int, std::vector<std::string>>&
  GetElementDataLabels() const
  {
    return element_component_labels_;
  }

  const std::map<int, std::vector<std::string>>&
  GetElementDataLabelsForOutput() const
  {
    return output_element_component_labels_;
  }

  const std::map<int, std::vector<std::string>>&
  GetDerivedElementDataLabelsForOutput() const
  {
    return derived_output_element_data_labels_;
  }

  /// \brief Set the use of displacement fluctuations instead of displacement.
  ///
  void
  SetUseDisplacementFluctuations()
  {
    use_displacement_fluctuations_ = true;
  }

 protected:
  //! Problem dimension, either 2 or 3.
  int dim_ = 3;

  //! Critical time step
  double critical_time_step_ = 0.0;

  //! Output labels for node data that will be written to disk
  std::vector<std::string> output_node_component_labels_;

  //! Map key is the block_id, vector contains component-wise label for each
  //! scalar entry in the data array.
  std::map<int, std::vector<std::string>> element_component_labels_;

  //! Output labels for element data that will be written to disk.
  std::map<int, std::vector<std::string>> output_element_component_labels_;

  //! Output labels for derived element data that will be written to disk.
  std::map<int, std::vector<std::string>> derived_output_element_data_labels_;

  //! Flag to use displacement fluctuations
  bool use_displacement_fluctuations_ = false;
};

}  // namespace nimble

#endif
