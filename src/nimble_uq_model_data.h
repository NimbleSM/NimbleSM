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

#ifndef NIMBLESM_NIMBLE_UQ_MODEL_DATA_H
#define NIMBLESM_NIMBLE_UQ_MODEL_DATA_H

#ifdef NIMBLE_HAVE_UQ

#include "nimble_data_manager.h"
#include "nimble_model_data.h"
#include "nimble_view.h"

namespace nimble_uq {

class UqModel;

class ModelData : public nimble::ModelData
{
 public:
  ModelData() = default;

  ~ModelData() override = default;

  /// \brief Initialize the different blocks in the mesh
  ///
  /// \param data_manager Reference to the data manager
  /// \param material_factory_base Shared pointer to the material factory
  void
  InitializeBlocks(nimble::DataManager& data_manager, const std::shared_ptr<MaterialFactoryType>& material_factory_base)
      override;

  /// \brief Write output of simulation in Exodus format
  ///
  /// \param[in] data_manager Reference to data manager
  /// \param[in] time_current Time value
  void
  WriteExodusOutput(nimble::DataManager& data_manager, double time_current) override;

  /// \brief Apply initial conditions
  void
  ApplyInitialConditions(nimble::DataManager& data_manager) override;

  /// \brief Apply kinematic conditions
  void
  ApplyKinematicConditions(DataManager& data_manager, double time_current, double time_previous) override;

  /// \brief Update model with new velocity
  ///
  /// \param[in] data_manager Reference to the data manager
  /// \param[in] dt Current time step
  ///
  /// \note This routine is usually empty.
  ///       The UQ model data is one case using this routine.
  void
  UpdateWithNewVelocity(nimble::DataManager& data_manager, double dt) override;

  /// \brief Update model with new displacement
  ///
  /// \param[in] data_manager Reference to the data manager
  /// \param[in] dt Current time step
  ///
  /// \note This routine is usually empty.
  ///       The UQ model data is one case using this routine.
  void
  UpdateWithNewDisplacement(nimble::DataManager& data_manager, double dt) override;

 protected:
  std::shared_ptr<nimble::UqModel> uq_model_;

  std::vector<nimble::Viewify<2>> bc_offnom_velocity_views_;
};

}  // namespace nimble_uq

#endif

#endif  // NIMBLESM_NIMBLE_UQ_MODEL_DATA_H
