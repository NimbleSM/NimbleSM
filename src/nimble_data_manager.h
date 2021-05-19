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

#include "nimble_block_material_interface_factory_base.h"
#include "nimble_data_utils.h"
#include "nimble_defs.h"
#include "nimble_exodus_output.h"
#include "nimble_model_data.h"

#ifdef NIMBLE_HAVE_DARMA
#include "darma.h"
#else
#include <map>
#include <string>
#include <vector>
#endif

namespace nimble {

class BoundaryConditionManager;
class GenesisMesh;
class ModelDataBase;
class Parser;
class VectorCommunicator;

struct RVEData
{
  nimble::ModelData                               model_data_;
  std::vector<double>                             residual_vector_;
  std::vector<double>                             linear_solver_solution_;
  CRSMatrixContainer                              tangent_stiffness_;
  bool                                            write_exodus_output_;
  ExodusOutput                                    exodus_output_;
  std::map<int, std::vector<std::vector<double>>> derived_elem_data_;

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    ar | model_data_ | residual_vector_ | linear_solver_solution_ | tangent_stiffness_ | write_exodus_output_ |
        exodus_output_ | derived_elem_data_;
  }
#endif
};

class DataManager
{
 public:
  /// \brief Constructor
  ///
  /// \param parser Reference to parser information
  /// \param mesh Reference to mesh
  /// \param rve_mesh Reference to RVE mesh
  ///
  /// \note The RVE mesh could be empty.
  DataManager(const nimble::Parser& parser, const nimble::GenesisMesh& mesh, const nimble::GenesisMesh& rve_mesh);

  /// \brief Destructor
  ~DataManager() = default;

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    ar | macroscale_data_ | rve_data_;
  }
#endif

  /// \brief Allocate RVE data for an integration point in an element
  ///
  /// \param global_element_id Global ID to element
  /// \param integration_point_id Integration point ID
  /// \return Reference to RVE data
  RVEData&
  AllocateRVEData(int global_element_id, int integration_point_id);

  /// \brief Get RVE data for an integration point in an element
  ///
  /// \param global_element_id Global ID to element
  /// \param integration_point_id Integration point ID
  /// \return Reference to RVE data
  RVEData&
  GetRVEData(int global_element_id, int integration_point_id);

  /// \brief Initialize the data for Exodus output
  ///
  /// \param filename File name for the output files
  void
  InitializeOutput(const std::string& filename);

  /// \brief Return constant reference to parser information
  ///
  /// \return Reference to parser information
  const nimble::Parser&
  GetParser() const
  {
    return parser_;
  }

  /// \brief Return constant reference to mesh
  ///
  /// \return Reference to mesh
  const nimble::GenesisMesh&
  GetMesh() const
  {
    return mesh_;
  }

  /// \brief Return constant reference to RVE mesh
  ///
  /// \return Reference to RVE mesh
  const nimble::GenesisMesh&
  GetRVEMesh() const
  {
    return rve_mesh_;
  }

  /// \brief Return shared pointer to ModelData objet
  ///
  /// \return Shared pointer
  std::shared_ptr<nimble::ModelDataBase>
  GetMacroScaleData()
  {
    return macroscale_data_;
  }

  /// \brief Return a const reference to the field IDs
  ///
  /// \return Field IDs
  const nimble::FieldIds&
  GetFieldIDs() const
  {
    return field_ids_;
  }

  /// \brief Return reference to the field IDs
  ///
  /// \return Field IDs
  nimble::FieldIds&
  GetFieldIDs()
  {
    return field_ids_;
  }

  /// \brief Return shared pointer to VectorCommunicator objet
  ///
  /// \return Shared pointer
  std::shared_ptr<nimble::VectorCommunicator>
  GetVectorCommunicator()
  {
    return vector_communicator_;
  }

  /// \brief Return reference to RVE-macroscale deformation gradient
  ///
  /// \return Reference
  std::vector<double>&
  GetRVEDeformationGradient()
  {
    return rve_macroscale_deformation_gradient_;
  }

  /// \brief Write output of simulations
  ///
  /// \param[in] time_current Time value
  void
  WriteOutput(double time_current);

  std::shared_ptr<nimble::ExodusOutput>
  GetExodusOutput()
  {
    return exodus_output_;
  }

  /// \brief Set BlockMaterialInterfaceFactoryBase object and initialize
  ///        block data information.
  ///
  void
  SetBlockMaterialInterfaceFactory(
      const std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase>& block_material_factory);

  /// \brief Return shared pointer to BlockMaterialInterfaceFactoryBase object
  ///
  /// \return Shared pointer
  const std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase>&
  GetBlockMaterialInterfaceFactory() const;

  /// \brief Return shared pointer to BoundaryConditionManager object
  ///
  /// \return Shared pointer
  std::shared_ptr<nimble::BoundaryConditionManager>
  GetBoundaryConditionManager()
  {
    return boundary_condition_;
  }

 protected:
  /// \brief Initialize data for simulation
  void
  Initialize();

 protected:
  const nimble::Parser&                  parser_;
  const nimble::GenesisMesh&             mesh_;
  const nimble::GenesisMesh&             rve_mesh_;
  std::shared_ptr<nimble::ModelDataBase> macroscale_data_ = nullptr;

  std::map<std::pair<int, int>, RVEData> rve_data_;
  std::vector<double>                    rve_macroscale_deformation_gradient_;

  nimble::FieldIds field_ids_;

  std::shared_ptr<nimble::VectorCommunicator> vector_communicator_ = nullptr;
  std::shared_ptr<nimble::ExodusOutput>       exodus_output_       = nullptr;

  std::shared_ptr<nimble::BlockMaterialInterfaceFactoryBase> block_material_factory_ = nullptr;

  std::shared_ptr<nimble::BoundaryConditionManager> boundary_condition_ = nullptr;
};

}  // namespace nimble

#endif
