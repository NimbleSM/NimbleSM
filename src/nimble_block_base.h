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

#ifndef NIMBLE_BLOCK_BASE_H
#define NIMBLE_BLOCK_BASE_H

#include <memory>

#include "nimble_genesis_mesh.h"
#include "nimble_material.h"
#include "nimble_view.h"

namespace nimble {

class Element;

class BlockBase
{
 public:
  BlockBase() = default;

  virtual ~BlockBase() = default;

  virtual void
  InstantiateElement() = 0;

  double
  GetDensity() const
  {
    return material_->GetDensity();
  }

  double
  GetBulkModulus() const
  {
    return material_->GetBulkModulus();
  }

  double
  GetShearModulus() const
  {
    return material_->GetShearModulus();
  }

  std::shared_ptr<Material>
  GetMaterialPointer() const
  {
    return material_;
  }

  std::shared_ptr<Element>
  GetElementPointer() const
  {
    return element_;
  }

  virtual double
  ComputeCriticalTimeStep(
      const nimble::Viewify<2>& node_reference_coordinates,
      const nimble::Viewify<2>& node_displacements,
      int                       num_elem,
      const int*                elem_conn) const;

  template <typename MatT>
  void
  ComputeTangentStiffnessMatrix(
      int                 num_global_unknowns,
      const double* const reference_coordinates,
      const double* const displacement,
      int                 num_elem,
      const int* const    elem_conn,
      const int* const    global_node_ids,
      MatT&               tangent_stiffness) const;

 protected:
  std::string                macro_material_parameters_ = "none";
  std::map<int, std::string> rve_material_parameters_;
  std::string                rve_boundary_condition_strategy_ = "none";
  std::vector<int>           rve_output_global_elem_ids_;
  // todo: can we avoid carrying the rve_mesh around?
  GenesisMesh rve_mesh_;

  std::shared_ptr<Element>  element_  = nullptr;
  std::shared_ptr<Material> material_ = nullptr;
};

}  // namespace nimble

#endif  // NIMBLE_BLOCK_BASE_H
