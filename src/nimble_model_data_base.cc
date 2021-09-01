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

#include "nimble_model_data_base.h"

#include "nimble_boundary_condition_manager.h"
#include "nimble_data_manager.h"

namespace nimble {

void
ModelDataBase::SetDimension(int dim)
{
  if (dim != 2 && dim != 3) throw std::invalid_argument("\nError:  Invalid dimension in ModelData\n");
  dim_ = dim;
}

void
ModelDataBase::SetReferenceCoordinates(const nimble::GenesisMesh& mesh)
{
  const double* const ref_coord_x = mesh.GetCoordinatesX();
  const double* const ref_coord_y = mesh.GetCoordinatesY();
  const double* const ref_coord_z = mesh.GetCoordinatesZ();
  //
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  int  num_nodes            = static_cast<int>(mesh.GetNumNodes());
  if (dim_ == 2) {
    for (int i = 0; i < num_nodes; i++) {
      reference_coordinate(i, 0) = ref_coord_x[i];
      reference_coordinate(i, 1) = ref_coord_y[i];
    }
  } else if (dim_ == 3) {
    for (int i = 0; i < num_nodes; i++) {
      reference_coordinate(i, 0) = ref_coord_x[i];
      reference_coordinate(i, 1) = ref_coord_y[i];
      reference_coordinate(i, 2) = ref_coord_z[i];
    }
  } else {
    throw std::invalid_argument(" -- Inappropriate Spatial Dimension");
  }
}

void
ModelDataBase::ApplyInitialConditions(nimble::DataManager& data_manager)
{
  auto bc                   = data_manager.GetBoundaryConditionManager();
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto velocity             = GetVectorNodeData("velocity");
  bc->ApplyInitialConditions(reference_coordinate, velocity);
}

void
ModelDataBase::ApplyKinematicConditions(nimble::DataManager& data_manager, double time_current, double time_previous)
{
  auto bc                   = data_manager.GetBoundaryConditionManager();
  auto reference_coordinate = GetVectorNodeData("reference_coordinate");
  auto velocity             = GetVectorNodeData("velocity");
  auto displacement         = GetVectorNodeData("displacement");
  bc->ApplyKinematicBC(time_current, time_previous, reference_coordinate, displacement, velocity);
}

}  // namespace nimble
