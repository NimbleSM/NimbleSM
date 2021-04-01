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

#ifndef SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_H_
#define SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_H_

#include "nimble_block_material_interface_base.h"
#include "nimble_data_utils.h"
#include "nimble_kokkos_block.h"
#include "nimble_kokkos_defs.h"
#include "nimble_utils.h"

namespace nimble {

class ModelDataBase;
class FieldIds;

}  // namespace nimble

namespace nimble_kokkos {

class ModelData;

using ElemPointRangePolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
inline ElemPointRangePolicy
make_elem_point_range_policy(const int num_block_elems, const int num_points_per_elem)
{
  return ElemPointRangePolicy({0, 0}, {num_block_elems, num_points_per_elem});
}

class BlockMaterialInterface : public nimble::BlockMaterialInterfaceBase
{
 public:
  BlockMaterialInterface(
      double                                time_n_,
      double                                time_np1_,
      const nimble::FieldIds&               field_ids_,
      const std::vector<nimble::BlockData>& blocks_,
      nimble::ModelDataBase*                model_data_);

  ~BlockMaterialInterface() override = default;

  void
  ComputeStress() const override;

 protected:
  const double                   time_n;
  const double                   time_np1;
  const nimble::FieldIds&        field_ids;
  nimble_kokkos::ModelData*      model_data;
  std::vector<nimble::BlockData> blocks;
};

}  // namespace nimble_kokkos

#endif /* SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_H_ */
