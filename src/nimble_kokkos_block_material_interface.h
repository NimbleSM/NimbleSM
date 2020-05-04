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

#include <nimble_data_utils.h>
#include <nimble_kokkos_block.h>
#include <nimble_kokkos_data_manager.h>
#include <nimble_kokkos_defs.h>
#include <nimble_utils.h>

namespace nimble_kokkos
{

using ElemPointRangePolicy = Kokkos::MDRangePolicy<Kokkos::Rank<2> >;
inline ElemPointRangePolicy make_elem_point_range_policy(const int num_block_elems, const int num_points_per_elem)
{
  return ElemPointRangePolicy( { 0, 0 }, { num_block_elems, num_points_per_elem });
}

struct BlockData {
  BlockData(nimble_kokkos::Block &block_, nimble::Material* material_d_, const int block_id_,
            const int num_block_elems_, const int num_points_per_block_elem_)
      :
      block(block_),
      material_device(material_d_),
      id(block_id_),
      num_elems(num_block_elems_),
      num_points_per_elem(num_points_per_block_elem_) {
  }

  nimble_kokkos::Block &block;
  nimble::Material *material_device;
  int id;
  int num_elems;
  int num_points_per_elem;
};

class BlockMaterialInterface {
 public:
  BlockMaterialInterface(const double time_n_, const double time_np1_, const FieldIds &field_ids_,
                         const std::vector<BlockData>& blocks_,
                         nimble_kokkos::ModelData &model_data_)
      :
      time_n(time_n_),
      time_np1(time_np1_),
      field_ids(field_ids_),
      model_data(model_data_),
      blocks(blocks_) {
  }

  virtual ~BlockMaterialInterface() = default;

  virtual void ComputeStress() const;

 protected:
  const double time_n;
  const double time_np1;
  const FieldIds &field_ids;
  nimble_kokkos::ModelData &model_data;
  std::vector<BlockData> blocks;
};

}

#endif /* SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_H_ */
