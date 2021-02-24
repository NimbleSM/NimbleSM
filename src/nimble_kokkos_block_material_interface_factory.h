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

#ifndef SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_FACTORY_H_
#define SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_FACTORY_H_

#include <memory>
#include <vector>

#include "nimble_block_material_interface_factory_base.h"

namespace nimble_kokkos {

class BlockMaterialInterfaceFactory : public nimble::BlockMaterialInterfaceFactoryBase
{
 public:

  BlockMaterialInterfaceFactory() = default;
  ~BlockMaterialInterfaceFactory() override = default;

  std::shared_ptr<nimble::BlockMaterialInterfaceBase> create(double time_n, double time_np1,
                                                 const nimble::FieldIds &field_ids,
                                                 const std::vector<nimble::BlockData> &blocks,
                                                 nimble::ModelDataBase *model_data_ptr) const override;
};

}

#endif /* SRC_NIMBLE_KOKKOS_BLOCK_MATERIAL_INTERFACE_FACTORY_H_ */
