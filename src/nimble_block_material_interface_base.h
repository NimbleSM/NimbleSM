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

#ifndef SRC_NIMBLE_BLOCK_MATERIAL_INTERFACE_BASE_H
#define SRC_NIMBLE_BLOCK_MATERIAL_INTERFACE_BASE_H

namespace nimble {

class BlockBase;
class Material;

struct BlockData {
  BlockData(nimble::BlockBase *block_, nimble::Material* material_d_, const int block_id_,
            const int num_block_elems_, const int num_points_per_block_elem_)
      :
      block(block_),
      material_device(material_d_),
      id(block_id_),
      num_elems(num_block_elems_),
      num_points_per_elem(num_points_per_block_elem_) {
  }

  nimble::BlockBase *block = nullptr;
  nimble::Material *material_device = nullptr;
  int id = 0;
  int num_elems = 0;
  int num_points_per_elem = 0;
};

class BlockMaterialInterfaceBase {

public:

  BlockMaterialInterfaceBase() = default;

  virtual ~BlockMaterialInterfaceBase() = default;

  virtual void ComputeStress() const = 0;

};

}


#endif // SRC_NIMBLE_BLOCK_MATERIAL_INTERFACE_BASE_H
