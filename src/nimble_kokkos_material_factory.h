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

#ifndef SRC_NIMBLE_KOKKOS_MATERIAL_FACTORY_H_
#define SRC_NIMBLE_KOKKOS_MATERIAL_FACTORY_H_

#include <memory>
#include <string>

#include "nimble_material_factory_base.h"

namespace nimble {
struct NGPLAMEData;
}

namespace nimble_kokkos {

class MaterialFactory : public nimble::MaterialFactoryBase
{
 public:
  explicit MaterialFactory();
  ~MaterialFactory() override = default;

  //--- Specific routines to Kokkos-based environment

  inline std::shared_ptr<nimble::Material>
  get_material_host() const
  {
    return material;
  }

  inline nimble::Material*
  get_material_device() const
  {
    return material_device;
  }

  inline std::shared_ptr<nimble::NGPLAMEData>
  get_ngp_lame_data() const
  {
    return ngp_lame_data;
  }

 protected:
  void
  create() override;

 protected:
  ///
  /// \note The variable material is used as 'material_host'
  /// (i.e. a variable defined on the host memory space).
  ///

  nimble::Material*                    material_device;
  std::shared_ptr<nimble::NGPLAMEData> ngp_lame_data;
};

}  // namespace nimble_kokkos

#endif /* SRC_NIMBLE_KOKKOS_MATERIAL_FACTORY_H_ */