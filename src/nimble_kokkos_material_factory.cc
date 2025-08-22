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

#include <nimble_kokkos_defs.h>
#include <nimble_kokkos_material_factory.h>
#include <nimble_material.h>

#include <stdexcept>
#include <tuple>
#include <utility>

namespace nimble_kokkos {

using nimble::Material;

MaterialFactory::MaterialFactory() : MaterialFactoryBase(), material_device(nullptr) {}

template <typename MatType>
inline std::pair<std::shared_ptr<MatType>, MatType*>
allocate_material_on_host_and_device(const nimble::MaterialParameters& mat_params_struct)
{
  auto     mat_host     = std::make_shared<MatType>(mat_params_struct);
  auto     mat_device   = static_cast<MatType*>(Kokkos::kokkos_malloc<>("Material", sizeof(MatType)));
  MatType& mat_host_ref = *mat_host;
  Kokkos::parallel_for(
      1, KOKKOS_LAMBDA(int) { new (mat_device) MatType(mat_host_ref); });
  Kokkos::fence();
  return std::make_pair(mat_host, mat_device);
}

void
MaterialFactory::create()
{
  auto name_string = material_params->GetMaterialName(false);
  if (name_string == "neohookean") {
    std::tie(material, material_device) =
        allocate_material_on_host_and_device<nimble::NeohookeanMaterial>(*material_params);
  } else if (name_string == "elastic") {
    std::tie(material, material_device) =
        allocate_material_on_host_and_device<nimble::ElasticMaterial>(*material_params);
  } else if (name_string == "j2-plasticity") {
    std::tie(material, material_device) =
        allocate_material_on_host_and_device<nimble::J2PlasticityMaterial>(*material_params);
  } else {
    throw std::invalid_argument(
        "\nError in Block::InstantiateMaterialModel(), invalid material model "
        "name.\n");
  }
}

}  // namespace nimble_kokkos
