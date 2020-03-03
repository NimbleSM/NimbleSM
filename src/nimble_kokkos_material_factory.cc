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

#include <stdexcept>
#include <tuple>
#include <utility>
#include <nimble_kokkos_defs.h>
#include <nimble_kokkos_material_factory.h>
#include <nimble_material_factory_util.h>

namespace nimble_kokkos {

using nimble::Material;

MaterialFactory::MaterialFactory()
    : material_device(nullptr) {
}

void MaterialFactory::parse_and_create(const std::string& mat_params, const int num_points) {
  material_params = nimble::ParseMaterialParametersString(mat_params.c_str(), num_points);
  create();
}

template <typename MatType>
inline std::pair<std::shared_ptr<Material>, Material*> allocate_material_on_host_and_device(const nimble::MaterialParameters& mat_params_struct) {
  auto mat_host = std::make_shared<MatType>(mat_params_struct);
  auto mat_device = static_cast<Material*>(Kokkos::kokkos_malloc<>("Material", sizeof(MatType)));
  nimble::Material* pointer_that_lives_on_the_stack = mat_device;
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
    new (pointer_that_lives_on_the_stack) MatType(mat_params_struct);
  });
  return std::make_pair(mat_host, mat_device);
}

void MaterialFactory::create() {
  char name[nimble::MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  material_params->GetMaterialName(name, false);
  std::string name_string(name);
  if (nimble::StringsAreEqual(name_string.c_str(), "neohookean")) {
    std::tie(material_host, material_device) = allocate_material_on_host_and_device<nimble::NeohookeanMaterial>(
        *material_params);
  } else if (nimble::StringsAreEqual(name_string.c_str(), "elastic")) {
    std::tie(material_host, material_device) = allocate_material_on_host_and_device<nimble::ElasticMaterial>(
        *material_params);
  } else {
    throw std::logic_error("\nError in Block::InstantiateMaterialModel(), invalid material model name.\n");
  }
}

}
