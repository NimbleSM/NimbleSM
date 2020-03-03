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

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <nimble_material.h>
#include <nimble_material_factory.h>
#include <nimble_material_factory_util.h>
#include <stddef.h>

namespace nimble {

MaterialFactory::MaterialFactory() {
}

void MaterialFactory::parse_and_create(const std::string &mat_params) {
  material_params = ParseMaterialParametersString(mat_params.c_str());
  create();
}

void MaterialFactory::create() {
  char name[nimble::MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  material_params->GetMaterialName(name, false);
  std::string name_string(name);
  if (StringsAreEqual(name_string.c_str(), "neohookean")) {
    material = std::make_shared<NeohookeanMaterial>(*material_params);
  } else if (StringsAreEqual(name_string.c_str(), "elastic")) {
    material = std::make_shared<ElasticMaterial>(*material_params);
  } else {
    throw std::logic_error("\nError in Block::InstantiateMaterialModel(), invalid material model name.\n");
  }
}

}
