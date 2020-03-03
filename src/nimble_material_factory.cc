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

#ifdef NIMBLE_HAVE_EXTRAS
#include <nimble_lame_material.h>
#endif

namespace nimble {

void MaterialFactory::create() {
  // the first entry in the material parameters string is the material model name
  size_t space_pos = material_params.find(" ");
  const std::string name = material_params.substr(0, space_pos);

#ifdef NIMBLE_HAVE_EXTRAS
  // LAME material models are designated with lame_
  const bool is_lame_model = name.size() > 5 && name.substr(0,5) == "lame_";
#endif

  char material_name[MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  int num_material_parameters;
  char material_parameter_names[MaterialParameters::MAX_NUM_MAT_PARAM][MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  double material_parameter_values[MaterialParameters::MAX_NUM_MAT_PARAM];
  ParseMaterialParametersString(material_params.c_str(), material_name, num_material_parameters, material_parameter_names, material_parameter_values);
  MaterialParameters material_parameters_struct(material_name, num_material_parameters, material_parameter_names, material_parameter_values);

  if (StringsAreEqual(material_name, "neohookean")) {
   material = std::make_shared<NeohookeanMaterial>(material_parameters_struct);
  }
  else if (StringsAreEqual(material_name, "elastic")) {
    material = std::make_shared<ElasticMaterial>(material_parameters_struct);
  }
#ifdef NIMBLE_HAVE_EXTRAS
  else if (is_lame_model) {
    material = std::make_shared<LAMEMaterial>(material_parameters_struct);
  }
#endif
  else {
    throw std::logic_error("\nError in Block::InstantiateMaterialModel(), invalid material model name.\n");
  }

}

}
