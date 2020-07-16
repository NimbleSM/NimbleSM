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
#include <stddef.h>

namespace nimble {

MaterialFactoryBase::MaterialFactoryBase()
{
  NeohookeanMaterial::register_supported_material_parameters(*this);
  ElasticMaterial::register_supported_material_parameters(*this);
}

std::shared_ptr<MaterialParameters> MaterialFactoryBase::ParseMaterialParametersString(const char *material_parameters,
                                                                                       const int num_material_points) const {

  char material_name[MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  int num_material_parameters = 0;
  char all_material_parameter_names[MaterialParameters::MAX_NUM_MAT_PARAM][MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  char material_double_parameter_names[MaterialParameters::MAX_NUM_MAT_PARAM][MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  double material_double_parameter_values[MaterialParameters::MAX_NUM_MAT_PARAM];
  char material_string_parameter_names[MaterialParameters::MAX_NUM_MAT_PARAM][MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  char material_string_parameter_values[MaterialParameters::MAX_NUM_MAT_PARAM][MaterialParameters::MAX_MAT_MODEL_STR_LEN];

  // The first string is the material name, followed by the material properties (key-value pairs)

  char all_material_parameter_value_strings[MaterialParameters::MAX_NUM_MAT_PARAM][MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  size_t string_length = strlen(material_parameters);
  size_t sub_string_length = 0;
  bool last_char_was_space = false;
  bool is_material_name = true;
  bool is_param_name = false;
  bool is_param_value = false;

  memset(material_name, 0, MaterialParameters::MAX_MAT_MODEL_STR_LEN);
  for (int i = 0; i < MaterialParameters::MAX_NUM_MAT_PARAM; ++i) {
    memset(all_material_parameter_value_strings[i], 0, MaterialParameters::MAX_MAT_MODEL_STR_LEN);
    memset(all_material_parameter_names[i], 0, MaterialParameters::MAX_MAT_MODEL_STR_LEN);
    memset(material_double_parameter_names[i], 0, MaterialParameters::MAX_MAT_MODEL_STR_LEN);
    material_double_parameter_values[i] = 0.0;
    memset(material_string_parameter_names[i], 0, MaterialParameters::MAX_MAT_MODEL_STR_LEN);
    memset(material_string_parameter_values[i], 0, MaterialParameters::MAX_MAT_MODEL_STR_LEN);
  }

  for (int i = 0; i < string_length; ++i) {
    char val = material_parameters[i];
    if (val == ' ') {
      if (!last_char_was_space) {
        if (is_material_name) {
          is_material_name = false;
          is_param_name = true;
          sub_string_length = 0;
          num_material_parameters += 1;
        } else if (is_param_name) {
          is_param_name = false;
          is_param_value = true;
          sub_string_length = 0;
        } else if (is_param_value) {
          is_param_name = true;
          is_param_value = false;
          sub_string_length = 0;
          num_material_parameters += 1;
        }
      }
      last_char_was_space = true;
    } else {
      if (is_material_name) {
        material_name[sub_string_length++] = val;
      } else if (is_param_name) {
        all_material_parameter_names[num_material_parameters - 1][sub_string_length++] = val;
      } else if (is_param_value) {
        all_material_parameter_value_strings[num_material_parameters - 1][sub_string_length++] = val;
      }
      last_char_was_space = false;
    }
  }

  int num_material_double_parameters = 0;
  int num_material_string_parameters = 0;
  for (int i = 0; i < num_material_parameters; ++i) {
    int name_len = StringLength(all_material_parameter_names[i]);
    if (std::find(valid_double_parameter_names.begin(), valid_double_parameter_names.end(),
                  all_material_parameter_names[i]) != valid_double_parameter_names.end()) {
      for (int j = 0; j < name_len; ++j) {
        material_double_parameter_names[num_material_double_parameters][j] = all_material_parameter_names[i][j];
      }
      material_double_parameter_values[num_material_double_parameters] = atof(all_material_parameter_value_strings[i]);
      ++num_material_double_parameters;
    } else if (std::find(valid_string_parameter_names.begin(), valid_string_parameter_names.end(),
                         all_material_parameter_names[i]) != valid_string_parameter_names.end()) {
      int len = StringLength(all_material_parameter_value_strings[i]);
      for (int j = 0; j < name_len; ++j) {
        material_string_parameter_names[num_material_string_parameters][j] = all_material_parameter_names[i][j];
      }
      for (int j = 0; j < len; ++j) {
        material_string_parameter_values[num_material_string_parameters][j] = all_material_parameter_value_strings[i][j];
      }
      ++num_material_string_parameters;
    } else {
      std::string errMsg = "Invalid material parameter encountered: '"
          + std::string(all_material_parameter_names[i]) + "'";
      std::runtime_error err(errMsg);
      throw err;
    }
  }

  assert(num_material_double_parameters + num_material_string_parameters == num_material_parameters);

  return std::make_shared<MaterialParameters>(material_name, num_material_double_parameters,
                                              material_double_parameter_names, material_double_parameter_values,
                                              num_material_string_parameters, material_string_parameter_names,
                                              material_string_parameter_values, num_material_points);
}

MaterialFactory::MaterialFactory()
    :
    MaterialFactoryBase() {
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
