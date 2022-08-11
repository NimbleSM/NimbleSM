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

#include "nimble_material_factory_base.h"

#include <stdexcept>
#include <tuple>
#include <vector>

#include "nimble_macros.h"
#include "nimble_material.h"
#include "nimble_parser_util.h"

namespace nimble {

MaterialFactoryBase::MaterialFactoryBase()
{
    NeohookeanMaterial::register_supported_material_parameters(*this);
    ElasticMaterial::register_supported_material_parameters(*this);
    J2PlasticityMaterial::register_supported_material_parameters(*this);
}

std::shared_ptr<MaterialParameters>
MaterialFactoryBase::ParseMaterialParametersString(const std::string& material_parameters, int num_material_points)
    const
{
  auto tokens = nimble::tokenize_string(material_parameters);

  // The first string is the material name, followed by the material properties
  // (key-value pairs)

  NIMBLE_ASSERT(tokens.size() > 1);

  const std::string material_name = tokens.front();
  auto              token         = tokens.cbegin() + 1;
  auto              tokens_end    = tokens.cend();

  std::map<std::string, std::string> material_string_parameters;
  std::map<std::string, double>      material_double_parameters;
  for (; token != tokens_end; token += 2) {
    auto&& key = *token;
    auto&& val = *(token + 1);
    if (std::find(valid_double_parameter_names.begin(), valid_double_parameter_names.end(), key) !=
        valid_double_parameter_names.end()) {
      double double_val = nimble::string_to_double(val);
      material_double_parameters.insert(std::make_pair(key, double_val));
    } else if (
        std::find(valid_string_parameter_names.begin(), valid_string_parameter_names.end(), key) !=
        valid_string_parameter_names.end()) {
      material_string_parameters.insert(std::make_pair(key, val));
    } else {
      std::string errMsg = "Invalid material parameter encountered: '" + key + "'";
      throw std::invalid_argument(errMsg);
    }
  }

  return std::make_shared<MaterialParameters>(
      material_name, material_string_parameters, material_double_parameters, num_material_points);
}

std::map<std::string, double>
MaterialFactoryBase::ParseMaterialParamsStringToMap(const std::string& material_parameters) const
{
  auto tokens = nimble::tokenize_string(material_parameters);

  // The first string is the material name, followed by the material properties (key-value pairs)
  NIMBLE_ASSERT(tokens.size() > 1);

  const std::string material_name = tokens.front();
  auto              token         = tokens.cbegin() + 1;
  auto              tokens_end    = tokens.cend();

  std::map<std::string, std::string> material_string_parameters;
  std::map<std::string, double>      material_double_parameters;
  for (; token != tokens_end; token += 2) {
    auto&& key = *token;
    auto&& val = *(token + 1);
    if (std::find(valid_double_parameter_names.begin(), valid_double_parameter_names.end(), key) !=
        valid_double_parameter_names.end()) {
      double double_val = nimble::string_to_double(val);
      material_double_parameters.insert(std::make_pair(key, double_val));
    } else if (
        std::find(valid_string_parameter_names.begin(), valid_string_parameter_names.end(), key) !=
        valid_string_parameter_names.end()) {
      material_string_parameters.insert(std::make_pair(key, val));
    } else {
      std::string errMsg = "Invalid material parameter encountered: '" + key + "'";
      throw std::invalid_argument(errMsg);
    }
  }

  return material_double_parameters;
}

}  // namespace nimble
