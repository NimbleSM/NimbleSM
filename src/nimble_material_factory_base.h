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

#ifndef SRC_NIMBLE_MATERIAL_FACTORY_BASE_H_
#define SRC_NIMBLE_MATERIAL_FACTORY_BASE_H_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace nimble {
class Material;
class MaterialParameters;
}

namespace nimble {

class MaterialFactoryBase {
 private:
  static inline void find_or_insert_string_in_vector(const std::string &str, std::vector<std::string> &vec) {
    if (std::find(vec.begin(), vec.end(), str) == vec.end()) {
      vec.push_back(str);
    }
  }

 public:
  MaterialFactoryBase();
  virtual ~MaterialFactoryBase() = default;

  inline void add_valid_double_parameter_name(const char *name) {
    find_or_insert_string_in_vector(std::string(name), valid_double_parameter_names);
  }

  inline void add_valid_string_parameter_name(const char *name) {
    find_or_insert_string_in_vector(std::string(name), valid_string_parameter_names);
  }

  virtual std::shared_ptr<nimble::Material> get_material() const
  { return material; }

  virtual void parse_and_create(const std::string& mat_params,
                                int num_points)
  {
    material_params = ParseMaterialParametersString(mat_params, num_points);
    create();
  }

  virtual void parse_and_create(const std::string& mat_params)
  { parse_and_create(mat_params, 0); }

protected:

  virtual void create() = 0;

protected:

  std::shared_ptr<nimble::MaterialParameters> ParseMaterialParametersString(const std::string& material_parameters,
                                                                            int num_material_points = 0) const;

  //
  //--- Protected Variables
  //

  std::shared_ptr< nimble::Material > material = nullptr;
  std::shared_ptr<const nimble::MaterialParameters> material_params;

private:
  std::vector<std::string> valid_double_parameter_names;
  std::vector<std::string> valid_string_parameter_names;
};

}

#endif /* SRC_NIMBLE_MATERIAL_FACTORY_BASE_H_ */
