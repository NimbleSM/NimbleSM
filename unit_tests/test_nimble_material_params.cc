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

#include <gtest/gtest.h>
#include <nimble_material.h>
#include <nimble_material_factory.h>
#include <nimble_material_factory_util.h>
#include <memory>
#include <string>

namespace nimble {

TEST(nimble_material_params, construct) {
  nimble::MaterialParameters params;
}

TEST(nimble_material_params, add_double_parameter) {
  nimble::MaterialParameters params;
  params.AddParameter("testParam", 1.0);

  ASSERT_TRUE(params.IsParameter("testParam"));
  ASSERT_DOUBLE_EQ(params.GetParameterValue("testParam"), 1.0);
}

class TestMaterialFactory : public MaterialFactoryBase {
 public:
  TestMaterialFactory()
      :
      MaterialFactoryBase() {
  }
  ~TestMaterialFactory() = default;

  void register_test_double_parameter(const char* name) {
    add_valid_double_parameter_name(name);
  }

  void register_test_string_parameter(const char* name) {
    add_valid_string_parameter_name(name);
  }

  std::shared_ptr<MaterialParameters> parse_string(const char *params) const {
    return ParseMaterialParametersString(params);
  }

  void create() override {}
};

TEST(nimble_material_params, invalid_parameter_throws) {
  const std::string material_string = "neohookean some_stuff 1.0e6 stuff 0.27 density 1.0e3";
  EXPECT_ANY_THROW(
      TestMaterialFactory().parse_string(material_string.c_str()));

}

TEST(nimble_material_params, parse_double_parameters) {
  const std::string material_string = "neohookean bulk_modulus 1.0e6 shear_modulus 5.e5 density 1.0e3";
  auto params = TestMaterialFactory().parse_string(material_string.c_str());

  char matName[64];
  params->GetMaterialName(matName, false);
  std::string stringMatName(matName);
  ASSERT_EQ(stringMatName, "neohookean");

  ASSERT_EQ(params->GetNumParameters(), 3);
  ASSERT_EQ(params->GetNumStringParameters(), 0);

  ASSERT_FALSE(params->IsParameter("testParam"));

  ASSERT_TRUE(params->IsParameter("bulk_modulus"));
  ASSERT_TRUE(params->IsParameter("density"));
  ASSERT_TRUE(params->IsParameter("shear_modulus"));

  EXPECT_DOUBLE_EQ(params->GetParameterValue("bulk_modulus"), 1.0e6);
  EXPECT_DOUBLE_EQ(params->GetParameterValue("shear_modulus"), 5.e5);
  EXPECT_DOUBLE_EQ(params->GetParameterValue("density"), 1.0e3);
}

TEST(nimble_material_params, register_new_test_property) {
  const std::string material_string = "neohookean bulk_modulus 1.0e6 shear_modulus 5.e5 density 1.0e3 test_property 2.0";
  TestMaterialFactory fact;
  fact.register_test_double_parameter("test_property");
  auto params = fact.parse_string(material_string.c_str());

  ASSERT_TRUE(params->IsParameter("test_property"));
  EXPECT_DOUBLE_EQ(params->GetParameterValue("test_property"), 2.0);
}

TEST(nimble_material_params, register_new_test_string_property) {
  const std::string material_string = "neohookean bulk_modulus 1.0e6 shear_modulus 5.e5 density 1.0e3 test_property custom_property_val";
  TestMaterialFactory fact;
  fact.register_test_string_parameter("test_property");
  auto params = fact.parse_string(material_string.c_str());

  ASSERT_TRUE(params->IsStringParameter("test_property"));
  auto c = params->GetStringParameterValue("test_property");
  EXPECT_EQ(std::string(params->GetStringParameterValue("test_property")), "custom_property_val");
}

}

