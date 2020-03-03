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
#include <vector>
#include <nimble_kokkos_defs.h>
#include <nimble_kokkos_material_factory.h>
#include <nimble_material.h>
#include <nimble_material_factory_util.h>
#include <stddef.h>

#ifdef NIMBLE_HAVE_EXTRAS
#include <nimble_ngp_lame_material.h>
#endif

namespace nimble_kokkos {

using nimble::Material;

template <typename MatType>
inline std::pair<std::shared_ptr<Material>, Material*> allocate_material_on_host_and_device(nimble::MaterialParameters& mat_params_struct) {
  auto mat_host = std::make_shared<MatType>(mat_params_struct);
  auto mat_device = static_cast<Material*>(Kokkos::kokkos_malloc<>("Material", sizeof(MatType)));
  nimble::Material* pointer_that_lives_on_the_stack = mat_device;
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(int) {
    new (pointer_that_lives_on_the_stack) MatType(mat_params_struct);
  });
  return std::make_pair(mat_host, mat_device);
}

void MaterialFactory::create() {
  // the first entry in the material parameters string is the material model name
  size_t space_pos = material_params.find(" ");
  std::string name = material_params.substr(0, space_pos);

#ifdef NIMBLE_HAVE_EXTRAS
  // NGP LAME material models are designated with ngp_lame_
  const bool is_ngp_lame_model = name.size() > 9 && name.substr(0,9) == "ngp_lame_";
#endif

  char material_name[nimble::MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  int num_material_parameters;
  char material_parameter_names[nimble::MaterialParameters::MAX_NUM_MAT_PARAM][nimble::MaterialParameters::MAX_MAT_MODEL_STR_LEN];
  double material_parameter_values[nimble::MaterialParameters::MAX_NUM_MAT_PARAM];
  nimble::ParseMaterialParametersString(material_params.c_str(), material_name, num_material_parameters, material_parameter_names, material_parameter_values);
  nimble::MaterialParameters material_parameters_struct(material_name, num_material_parameters, material_parameter_names, material_parameter_values, num_material_points);
  if (nimble::StringsAreEqual(material_name, "neohookean")) {
    std::tie(material_host, material_device) = allocate_material_on_host_and_device<nimble::NeohookeanMaterial>(
        material_parameters_struct);
  } else if (nimble::StringsAreEqual(material_name, "elastic")) {
    std::tie(material_host, material_device) = allocate_material_on_host_and_device<nimble::ElasticMaterial>(
        material_parameters_struct);
  }
#ifdef NIMBLE_HAVE_EXTRAS
  else if (is_ngp_lame_model) {

    lame::ngp::MaterialType ngp_lame_material_type;
    int num_state_data;
    std::vector<double> parameters_vec;
    std::vector<lame::ngp::UserFunction> mat_functions;

    if (nimble::StringsAreEqual(material_name, "ngp_lame_hypoelastic")) {
      ngp_lame_material_type = lame::ngp::MaterialType::HYPOELASTIC;
      parameters_vec.push_back(material_parameters_struct.GetParameterValue("bulk_modulus"));
      parameters_vec.push_back(material_parameters_struct.GetParameterValue("shear_modulus"));
      num_state_data = 0;
    }
    else if(nimble::StringsAreEqual(material_name, "ngp_lame_neohookean")) {
      ngp_lame_material_type = lame::ngp::MaterialType::NEOHOOKEAN;
      parameters_vec.push_back(material_parameters_struct.GetParameterValue("bulk_modulus"));
      parameters_vec.push_back(material_parameters_struct.GetParameterValue("shear_modulus"));
      num_state_data = 0;
    }
    else if(nimble::StringsAreEqual(material_name, "ngp_lame_j2_plasticity")) {
      ngp_lame_material_type = lame::ngp::MaterialType::J2_PLASTICITY;
      parameters_vec.push_back(material_parameters_struct.GetParameterValue("youngs_modulus"));
      parameters_vec.push_back(material_parameters_struct.GetParameterValue("poisson_ratio"));
      parameters_vec.push_back(material_parameters_struct.GetParameterValue("yield_stress"));
      parameters_vec.push_back(material_parameters_struct.GetParameterValue("beta"));
      parameters_vec.push_back(material_parameters_struct.GetParameterValue("hardening_modulus"));
      num_state_data = 12;
    }
    else {
      throw std::logic_error("\nError in Block::InstantiateMaterialModel(), invalid NGP LAME material model name.\n");
    }

    lame::ngp::MaterialProps mat_props(parameters_vec);
    lame::ngp::MaterialAllocator* ngp_lame_material_allocator = lame::ngp::MaterialAllocator::create(ngp_lame_material_type,
                                                                                                     parameters_vec,
                                                                                                     mat_functions);
    ngp_lame_data = std::make_shared<nimble::NGPLAMEData>(num_material_points, num_state_data);
    double dt = 1.0;
    lame::ngp::MaterialParams mat_params(ngp_lame_data->disp_grad_,
                                         ngp_lame_data->velo_grad_,
                                         ngp_lame_data->stress_old_,
                                         ngp_lame_data->stress_new_,
                                         ngp_lame_data->state_old_,
                                         ngp_lame_data->state_new_,
                                         dt);
    mtk_ngp::DevicePtr<lame::ngp::Material> ngp_lame_material_d =
        ngp_lame_material_allocator->allocate_device_material(mat_params, num_material_points);
    material_host = std::make_shared<nimble::NGPLAMEMaterial>(material_parameters_struct, ngp_lame_material_d);
  }
#endif
  else {
    throw std::logic_error("\nError in Block::InstantiateMaterialModel(), invalid material model name.\n");
  }
}

}
