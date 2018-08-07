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

#include "nimble_material.h"
#include "nimble_utils.h"
#include <sstream>

namespace nimble {

    void ParseMaterialParametersString(const char* material_parameters,
                                       char material_name[MaterialParameters::MAX_MAT_MODEL_STR_LEN],
                                       int& num_material_parameters,
                                       char material_parameter_names[MaterialParameters::MAX_NUM_MAT_PARAM][MaterialParameters::MAX_MAT_MODEL_STR_LEN],
                                       double material_parameter_values[MaterialParameters::MAX_NUM_MAT_PARAM]) {

    num_material_parameters = 0;

    // The first string is the material name, followed by the material properties (key-value pairs)

    char material_parameter_value_strings[MaterialParameters::MAX_NUM_MAT_PARAM][MaterialParameters::MAX_MAT_MODEL_STR_LEN];
    size_t string_length = strlen(material_parameters);
    size_t sub_string_length = 0;
    bool last_char_was_space = false;
    bool is_material_name = true;
    bool is_param_name = false;
    bool is_param_value = false;

    memset(material_name, 0, sizeof(material_name));
    for (int i=0 ; i<MaterialParameters::MAX_NUM_MAT_PARAM ; ++i) {
      memset(material_parameter_names[i], 0, sizeof(material_parameter_names[i]));
      memset(material_parameter_value_strings[i], 0, sizeof(material_parameter_value_strings[i]));
      material_parameter_values[i] = 0.0;
    }

    for (int i=0 ; i<string_length ; ++i) {
      char val = material_parameters[i];
      if (val == ' ') {
        if (!last_char_was_space) {
          if (is_material_name) {
            is_material_name = false;
            is_param_name = true;
            sub_string_length = 0;
            num_material_parameters += 1;
          }
          else if (is_param_name) {
            is_param_name = false;
            is_param_value = true;
            sub_string_length = 0;
          }
          else if (is_param_value) {
            is_param_name = true;
            is_param_value = false;
            sub_string_length = 0;
            num_material_parameters += 1;
          }
        }
        last_char_was_space = true;
      }
      else {
        if (is_material_name) {
          material_name[sub_string_length++] = val;
        }
        else if (is_param_name) {
          material_parameter_names[num_material_parameters - 1][sub_string_length++] = val;
        }
        else if (is_param_value) {
          material_parameter_value_strings[num_material_parameters - 1][sub_string_length++] = val;
        }
        last_char_was_space = false;
      }
    }
    for (int i=0 ; i<num_material_parameters; ++i) {
      material_parameter_values[i] = atof(material_parameter_value_strings[i]);
    }
  }

  NeohookeanMaterial::NeohookeanMaterial(MaterialParameters const & material_parameters)
    : Material(material_parameters), num_state_variables_(0), dim_(0), density_(0.0), bulk_modulus_(0.0), shear_modulus_(0.0) {
    dim_ = 3;
    density_ = material_parameters_.GetParameterValue("density");
    bulk_modulus_ =  material_parameters_.GetParameterValue("bulk_modulus");
    shear_modulus_ =  material_parameters_.GetParameterValue("shear_modulus");
  }

  void NeohookeanMaterial::GetStress(int elem_id,
                                     int num_pts,
                                     double time_previous,
                                     double time_current,
                                     const double * const deformation_gradient_n,
                                     const double * const deformation_gradient_np1,
                                     const double * const unrotated_stress_n,
                                     double* unrotated_stress_np1,
                                     const double * const state_data_n,
                                     double* state_data_np1,
                                     DataManager& data_manager,
                                     bool is_output_step) {

    double xj,fac,pressure,bxx,byy,bzz,bxy,byz,bzx,trace;
    double sxx,syy,szz,sxy,syz,szx,syx,szy,sxz;

    // left stretch and rotation
    double v[6], r[9];
    // Cauchy stress
    double* sig = unrotated_stress_np1;

    for (int pt = 0 ; pt < num_pts ; pt++){

      Polar_Decomp(&deformation_gradient_np1[9*pt], v, r);

      CheckVectorSanity(9, &deformation_gradient_np1[9*pt], "neohookean deformation_gradient_np1");
      CheckVectorSanity(6, v, "neohookean v");
      CheckVectorSanity(9, r, "neohookean r");

      xj =     v[K_S_XX]*v[K_S_YY]*v[K_S_ZZ]
 	 + 2.0*v[K_S_XY]*v[K_S_YZ]*v[K_S_ZX]
	 -     v[K_S_XX]*v[K_S_YZ]*v[K_S_YZ]
	 -     v[K_S_YY]*v[K_S_ZX]*v[K_S_ZX]
	 -     v[K_S_ZZ]*v[K_S_XY]*v[K_S_XY];

      double cbrt_xj = std::cbrt(xj);
      fac = 1.0 / (cbrt_xj * cbrt_xj);

      pressure = 0.5*bulk_modulus_*(xj - 1.0/xj);

      bxx = v[K_S_XX]*v[K_S_XX]
          + v[K_S_XY]*v[K_S_YX]
          + v[K_S_XZ]*v[K_S_ZX];

      byy = v[K_S_YX]*v[K_S_XY]
          + v[K_S_YY]*v[K_S_YY]
          + v[K_S_YZ]*v[K_S_ZY];

      bzz = v[K_S_ZX]*v[K_S_XZ]
          + v[K_S_ZY]*v[K_S_YZ]
          + v[K_S_ZZ]*v[K_S_ZZ];

      bxy = v[K_S_XX]*v[K_S_XY]
          + v[K_S_XY]*v[K_S_YY]
          + v[K_S_XZ]*v[K_S_ZY];

      byz = v[K_S_YX]*v[K_S_XZ]
          + v[K_S_YY]*v[K_S_YZ]
          + v[K_S_YZ]*v[K_S_ZZ];

      bzx = v[K_S_ZX]*v[K_S_XX]
          + v[K_S_ZY]*v[K_S_YX]
          + v[K_S_ZZ]*v[K_S_ZX];

      bxx = fac*bxx;
      byy = fac*byy;
      bzz = fac*bzz;
      bxy = fac*bxy;
      byz = fac*byz;
      bzx = fac*bzx;

      trace = bxx + byy + bzz;

      bxx = bxx - trace/3.0;
      byy = byy - trace/3.0;
      bzz = bzz - trace/3.0;

      sig[K_S_XX] = pressure + shear_modulus_*bxx/xj;
      sig[K_S_YY] = pressure + shear_modulus_*byy/xj;
      sig[K_S_ZZ] = pressure + shear_modulus_*bzz/xj;
      sig[K_S_XY] =            shear_modulus_*bxy/xj;
      sig[K_S_YZ] =            shear_modulus_*byz/xj;
      sig[K_S_ZX] =            shear_modulus_*bzx/xj;

      sxx = sig[K_S_XX]*r[K_F_XX] + sig[K_S_XY]*r[K_F_YX] + sig[K_S_XZ]*r[K_F_ZX];
      syx = sig[K_S_YX]*r[K_F_XX] + sig[K_S_YY]*r[K_F_YX] + sig[K_S_YZ]*r[K_F_ZX];
      szx = sig[K_S_ZX]*r[K_F_XX] + sig[K_S_ZY]*r[K_F_YX] + sig[K_S_ZZ]*r[K_F_ZX];
      sxy = sig[K_S_XX]*r[K_F_XY] + sig[K_S_XY]*r[K_F_YY] + sig[K_S_XZ]*r[K_F_ZY];
      syy = sig[K_S_YX]*r[K_F_XY] + sig[K_S_YY]*r[K_F_YY] + sig[K_S_YZ]*r[K_F_ZY];
      szy = sig[K_S_ZX]*r[K_F_XY] + sig[K_S_ZY]*r[K_F_YY] + sig[K_S_ZZ]*r[K_F_ZY];
      sxz = sig[K_S_XX]*r[K_F_XZ] + sig[K_S_XY]*r[K_F_YZ] + sig[K_S_XZ]*r[K_F_ZZ];
      syz = sig[K_S_YX]*r[K_F_XZ] + sig[K_S_YY]*r[K_F_YZ] + sig[K_S_YZ]*r[K_F_ZZ];
      szz = sig[K_S_ZX]*r[K_F_XZ] + sig[K_S_ZY]*r[K_F_YZ] + sig[K_S_ZZ]*r[K_F_ZZ];

      sig[K_S_XX] = r[K_F_XX]*sxx + r[K_F_YX]*syx + r[K_F_ZX]*szx;
      sig[K_S_YY] = r[K_F_XY]*sxy + r[K_F_YY]*syy + r[K_F_ZY]*szy;
      sig[K_S_ZZ] = r[K_F_XZ]*sxz + r[K_F_YZ]*syz + r[K_F_ZZ]*szz;
      sig[K_S_XY] = r[K_F_XX]*sxy + r[K_F_YX]*syy + r[K_F_ZX]*szy;
      sig[K_S_YZ] = r[K_F_XY]*sxz + r[K_F_YY]*syz + r[K_F_ZY]*szz;
      sig[K_S_ZX] = r[K_F_XZ]*sxx + r[K_F_YZ]*syx + r[K_F_ZZ]*szx;

      sig += 6;
    }
  }

#ifdef NIMBLE_HAVE_KOKKOS
  void NeohookeanMaterial::GetStress(double time_previous,
                                     double time_current,
                                     nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_n,
                                     nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_np1,
                                     nimble_kokkos::DeviceSymTensorIntPtSingleEntryView unrotated_stress_n,
                                     nimble_kokkos::DeviceSymTensorIntPtSingleEntryView unrotated_stress_np1) {

    double xj,fac,pressure,bxx,byy,bzz,bxy,byz,bzx,trace;
    double sxx,syy,szz,sxy,syz,szx,syx,szy,sxz;

    // deformation gradient, left stretch, and rotation
    double def_grad[9], v[6], r[9];
    for (int i=0 ; i<9 ; i++) {
      def_grad[i] = deformation_gradient_np1[i];
    }
    // Cauchy stress
    nimble_kokkos::DeviceSymTensorIntPtSingleEntryView& sig = unrotated_stress_np1;

    Polar_Decomp(def_grad, v, r);

    CheckVectorSanity(9, def_grad, "neohookean deformation_gradient_np1");
    CheckVectorSanity(6, v, "neohookean v");
    CheckVectorSanity(9, r, "neohookean r");

    xj =     v[K_S_XX_]*v[K_S_YY_]*v[K_S_ZZ_]
       + 2.0*v[K_S_XY_]*v[K_S_YZ_]*v[K_S_ZX_]
       -     v[K_S_XX_]*v[K_S_YZ_]*v[K_S_YZ_]
       -     v[K_S_YY_]*v[K_S_ZX_]*v[K_S_ZX_]
       -     v[K_S_ZZ_]*v[K_S_XY_]*v[K_S_XY_];

    double cbrt_xj = std::cbrt(xj);
    fac = 1.0 / (cbrt_xj * cbrt_xj);

    pressure = 0.5*bulk_modulus_*(xj - 1.0/xj);

    bxx = v[K_S_XX_]*v[K_S_XX_]
        + v[K_S_XY_]*v[K_S_YX_]
        + v[K_S_XZ_]*v[K_S_ZX_];

    byy = v[K_S_YX_]*v[K_S_XY_]
        + v[K_S_YY_]*v[K_S_YY_]
        + v[K_S_YZ_]*v[K_S_ZY_];

    bzz = v[K_S_ZX_]*v[K_S_XZ_]
        + v[K_S_ZY_]*v[K_S_YZ_]
        + v[K_S_ZZ_]*v[K_S_ZZ_];

    bxy = v[K_S_XX_]*v[K_S_XY_]
        + v[K_S_XY_]*v[K_S_YY_]
        + v[K_S_XZ_]*v[K_S_ZY_];

    byz = v[K_S_YX_]*v[K_S_XZ_]
        + v[K_S_YY_]*v[K_S_YZ_]
        + v[K_S_YZ_]*v[K_S_ZZ_];

    bzx = v[K_S_ZX_]*v[K_S_XX_]
        + v[K_S_ZY_]*v[K_S_YX_]
        + v[K_S_ZZ_]*v[K_S_ZX_];

    bxx = fac*bxx;
    byy = fac*byy;
    bzz = fac*bzz;
    bxy = fac*bxy;
    byz = fac*byz;
    bzx = fac*bzx;

    trace = bxx + byy + bzz;

    bxx = bxx - trace/3.0;
    byy = byy - trace/3.0;
    bzz = bzz - trace/3.0;

    sig(K_S_XX_) = pressure + shear_modulus_*bxx/xj;
    sig(K_S_YY_) = pressure + shear_modulus_*byy/xj;
    sig(K_S_ZZ_) = pressure + shear_modulus_*bzz/xj;
    sig(K_S_XY_) =            shear_modulus_*bxy/xj;
    sig(K_S_YZ_) =            shear_modulus_*byz/xj;
    sig(K_S_ZX_) =            shear_modulus_*bzx/xj;

    sxx = sig(K_S_XX_)*r[K_F_XX_] + sig(K_S_XY_)*r[K_F_YX_] + sig(K_S_XZ_)*r[K_F_ZX_];
    syx = sig(K_S_YX_)*r[K_F_XX_] + sig(K_S_YY_)*r[K_F_YX_] + sig(K_S_YZ_)*r[K_F_ZX_];
    szx = sig(K_S_ZX_)*r[K_F_XX_] + sig(K_S_ZY_)*r[K_F_YX_] + sig(K_S_ZZ_)*r[K_F_ZX_];
    sxy = sig(K_S_XX_)*r[K_F_XY_] + sig(K_S_XY_)*r[K_F_YY_] + sig(K_S_XZ_)*r[K_F_ZY_];
    syy = sig(K_S_YX_)*r[K_F_XY_] + sig(K_S_YY_)*r[K_F_YY_] + sig(K_S_YZ_)*r[K_F_ZY_];
    szy = sig(K_S_ZX_)*r[K_F_XY_] + sig(K_S_ZY_)*r[K_F_YY_] + sig(K_S_ZZ_)*r[K_F_ZY_];
    sxz = sig(K_S_XX_)*r[K_F_XZ_] + sig(K_S_XY_)*r[K_F_YZ_] + sig(K_S_XZ_)*r[K_F_ZZ_];
    syz = sig(K_S_YX_)*r[K_F_XZ_] + sig(K_S_YY_)*r[K_F_YZ_] + sig(K_S_YZ_)*r[K_F_ZZ_];
    szz = sig(K_S_ZX_)*r[K_F_XZ_] + sig(K_S_ZY_)*r[K_F_YZ_] + sig(K_S_ZZ_)*r[K_F_ZZ_];

    sig(K_S_XX_) = r[K_F_XX_]*sxx + r[K_F_YX_]*syx + r[K_F_ZX_]*szx;
    sig(K_S_YY_) = r[K_F_XY_]*sxy + r[K_F_YY_]*syy + r[K_F_ZY_]*szy;
    sig(K_S_ZZ_) = r[K_F_XZ_]*sxz + r[K_F_YZ_]*syz + r[K_F_ZZ_]*szz;
    sig(K_S_XY_) = r[K_F_XX_]*sxy + r[K_F_YX_]*syy + r[K_F_ZX_]*szy;
    sig(K_S_YZ_) = r[K_F_XY_]*sxz + r[K_F_YY_]*syz + r[K_F_ZY_]*szz;
    sig(K_S_ZX_) = r[K_F_XZ_]*sxx + r[K_F_YZ_]*syx + r[K_F_ZZ_]*szx;

  }
#endif

  void NeohookeanMaterial::GetTangent(int num_pts,
                                      double* material_tangent) const {

    double lambda = bulk_modulus_ - 2.0*shear_modulus_/3.0;
    double mu = shear_modulus_;

    for (int int_pt=0 ; int_pt<num_pts ; int_pt++) {
      int offset = int_pt*36;
      material_tangent[offset]      = lambda + 2.0*mu;
      material_tangent[offset + 1]  = lambda;
      material_tangent[offset + 2]  = lambda;
      material_tangent[offset + 3]  = 0.0;
      material_tangent[offset + 4]  = 0.0;
      material_tangent[offset + 5]  = 0.0;
      material_tangent[offset + 6]  = lambda;
      material_tangent[offset + 7]  = lambda + 2.0*mu;
      material_tangent[offset + 8]  = lambda;
      material_tangent[offset + 9]  = 0.0;
      material_tangent[offset + 10] = 0.0;
      material_tangent[offset + 11] = 0.0;
      material_tangent[offset + 12] = lambda;
      material_tangent[offset + 13] = lambda;
      material_tangent[offset + 14] = lambda + 2.0*mu;
      material_tangent[offset + 15] = 0.0;
      material_tangent[offset + 16] = 0.0;
      material_tangent[offset + 17] = 0.0;
      material_tangent[offset + 18] = 0.0;
      material_tangent[offset + 19] = 0.0;
      material_tangent[offset + 20] = 0.0;
      material_tangent[offset + 21] = mu;
      material_tangent[offset + 22] = 0.0;
      material_tangent[offset + 23] = 0.0;
      material_tangent[offset + 24] = 0.0;
      material_tangent[offset + 25] = 0.0;
      material_tangent[offset + 26] = 0.0;
      material_tangent[offset + 27] = 0.0;
      material_tangent[offset + 28] = mu;
      material_tangent[offset + 29] = 0.0;
      material_tangent[offset + 30] = 0.0;
      material_tangent[offset + 31] = 0.0;
      material_tangent[offset + 32] = 0.0;
      material_tangent[offset + 33] = 0.0;
      material_tangent[offset + 34] = 0.0;
      material_tangent[offset + 35] = mu;
    }
  }
}
