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

  ElasticMaterial::ElasticMaterial(MaterialParameters const & material_parameters)
    : Material(material_parameters), num_state_variables_(0), dim_(0), density_(0.0), bulk_modulus_(0.0), shear_modulus_(0.0) {
    dim_ = 3;
    density_ = material_parameters_.GetParameterValue("density");
    bulk_modulus_ =  material_parameters_.GetParameterValue("bulk_modulus");
    shear_modulus_ =  material_parameters_.GetParameterValue("shear_modulus");
  }

  void ElasticMaterial::GetStress(int elem_id,
                                  int num_pts,
                                  double time_previous,
                                  double time_current,
                                  const double * const deformation_gradient_n,
                                  const double * const deformation_gradient_np1,
                                  const double * const stress_n,
                                  double* stress_np1,
                                  const double * const state_data_n,
                                  double* state_data_np1,
                                  DataManager& data_manager,
                                  bool is_output_step) {

    double* stress = stress_np1;
    const double * def_grad = deformation_gradient_np1;
    double strain[6];
    double trace_strain;

    double two_mu = 2.0 * shear_modulus_;
    double lambda = bulk_modulus_ - 2.0*shear_modulus_/3.0;

    for (int pt = 0 ; pt < num_pts ; pt++, def_grad+=9, stress+=6){
      strain[K_S_XX] = def_grad[K_F_XX] - 1.0;
      strain[K_S_YY] = def_grad[K_F_YY] - 1.0;
      strain[K_S_ZZ] = def_grad[K_F_ZZ] - 1.0;
      strain[K_S_XY] = 0.5*(def_grad[K_F_XY] + def_grad[K_F_YX]);
      strain[K_S_YZ] = 0.5*(def_grad[K_F_YZ] + def_grad[K_F_ZY]);
      strain[K_S_ZX] = 0.5*(def_grad[K_F_ZX] + def_grad[K_F_XZ]);

      trace_strain = strain[K_S_XX] + strain[K_S_YY] + strain[K_S_ZZ];

      stress[K_S_XX] = two_mu * strain[K_S_XX] + lambda * trace_strain;
      stress[K_S_YY] = two_mu * strain[K_S_YY] + lambda * trace_strain;
      stress[K_S_ZZ] = two_mu * strain[K_S_ZZ] + lambda * trace_strain;
      stress[K_S_XY] = two_mu * strain[K_S_XY];
      stress[K_S_YZ] = two_mu * strain[K_S_YZ];
      stress[K_S_ZX] = two_mu * strain[K_S_ZX];
    }
    // TODO rotate stress?
  }

#ifdef NIMBLE_HAVE_KOKKOS
  void ElasticMaterial::GetStress(double time_previous,
                                  double time_current,
                                  nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_n,
                                  nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_np1,
                                  nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_n,
                                  nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_np1) {

    nimble_kokkos::DeviceFullTensorIntPtSingleEntryView& def_grad = deformation_gradient_np1;
    nimble_kokkos::DeviceSymTensorIntPtSingleEntryView& stress = stress_np1;
    double strain[6];
    double trace_strain;

    double two_mu = 2.0 * shear_modulus_;
    double lambda = bulk_modulus_ - 2.0*shear_modulus_/3.0;

    strain[K_S_XX_] = def_grad[K_F_XX_] - 1.0;
    strain[K_S_YY_] = def_grad[K_F_YY_] - 1.0;
    strain[K_S_ZZ_] = def_grad[K_F_ZZ_] - 1.0;
    strain[K_S_XY_] = 0.5*(def_grad[K_F_XY_] + def_grad[K_F_YX_]);
    strain[K_S_YZ_] = 0.5*(def_grad[K_F_YZ_] + def_grad[K_F_ZY_]);
    strain[K_S_ZX_] = 0.5*(def_grad[K_F_ZX_] + def_grad[K_F_XZ_]);

    trace_strain = strain[K_S_XX] + strain[K_S_YY] + strain[K_S_ZZ];

    stress[K_S_XX_] = two_mu * strain[K_S_XX_] + lambda * trace_strain;
    stress[K_S_YY_] = two_mu * strain[K_S_YY_] + lambda * trace_strain;
    stress[K_S_ZZ_] = two_mu * strain[K_S_ZZ_] + lambda * trace_strain;
    stress[K_S_XY_] = two_mu * strain[K_S_XY_];
    stress[K_S_YZ_] = two_mu * strain[K_S_YZ_];
    stress[K_S_ZX_] = two_mu * strain[K_S_ZX_];

    // TODO rotate stress?
  }
#endif

#ifdef NIMBLE_HAVE_UQ
  void ElasticMaterial::GetOffNominalStress(const std::vector<double> & params_this_sample,
                                            const int & bulk_mod_idx,
                                            const int & shear_mod_idx,
                                            int num_pts,
                                            const double * const deformation_gradient_np1,
                                            double* stress_np1) {

    double bulk_mod = (bulk_mod_idx == -1) ? bulk_modulus_ : params_this_sample[bulk_mod_idx];
    double shear_mod = (shear_mod_idx == -1) ? shear_modulus_ : params_this_sample[shear_mod_idx];

    double* stress = stress_np1;
    const double * def_grad = deformation_gradient_np1;
    double strain[6];
    double trace_strain;

    double two_mu = 2.0 * shear_modulus_;
    double lambda = bulk_mod - 2.0*shear_mod/3.0;

    for (int pt = 0 ; pt < num_pts ; pt++, def_grad+=9, stress+=6){
      strain[K_S_XX] = def_grad[K_F_XX] - 1.0;
      strain[K_S_YY] = def_grad[K_F_YY] - 1.0;
      strain[K_S_ZZ] = def_grad[K_F_ZZ] - 1.0;
      strain[K_S_XY] = 0.5*(def_grad[K_F_XY] + def_grad[K_F_YX]);
      strain[K_S_YZ] = 0.5*(def_grad[K_F_YZ] + def_grad[K_F_ZY]);
      strain[K_S_ZX] = 0.5*(def_grad[K_F_ZX] + def_grad[K_F_XZ]);

      trace_strain = strain[K_S_XX] + strain[K_S_YY] + strain[K_S_ZZ];

      stress[K_S_XX] = two_mu * strain[K_S_XX] + lambda * trace_strain;
      stress[K_S_YY] = two_mu * strain[K_S_YY] + lambda * trace_strain;
      stress[K_S_ZZ] = two_mu * strain[K_S_ZZ] + lambda * trace_strain;
      stress[K_S_XY] = two_mu * strain[K_S_XY];
      stress[K_S_YZ] = two_mu * strain[K_S_YZ];
      stress[K_S_ZX] = two_mu * strain[K_S_ZX];
    }


  }
#endif

  void ElasticMaterial::GetTangent(int num_pts,
                                   double* material_tangent) const {

    double lambda = bulk_modulus_ - 2.0*shear_modulus_/3.0;
    double mu = shear_modulus_;
    double two_mu = 2.0*shear_modulus_;

    for (int int_pt=0 ; int_pt<num_pts ; int_pt++) {
      int offset = int_pt*36;
      material_tangent[offset]      = lambda + two_mu;
      material_tangent[offset + 1]  = lambda;
      material_tangent[offset + 2]  = lambda;
      material_tangent[offset + 3]  = 0.0;
      material_tangent[offset + 4]  = 0.0;
      material_tangent[offset + 5]  = 0.0;
      material_tangent[offset + 6]  = lambda;
      material_tangent[offset + 7]  = lambda + two_mu;
      material_tangent[offset + 8]  = lambda;
      material_tangent[offset + 9]  = 0.0;
      material_tangent[offset + 10] = 0.0;
      material_tangent[offset + 11] = 0.0;
      material_tangent[offset + 12] = lambda;
      material_tangent[offset + 13] = lambda;
      material_tangent[offset + 14] = lambda + two_mu;
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
                                     const double * const stress_n,
                                     double* stress_np1,
                                     const double * const state_data_n,
                                     double* state_data_np1,
                                     DataManager& data_manager,
                                     bool is_output_step) {

    double xj,fac,pressure,bxx,byy,bzz,bxy,byz,bzx,trace;
    double sxx,syy,szz,sxy,syz,szx,syx,szy,sxz;

    // left stretch and rotation
    double v[6], r[9];
    // Cauchy stress
    double* sig = stress_np1;

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

      sig += 6;
    }
  }

#ifdef NIMBLE_HAVE_UQ
  void NeohookeanMaterial::GetOffNominalStress(const std::vector<double> & params_this_sample,
                                               const int & bulk_mod_idx,
                                               const int & shear_mod_idx,
                                               int num_pts,
                                               const double * const deformation_gradient_np1,
                                               double* stress_np1) {

    double bulk_mod = (bulk_mod_idx == -1) ? bulk_modulus_ : params_this_sample[bulk_mod_idx];
    double shear_mod = (shear_mod_idx == -1) ? shear_modulus_ : params_this_sample[shear_mod_idx];

    double xj,fac,pressure,bxx,byy,bzz,bxy,byz,bzx,trace;
    double sxx,syy,szz,sxy,syz,szx,syx,szy,sxz;

    // left stretch and rotation
    double v[6], r[9];
    // Cauchy stress
    double* sig = stress_np1;

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

      pressure = 0.5*bulk_mod*(xj - 1.0/xj);

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

      sig[K_S_XX] = pressure + shear_mod*bxx/xj;
      sig[K_S_YY] = pressure + shear_mod*byy/xj;
      sig[K_S_ZZ] = pressure + shear_mod*bzz/xj;
      sig[K_S_XY] =            shear_mod*bxy/xj;
      sig[K_S_YZ] =            shear_mod*byz/xj;
      sig[K_S_ZX] =            shear_mod*bzx/xj;

      sig += 6;
    }

  }
#endif

#ifdef NIMBLE_HAVE_KOKKOS
  void NeohookeanMaterial::GetStress(double time_previous,
                                     double time_current,
                                     nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_n,
                                     nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_np1,
                                     nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_n,
                                     nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_np1) {

    double xj,fac,pressure,bxx,byy,bzz,bxy,byz,bzx,trace;
    double sxx,syy,szz,sxy,syz,szx,syx,szy,sxz;

    // deformation gradient, left stretch, and rotation
    double def_grad[9], v[6], r[9];
    for (int i=0 ; i<9 ; i++) {
      def_grad[i] = deformation_gradient_np1[i];
    }

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

    stress_np1(K_S_XX_) = pressure + shear_modulus_*bxx/xj;
    stress_np1(K_S_YY_) = pressure + shear_modulus_*byy/xj;
    stress_np1(K_S_ZZ_) = pressure + shear_modulus_*bzz/xj;
    stress_np1(K_S_XY_) =            shear_modulus_*bxy/xj;
    stress_np1(K_S_YZ_) =            shear_modulus_*byz/xj;
    stress_np1(K_S_ZX_) =            shear_modulus_*bzx/xj;
  }
#endif

  void NeohookeanMaterial::GetTangent(int num_pts,
                                      double* material_tangent) const {

    double lambda = bulk_modulus_ - 2.0*shear_modulus_/3.0;
    double mu = shear_modulus_;
    double two_mu = 2.0*shear_modulus_;

    for (int int_pt=0 ; int_pt<num_pts ; int_pt++) {
      int offset = int_pt*36;
      material_tangent[offset]      = lambda + two_mu;
      material_tangent[offset + 1]  = lambda;
      material_tangent[offset + 2]  = lambda;
      material_tangent[offset + 3]  = 0.0;
      material_tangent[offset + 4]  = 0.0;
      material_tangent[offset + 5]  = 0.0;
      material_tangent[offset + 6]  = lambda;
      material_tangent[offset + 7]  = lambda + two_mu;
      material_tangent[offset + 8]  = lambda;
      material_tangent[offset + 9]  = 0.0;
      material_tangent[offset + 10] = 0.0;
      material_tangent[offset + 11] = 0.0;
      material_tangent[offset + 12] = lambda;
      material_tangent[offset + 13] = lambda;
      material_tangent[offset + 14] = lambda + two_mu;
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
