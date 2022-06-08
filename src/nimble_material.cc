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

#include <nimble_material.h>
#include <nimble_material_factory.h>
#include <p3a_exp.hpp>

#include <cmath>

namespace nimble {

void
ElasticMaterial::register_supported_material_parameters(MaterialFactoryBase& factory)
{
  factory.add_valid_double_parameter_name("bulk_modulus");
  factory.add_valid_double_parameter_name("shear_modulus");
  factory.add_valid_double_parameter_name("density");
}

ElasticMaterial::ElasticMaterial(MaterialParameters const& material_parameters)
    : Material(), num_state_variables_(0), dim_(0), density_(0.0), bulk_modulus_(0.0), shear_modulus_(0.0)
{
  dim_           = 3;
  density_       = material_parameters.GetParameterValue("density");
  bulk_modulus_  = material_parameters.GetParameterValue("bulk_modulus");
  shear_modulus_ = material_parameters.GetParameterValue("shear_modulus");
}

void
ElasticMaterial::GetStress(
    int           elem_id,
    int           num_pts,
    double        time_previous,
    double        time_current,
    const double* deformation_gradient_n,
    const double* deformation_gradient_np1,
    const double* stress_n,
    double*       stress_np1,
    const double* state_data_n,
    double*       state_data_np1,
    DataManager&  data_manager,
    bool          is_output_step)
{
  double*       stress   = stress_np1;
  const double* def_grad = deformation_gradient_np1;

  nimble::Viewify<1, const double> null_view;

  for (int pt = 0; pt < num_pts; pt++, def_grad += 9, stress += 6) {
    nimble::Viewify<1, const double> def_grap_np1(def_grad, 9);
    GetStress(time_previous, time_current, null_view, def_grap_np1, null_view, {stress, 6});
  }
  // TODO rotate stress?
}

void
ElasticMaterial::GetStress(
    double                            time_previous,
    double                            time_current,
    nimble::Viewify<1, const double>& deformation_gradient_n,
    nimble::Viewify<1, const double>& deformation_gradient_np1,
    nimble::Viewify<1, const double>& stress_n,
    nimble::Viewify<1>                stress_np1) const
{
  auto   def_grad = deformation_gradient_np1.data();
  double strain[6];
  double trace_strain;

  double two_mu = 2.0 * shear_modulus_;
  double lambda = bulk_modulus_ - 2.0 * shear_modulus_ / 3.0;

  strain[K_S_XX] = def_grad[K_F_XX] - 1.0;
  strain[K_S_YY] = def_grad[K_F_YY] - 1.0;
  strain[K_S_ZZ] = def_grad[K_F_ZZ] - 1.0;
  strain[K_S_XY] = 0.5 * (def_grad[K_F_XY] + def_grad[K_F_YX]);
  strain[K_S_YZ] = 0.5 * (def_grad[K_F_YZ] + def_grad[K_F_ZY]);
  strain[K_S_ZX] = 0.5 * (def_grad[K_F_ZX] + def_grad[K_F_XZ]);

  trace_strain = strain[K_S_XX] + strain[K_S_YY] + strain[K_S_ZZ];

  stress_np1(K_S_XX) = two_mu * strain[K_S_XX] + lambda * trace_strain;
  stress_np1(K_S_YY) = two_mu * strain[K_S_YY] + lambda * trace_strain;
  stress_np1(K_S_ZZ) = two_mu * strain[K_S_ZZ] + lambda * trace_strain;
  stress_np1(K_S_XY) = two_mu * strain[K_S_XY];
  stress_np1(K_S_YZ) = two_mu * strain[K_S_YZ];
  stress_np1(K_S_ZX) = two_mu * strain[K_S_ZX];
}

#ifdef NIMBLE_HAVE_UQ
void
ElasticMaterial::GetOffNominalStress(
    const double&       bulk_mod,
    const double&       shear_mod,
    int                 num_pts,
    const double* const deformation_gradient_np1,
    double*             stress_np1)
{
  //--- Copy current values for bulk and shear moduli
  double old_bulk  = bulk_modulus_;
  double old_shear = shear_modulus_;

  bulk_modulus_  = bulk_mod;
  shear_modulus_ = shear_mod;

  // Cauchy stress
  double*                          sig = stress_np1;
  nimble::Viewify<1, const double> null_view;

  for (int pt = 0; pt < num_pts; pt++) {
    nimble::Viewify<1, const double> def_grad_view(&deformation_gradient_np1[9 * pt], 9);
    GetStress(0.0, 0.0, null_view, def_grad_view, null_view, {sig, 6});
    sig += 6;
  }

  //--- Reset bulk and shear moduli
  bulk_modulus_  = old_bulk;
  shear_modulus_ = old_shear;
}
#endif

void
ElasticMaterial::GetTangent(int num_pts, double* material_tangent) const
{
  double lambda = bulk_modulus_ - 2.0 * shear_modulus_ / 3.0;
  double mu     = shear_modulus_;
  double two_mu = 2.0 * shear_modulus_;

  for (int int_pt = 0; int_pt < num_pts; int_pt++) {
    int offset                    = int_pt * 36;
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

void
NeohookeanMaterial::register_supported_material_parameters(MaterialFactoryBase& factory)
{
  factory.add_valid_double_parameter_name("bulk_modulus");
  factory.add_valid_double_parameter_name("shear_modulus");
  factory.add_valid_double_parameter_name("density");
}

NeohookeanMaterial::NeohookeanMaterial(MaterialParameters const& material_parameters)
    : Material(),
      num_state_variables_(0),
      dim_(3),
      density_(material_parameters.GetParameterValue("density")),
      bulk_modulus_(material_parameters.GetParameterValue("bulk_modulus")),
      shear_modulus_(material_parameters.GetParameterValue("shear_modulus"))
{
}

void
NeohookeanMaterial::GetStress(
    int                 elem_id,
    int                 num_pts,
    double              time_previous,
    double              time_current,
    const double* const deformation_gradient_n,
    const double* const deformation_gradient_np1,
    const double* const stress_n,
    double*             stress_np1,
    const double* const state_data_n,
    double*             state_data_np1,
    DataManager&        data_manager,
    bool                is_output_step)
{
  // Cauchy stress
  double* sig = stress_np1;

  nimble::Viewify<1, const double> null_view;
  for (int pt = 0; pt < num_pts; pt++) {
    nimble::Viewify<1, const double> def_grad(&deformation_gradient_np1[9 * pt], 9);
    GetStress(time_previous, time_current, null_view, def_grad, null_view, {sig, 6});
    sig += 6;
  }
}

void
NeohookeanMaterial::GetStress(
    double                            time_previous,
    double                            time_current,
    nimble::Viewify<1, const double>& deformation_gradient_n,
    nimble::Viewify<1, const double>& deformation_gradient_np1,
    nimble::Viewify<1, const double>& stress_n,
    nimble::Viewify<1>                stress_np1) const
{
  double xj, fac, pressure, bxx, byy, bzz, bxy, byz, bzx, trace;

  // left stretch and rotation
  double v[6], r[9];
  Polar_Decomp(deformation_gradient_np1.data(), v, r);

  CheckVectorSanity(9, deformation_gradient_np1.data(), "neohookean deformation_gradient_np1");
  CheckVectorSanity(6, v, "neohookean v");
  CheckVectorSanity(9, r, "neohookean r");

  xj = v[K_S_XX] * v[K_S_YY] * v[K_S_ZZ] + 2.0 * v[K_S_XY] * v[K_S_YZ] * v[K_S_ZX] - v[K_S_XX] * v[K_S_YZ] * v[K_S_YZ] -
       v[K_S_YY] * v[K_S_ZX] * v[K_S_ZX] - v[K_S_ZZ] * v[K_S_XY] * v[K_S_XY];

  double cbrt_xj = std::cbrt(xj);
  fac            = 1.0 / (cbrt_xj * cbrt_xj);

  pressure = 0.5 * bulk_modulus_ * (xj - 1.0 / xj);

  bxx = v[K_S_XX] * v[K_S_XX] + v[K_S_XY] * v[K_S_YX] + v[K_S_XZ] * v[K_S_ZX];

  byy = v[K_S_YX] * v[K_S_XY] + v[K_S_YY] * v[K_S_YY] + v[K_S_YZ] * v[K_S_ZY];

  bzz = v[K_S_ZX] * v[K_S_XZ] + v[K_S_ZY] * v[K_S_YZ] + v[K_S_ZZ] * v[K_S_ZZ];

  bxy = v[K_S_XX] * v[K_S_XY] + v[K_S_XY] * v[K_S_YY] + v[K_S_XZ] * v[K_S_ZY];

  byz = v[K_S_YX] * v[K_S_XZ] + v[K_S_YY] * v[K_S_YZ] + v[K_S_YZ] * v[K_S_ZZ];

  bzx = v[K_S_ZX] * v[K_S_XX] + v[K_S_ZY] * v[K_S_YX] + v[K_S_ZZ] * v[K_S_ZX];

  bxx = fac * bxx;
  byy = fac * byy;
  bzz = fac * bzz;
  bxy = fac * bxy;
  byz = fac * byz;
  bzx = fac * bzx;

  trace = bxx + byy + bzz;

  bxx = bxx - trace / 3.0;
  byy = byy - trace / 3.0;
  bzz = bzz - trace / 3.0;

  stress_np1(K_S_XX) = pressure + shear_modulus_ * bxx / xj;
  stress_np1(K_S_YY) = pressure + shear_modulus_ * byy / xj;
  stress_np1(K_S_ZZ) = pressure + shear_modulus_ * bzz / xj;
  stress_np1(K_S_XY) = shear_modulus_ * bxy / xj;
  stress_np1(K_S_YZ) = shear_modulus_ * byz / xj;
  stress_np1(K_S_ZX) = shear_modulus_ * bzx / xj;
}

#ifdef NIMBLE_HAVE_UQ
void
NeohookeanMaterial::GetOffNominalStress(
    const double& bulk_mod,
    const double& shear_mod,
    int           num_pts,
    const double* deformation_gradient_np1,
    double*       stress_np1)
{
  //--- Copy current values for bulk and shear moduli
  double old_bulk  = bulk_modulus_;
  double old_shear = shear_modulus_;

  bulk_modulus_  = bulk_mod;
  shear_modulus_ = shear_mod;

  // Cauchy stress
  double*                          sig = stress_np1;
  nimble::Viewify<1, const double> null_view;

  for (int pt = 0; pt < num_pts; pt++) {
    nimble::Viewify<1, const double> def_grad_view(&deformation_gradient_np1[9 * pt], 9);
    GetStress(0.0, 0.0, null_view, def_grad_view, null_view, {sig, 6});
    sig += 6;
  }

  //--- Reset bulk and shear moduli
  bulk_modulus_  = old_bulk;
  shear_modulus_ = old_shear;
}
#endif

void
NeohookeanMaterial::GetTangent(int num_pts, double* material_tangent) const
{
  double lambda = bulk_modulus_ - 2.0 * shear_modulus_ / 3.0;
  double mu     = shear_modulus_;
  double two_mu = 2.0 * shear_modulus_;

  for (int int_pt = 0; int_pt < num_pts; int_pt++) {
    int offset                    = int_pt * 36;
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

void
J2PlasticityMaterial::register_supported_material_parameters(MaterialFactoryBase& factory)
{
  factory.add_valid_double_parameter_name("bulk_modulus");
  factory.add_valid_double_parameter_name("shear_modulus");
  factory.add_valid_double_parameter_name("density");
  factory.add_valid_double_parameter_name("yield_stress");
  factory.add_valid_double_parameter_name("hardening_exponent");
  factory.add_valid_double_parameter_name("reference_plastic_strain");
  factory.add_valid_double_parameter_name("reference_viscoplastic_stress");
  factory.add_valid_double_parameter_name("rate_dependence_exponent");
  factory.add_valid_double_parameter_name("reference_plastic_strain_rate");
  factory.add_valid_double_parameter_name("melting_temperature");
  factory.add_valid_double_parameter_name("reference_temperature");
  factory.add_valid_double_parameter_name("temperature_exponent");
  factory.add_valid_double_parameter_name("specific_heat");
  factory.add_valid_double_parameter_name("taylor_quinney_coefficient");
}

J2PlasticityMaterial::J2PlasticityMaterial(MaterialParameters const& material_parameters)
    : Material(),
      num_state_variables_(10),
      dim_(3)
{
  props_.rho0 = material_parameters.GetParameterValue("density");
  props_.kappa = material_parameters.GetParameterValue("bulk_modulus");
  props_.mu = material_parameters.GetParameterValue("shear_modulus");
}

void
J2PlasticityMaterial::GetStress(
  int           elem_id,
  int           num_pts,
  double        time_previous,
  double        time_current,
  const double* deformation_gradient_n,
  const double* deformation_gradient_np1,
  const double* stress_n,
  double*       stress_np1,
  const double* state_data_n,
  double*       state_data_np1,
  DataManager&  data_manager,
  bool          is_output_step)
{
}

void
J2PlasticityMaterial::GetTangent(int num_pts, double* material_tangent) const
{
}

#ifdef NIMBLE_HAVE_UQ
void
J2PlasticityMaterial::GetOffNominalStress(
  const double& bulk_mod,
  const double& shear_mod,
  int           num_pts,
  const double* deformation_gradient_np1,
  double*       stress_np1)
{
}
#endif

void
J2PlasticityMaterial::GetStress(
  double                            time_previous,
  double                            time_current,
  nimble::Viewify<1, const double>& deformation_gradient_n,
  nimble::Viewify<1, const double>& deformation_gradient_np1,
  nimble::Viewify<1, const double>& stress_n,
  nimble::Viewify<1>                stress_np1) const
{ 
}

}  // namespace nimble
