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
#include <p3a_log.hpp>

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

using Tensor = p3a::matrix3x3<double>;

void
variational_J2_update(
    j2::Properties const& props,
    double const          dt,
    Tensor const&         F,
    Tensor&               Fp,
    Tensor&               sigma,
    double&               eqps)
{
  auto const K    = props.kappa;
  auto const G    = props.mu;
  auto const J    = p3a::determinant(F);
  auto const Jm13 = 1.0 / std::cbrt(J);
  auto const Jm23 = Jm13 * Jm13;
  auto const logJ = std::log(J);
  auto const p    = K * logJ / J;

  auto       Fe_tr        = F * p3a::inverse(Fp);
  auto       dev_Ce_tr    = Jm23 * p3a::transpose(Fe_tr) * Fe_tr;
  auto       dev_Ee_tr    = 0.5 * p3a::log(dev_Ce_tr);
  auto const dev_M_tr     = 2.0 * G * dev_Ee_tr;
  auto const sigma_tr_eff = std::sqrt(1.5) * p3a::norm(dev_M_tr);
  auto       Np           = Tensor(0, 0, 0, 0, 0, 0, 0, 0, 0);
  if (sigma_tr_eff > 0) {
    Np = 1.5 * dev_M_tr / sigma_tr_eff;
  }

  auto       S0 = j2::FlowStrength(props, eqps);
  auto const r0 = sigma_tr_eff - S0;
  auto       r  = r0;

  auto       delta_eqps         = 0.0;
  auto const residual_tolerance = 1e-10;
  auto const deqps_tolerance    = 1e-10;
  if (r > residual_tolerance) {
    constexpr auto max_iters = 8;
    auto           iters     = 0;
    auto           merit_old = 1.0;
    auto           merit_new = 1.0;

    auto converged = false;
    while (!converged) {
      if (iters == max_iters) break;
      auto ls_is_finished = false;
      auto delta_eqps0    = delta_eqps;
      merit_old           = r * r;
      auto H              = j2::HardeningRate(props, eqps + delta_eqps) + j2::ViscoplasticHardeningRate(props, delta_eqps, dt);
      auto dr             = -3.0 * G - H;
      auto correction     = -r / dr;

      // line search
      auto       alpha                  = 1.0;
      auto       line_search_iterations = 0;
      auto const backtrack_factor       = 0.1;
      auto const decrease_factor        = 1e-5;
      while (!ls_is_finished) {
        if (line_search_iterations == 20) {
          // line search has failed to satisfactorily improve newton step
          // just take the full newton step and hope for the best
          alpha = 1;
          break;
        }
        ++line_search_iterations;
        delta_eqps = delta_eqps0 + alpha * correction;
        if (delta_eqps < 0) delta_eqps = 0;
        auto Yeq          = j2::FlowStrength(props, eqps + delta_eqps);
        auto Yvis         = j2::ViscoplasticStress(props, delta_eqps, dt);
        auto residual     = sigma_tr_eff - 3.0 * G * delta_eqps - (Yeq + Yvis);
        merit_new         = residual * residual;
        auto decrease_tol = 1.0 - 2.0 * alpha * decrease_factor;
        if (merit_new <= decrease_tol * merit_old) {
          merit_old      = merit_new;
          ls_is_finished = true;
        } else {
          auto alpha_old = alpha;
          alpha          = alpha_old * alpha_old * merit_old / (merit_new - merit_old + 2.0 * alpha_old * merit_old);
          if (backtrack_factor * alpha_old > alpha) {
            alpha = backtrack_factor * alpha_old;
          }
        }
      }
      auto S    = j2::FlowStrength(props, eqps + delta_eqps) + j2::ViscoplasticStress(props, delta_eqps, dt);
      r         = sigma_tr_eff - 3.0 * G * delta_eqps - S;
      converged = (std::abs(r / r0) < residual_tolerance) || (delta_eqps < deqps_tolerance);
      ++iters;
    }
    if (!converged) {
      printf("variational J2 did not converge to specified tolerance 1.0e-10\n");
      // TODO: handle non-convergence error
    }
    auto dFp = p3a::exp(delta_eqps * Np);
    Fp       = dFp * Fp;
    eqps += delta_eqps;
  }
  auto const Ee_correction = delta_eqps * Np;
  auto const sigma_dev     = 1.0 / J * p3a::transpose(p3a::inverse(Fe_tr)) * (dev_M_tr - 2.0 * G * Ee_correction) * p3a::transpose(Fe_tr);
  auto const sigma_vol     = p * Tensor(1, 0, 0, 0, 1, 0, 0, 0, 1);
  sigma = sigma_dev + sigma_vol;
}

void
J2PlasticityMaterial::register_supported_material_parameters(MaterialFactoryBase& factory)
{
  factory.add_valid_double_parameter_name("elastic_modulus");
  factory.add_valid_double_parameter_name("poissons_ratio");
  factory.add_valid_double_parameter_name("density");
  factory.add_valid_double_parameter_name("yield_stress");
  factory.add_valid_double_parameter_name("hardening_exponent");
  factory.add_valid_double_parameter_name("reference_plastic_strain");
  factory.add_valid_double_parameter_name("reference_viscoplastic_stress");
  factory.add_valid_double_parameter_name("rate_dependence_exponent");
  factory.add_valid_double_parameter_name("reference_plastic_strain_rate");
}

J2PlasticityMaterial::J2PlasticityMaterial(MaterialParameters const& material_parameters)
    : Material()
{
  props_.E = material_parameters.GetParameterValue("elastic_modulus");
  props_.nu = material_parameters.GetParameterValue("poissons_ratio");
  props_.rho0 = material_parameters.GetParameterValue("density");
  props_.Y0 = material_parameters.GetParameterValue("yield_stress");
  props_.n = material_parameters.GetParameterValue("hardening_exponent");
  props_.eps0 = material_parameters.GetParameterValue("reference_plastic_strain");
  props_.Svis0 = material_parameters.GetParameterValue("reference_viscoplastic_stress");
  props_.m = material_parameters.GetParameterValue("rate_dependence_exponent");
  props_.eps_dot0 = material_parameters.GetParameterValue("reference_plastic_strain_rate");
  props_.kappa = props_.E / (1.0 - 2.0 * props_.nu) / 3.0;
  props_.mu  = props_.E / (1.0 + props_.nu) / 2.0;
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
  auto const num_ivs = NumStateVariables();
  auto const EQPS_OFFSET = num_ivs - 1;
  auto const defgrad_dim = dim_ * dim_;
  auto const sigma_dim = 2 * dim_;
  for (auto point = 0; point < num_pts; ++point) {
    Tensor F(
      deformation_gradient_np1[defgrad_dim * point + K_F_XX],
      deformation_gradient_np1[defgrad_dim * point + K_F_XY],
      deformation_gradient_np1[defgrad_dim * point + K_F_XZ],
      deformation_gradient_np1[defgrad_dim * point + K_F_YX],
      deformation_gradient_np1[defgrad_dim * point + K_F_YY],
      deformation_gradient_np1[defgrad_dim * point + K_F_YZ],
      deformation_gradient_np1[defgrad_dim * point + K_F_ZX],
      deformation_gradient_np1[defgrad_dim * point + K_F_ZY],
      deformation_gradient_np1[defgrad_dim * point + K_F_ZZ]
    );
    Tensor Fp(
      state_data_n[num_ivs * point + K_F_XX],
      state_data_n[num_ivs * point + K_F_XY],
      state_data_n[num_ivs * point + K_F_XZ],
      state_data_n[num_ivs * point + K_F_YX],
      state_data_n[num_ivs * point + K_F_YY],
      state_data_n[num_ivs * point + K_F_YZ],
      state_data_n[num_ivs * point + K_F_ZX],
      state_data_n[num_ivs * point + K_F_ZY],
      state_data_n[num_ivs * point + K_F_ZZ]
    );

    double eqps = state_data_n[num_ivs * point + EQPS_OFFSET];

    Tensor sigma(0, 0, 0, 0, 0, 0, 0, 0, 0);

    auto const dt = time_current - time_previous;

    variational_J2_update(props_, dt, F, Fp, sigma, eqps);

    stress_np1[sigma_dim * point + K_S_XX] = sigma(0, 0);
    stress_np1[sigma_dim * point + K_S_YY] = sigma(1, 1);
    stress_np1[sigma_dim * point + K_S_ZZ] = sigma(2, 2);
    stress_np1[sigma_dim * point + K_S_XY] = sigma(0, 1);
    stress_np1[sigma_dim * point + K_S_YZ] = sigma(1, 2);
    stress_np1[sigma_dim * point + K_S_ZX] = sigma(2, 0);

    state_data_np1[num_ivs * point + K_F_XX] = Fp(0, 0);
    state_data_np1[num_ivs * point + K_F_XY] = Fp(0, 1);
    state_data_np1[num_ivs * point + K_F_XZ] = Fp(0, 2);
    state_data_np1[num_ivs * point + K_F_YX] = Fp(1, 0);
    state_data_np1[num_ivs * point + K_F_YY] = Fp(1, 1);
    state_data_np1[num_ivs * point + K_F_YZ] = Fp(1, 2);
    state_data_np1[num_ivs * point + K_F_ZX] = Fp(2, 0);
    state_data_np1[num_ivs * point + K_F_ZY] = Fp(2, 1);
    state_data_np1[num_ivs * point + K_F_ZZ] = Fp(2, 2);

    state_data_np1[num_ivs * point + EQPS_OFFSET] = eqps;
  }
}

void
J2PlasticityMaterial::GetTangent(int num_pts, double* material_tangent) const
{
  NIMBLE_ABORT("Tangent for J2 plasticity not implemented.");
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
