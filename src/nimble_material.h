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

#ifndef NIMBLE_MATERIAL_H
#define NIMBLE_MATERIAL_H

#include <algorithm>
#include <string>

#include "nimble_data_utils.h"
#include "nimble_defs.h"
#include "nimble_utils.h"
#include "nimble_view.h"

namespace nimble {

class DataManager;
class MaterialFactoryBase;
class MaterialFactory;

class MaterialParameters
{
 public:
  static const int MAX_NUM_MAT_PARAM     = 64;
  static const int MAX_MAT_MODEL_STR_LEN = 64;

  inline MaterialParameters() : num_material_points_(0) {}

  NIMBLE_INLINE_FUNCTION
  MaterialParameters(
      const std::string&                        material_name,
      const std::map<std::string, std::string>& string_params,
      const std::map<std::string, double>&      double_params,
      int                                       num_material_points = 0)
      : material_name_(material_name),
        material_string_parameters_(string_params),
        material_double_parameters_(double_params),
        num_material_points_(num_material_points)
  {
  }

  static inline void
  ConvertStringToUpperCase(std::string& s)
  {
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
  }

  inline ~MaterialParameters() = default;

  inline void
  AddParameter(const char* parameter_name, double parameter_value)
  {
    material_double_parameters_.insert(std::make_pair(std::string(parameter_name), parameter_value));
  }

  inline void
  AddStringParameter(const char* parameter_name, const char* parameter_value)
  {
    material_string_parameters_.insert(std::make_pair(std::string(parameter_name), std::string(parameter_value)));
  }

  inline bool
  IsParameter(const char* parameter_name) const
  {
    return material_double_parameters_.find(parameter_name) != material_double_parameters_.end();
  }

  inline bool
  IsStringParameter(const char* parameter_name) const
  {
    return material_string_parameters_.find(parameter_name) != material_string_parameters_.end();
  }

  inline std::string
  GetMaterialName(bool upper_case = false) const
  {
    std::string name = material_name_;
    if (upper_case) { ConvertStringToUpperCase(name); }
    return name;
  }

  inline int
  GetNumParameters() const
  {
    return static_cast<int>(material_double_parameters_.size());
  }

  inline int
  GetNumStringParameters() const
  {
    return static_cast<int>(material_string_parameters_.size());
  }

  inline double
  GetParameterValue(const char* parameter_name) const
  {
    try {
      return material_double_parameters_.at(parameter_name);
    } catch (...) {
      std::string errMsg = "Double parameter '" + std::string(parameter_name) + "' does not exist";
      throw std::invalid_argument(errMsg);
    }
  }

  inline const std::string&
  GetStringParameterValue(const char* parameter_name) const
  {
    try {
      return material_string_parameters_.at(parameter_name);
    } catch (...) {
      std::string errMsg = "String parameter '" + std::string(parameter_name) + "' does not exist";
      throw std::invalid_argument(errMsg);
    }
  }

  inline int
  GetNumMaterialPoints() const
  {
    return num_material_points_;
  }

  inline const std::map<std::string, double>&
  GetParameters() const
  {
    return material_double_parameters_;
  }

  inline const std::map<std::string, std::string>&
  GetStringParameters() const
  {
    return material_string_parameters_;
  }

  inline void
  Print() const
  {
    printf("\n--MaterialParameters\n");
    printf("  material name %s\n", material_name_.c_str());
    int num_material_double_parameters_ = static_cast<int>(material_double_parameters_.size());
    printf("  number of material double parameters %d\n", num_material_double_parameters_);
    for (const auto& p : material_double_parameters_) { printf("  %s %f\n", p.first.c_str(), p.second); }
    int num_material_string_parameters_ = static_cast<int>(material_string_parameters_.size());
    printf("  number of material string parameters %d\n", num_material_string_parameters_);
    for (const auto& p : material_string_parameters_) { printf("  %s %s\n", p.first.c_str(), p.second.c_str()); }
  }

 private:
  std::string material_name_;

  std::map<std::string, double>      material_double_parameters_;
  std::map<std::string, std::string> material_string_parameters_;

  int num_material_points_;
};

class Material
{
 public:
  NIMBLE_FUNCTION
  Material() = default;

  NIMBLE_FUNCTION
  Material(const Material& mat) = default;

  NIMBLE_FUNCTION
  virtual ~Material() = default;

  NIMBLE_FUNCTION
  virtual bool
  IsNGPLAMEModel() const
  {
    return false;
  }

  NIMBLE_FUNCTION
  virtual int
  NumStateVariables() const = 0;

  NIMBLE_FUNCTION
  virtual std::pair<std::string, nimble::Length>
  GetStateVariableLabelAndType(int index) const = 0;

  NIMBLE_FUNCTION
  virtual double
  GetStateVariableInitialValue(int index) const = 0;

  NIMBLE_FUNCTION
  virtual double
  GetDensity() const = 0;

  NIMBLE_FUNCTION
  virtual double
  GetBulkModulus() const = 0;

  NIMBLE_FUNCTION
  virtual double
  GetShearModulus() const = 0;

  NIMBLE_FUNCTION
  virtual void
  GetStress(
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
      bool          is_output_step) = 0;

#ifdef NIMBLE_HAVE_KOKKOS
  struct DeviceElemState
  {
    //-- Variable to store one scalar per integration point
    nimble_kokkos::DeviceScalarIntPtSingleEntryView scalar;
    //-- Variable to store one vector (dim 3) per integration point
    nimble_kokkos::DeviceVectorIntPtSingleEntryView vec3D;
    //-- Variable to store one symmetric tensor (6 entries) per integration point
    nimble_kokkos::DeviceSymTensorIntPtSingleEntryView sym_tensor;
    //-- Variable to store one full tensor (9 entries) per integration point
    nimble_kokkos::DeviceFullTensorIntPtSingleEntryView full_tensor;
  };

  NIMBLE_FUNCTION
  virtual void
  GetStress(
      double                                                     time_previous,
      double                                                     time_current,
      const nimble_kokkos::DeviceFullTensorIntPtSingleEntryView& deformation_gradient_n,
      const nimble_kokkos::DeviceFullTensorIntPtSingleEntryView& deformation_gradient_np1,
      const nimble_kokkos::DeviceSymTensorIntPtSingleEntryView&  stress_n,
      nimble_kokkos::DeviceSymTensorIntPtSingleEntryView         stress_np1,
      const nimble::Material::DeviceElemState&                   state_n,
      nimble::Material::DeviceElemState&                         state_np1) const
  {
    const int                        def_g_len = 9, stress_len = 6;
    nimble::Viewify<1, const double> def_g_n(deformation_gradient_n.data(), def_g_len);
    nimble::Viewify<1, const double> def_g_np1(deformation_gradient_np1.data(), def_g_len);
    nimble::Viewify<1, const double> s_n(stress_n.data(), stress_len);
    nimble::Viewify<1>               s_np1(stress_np1.data(), stress_len);
    GetStress(time_previous, time_current, def_g_n, def_g_np1, s_n, s_np1);
  }
#endif

  NIMBLE_FUNCTION
  virtual void
  GetTangent(int num_pts, double* material_tangent) const = 0;

#ifdef NIMBLE_HAVE_UQ
  NIMBLE_FUNCTION
  virtual void
  GetOffNominalStress(
      const double& bulk_mod,
      const double& shear_mod,
      int           num_pts,
      const double* deformation_gradient_np1,
      double*       stress_np1) = 0;
#endif

 protected:
  NIMBLE_FUNCTION
  virtual void
  GetStress(
      double                            time_previous,
      double                            time_current,
      nimble::Viewify<1, const double>& deformation_gradient_n,
      nimble::Viewify<1, const double>& deformation_gradient_np1,
      nimble::Viewify<1, const double>& stress_n,
      nimble::Viewify<1>                stress_np1) const = 0;
};

class ElasticMaterial : public Material
{
 public:
  static void
  register_supported_material_parameters(MaterialFactoryBase& factory);

  NIMBLE_FUNCTION
  explicit ElasticMaterial(MaterialParameters const& material_parameters);

  NIMBLE_FUNCTION
  ElasticMaterial(const ElasticMaterial& mat) = default;

  NIMBLE_FUNCTION
  int
  NumStateVariables() const override
  {
    return num_state_variables_;
  };

  NIMBLE_FUNCTION
  virtual std::pair<std::string, nimble::Length>
  GetStateVariableLabelAndType(int index) const override
  {
    printf(
        "\n**** Error, bad index in "
        "ElasticMaterial::GetStateVariableLabelAndType().\n");
    return std::make_pair("", nimble::LENGTH_0);
  }

  NIMBLE_FUNCTION
  double
  GetStateVariableInitialValue(int index) const override
  {
    printf(
        "\n**** Error, bad index in "
        "ElasticMaterial::GetStateVariableInitialValue().\n");
    return 0.0;
  }

  NIMBLE_FUNCTION
  double
  GetDensity() const override
  {
    return density_;
  }

  NIMBLE_FUNCTION
  double
  GetBulkModulus() const override
  {
    return bulk_modulus_;
  }

  NIMBLE_FUNCTION
  double
  GetShearModulus() const override
  {
    return shear_modulus_;
  }

  NIMBLE_FUNCTION
  void
  GetStress(
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
      bool          is_output_step) override;

  NIMBLE_FUNCTION
  void
  GetTangent(int num_pts, double* material_tangent) const override;

#ifdef NIMBLE_HAVE_UQ
  NIMBLE_FUNCTION
  void
  GetOffNominalStress(
      const double& bulk_mod,
      const double& shear_mod,
      int           num_pts,
      const double* deformation_gradient_np1,
      double*       stress_np1) override;
#endif

 protected:
  NIMBLE_FUNCTION
  void
  GetStress(
      double                            time_previous,
      double                            time_current,
      nimble::Viewify<1, const double>& deformation_gradient_n,
      nimble::Viewify<1, const double>& deformation_gradient_np1,
      nimble::Viewify<1, const double>& stress_n,
      nimble::Viewify<1>                stress_np1) const override;

 private:
  int    num_state_variables_;
  int    dim_;
  double density_;
  double bulk_modulus_;
  double shear_modulus_;
};

class NeohookeanMaterial : public Material
{
 public:
  static void
  register_supported_material_parameters(MaterialFactoryBase& factory);

  NIMBLE_FUNCTION
  NeohookeanMaterial(const NeohookeanMaterial& mat) = default;

  NIMBLE_FUNCTION
  explicit NeohookeanMaterial(MaterialParameters const& material_parameters);

  NIMBLE_FUNCTION
  int
  NumStateVariables() const override
  {
    return num_state_variables_;
  };

  NIMBLE_FUNCTION
  virtual std::pair<std::string, nimble::Length>
  GetStateVariableLabelAndType(int index) const override
  {
    printf(
        "\n**** Error, bad index in "
        "NeohookeanMaterial::GetStateVariableLabelAndType().\n");
    return std::make_pair("", nimble::LENGTH_0);
  }

  NIMBLE_FUNCTION
  double
  GetStateVariableInitialValue(int index) const override
  {
    printf(
        "\n**** Error, bad index in "
        "NeohookeanMaterial::GetStateVariableInitialValue().\n");
    return 0.0;
  }

  NIMBLE_FUNCTION
  double
  GetDensity() const override
  {
    return density_;
  }

  NIMBLE_FUNCTION
  double
  GetBulkModulus() const override
  {
    return bulk_modulus_;
  }

  NIMBLE_FUNCTION
  double
  GetShearModulus() const override
  {
    return shear_modulus_;
  }

  NIMBLE_FUNCTION
  void
  GetStress(
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
      bool          is_output_step) override;

  NIMBLE_FUNCTION
  void
  GetTangent(int num_pts, double* material_tangent) const override;

#ifdef NIMBLE_HAVE_UQ
  NIMBLE_FUNCTION
  void
  GetOffNominalStress(
      const double& bulk_mod,
      const double& shear_mod,
      int           num_pts,
      const double* deformation_gradient_np1,
      double*       stress_np1) override;
#endif

 protected:
  NIMBLE_FUNCTION
  void
  GetStress(
      double                            time_previous,
      double                            time_current,
      nimble::Viewify<1, const double>& deformation_gradient_n,
      nimble::Viewify<1, const double>& deformation_gradient_np1,
      nimble::Viewify<1, const double>& stress_n,
      nimble::Viewify<1>                stress_np1) const override;

 private:
  int    num_state_variables_{2};
  int    dim_{3};
  double density_{0.0};
  double bulk_modulus_{0.0};
  double shear_modulus_{0.0};
  double yield_stress_{0.0};
};

namespace j2 {

struct Properties
{
  double rho0{0.0};

  double E{0.0};
  double nu{0.0};

  double kappa{0.0};
  double mu{0.0};

  double Y0{0.0};
  double n{0.0};
  double eps0{0.0};

  double Svis0{0.0};
  double m{0.0};
  double eps_dot0{0.0};
};

inline double
HardeningPotential(Properties const props, double const eqps)
{
  double const& Y0   = props.Y0;
  double const& n    = props.n;
  double const& eps0 = props.eps0;

  if (n <= 0.0) return Y0 * eqps;

  double const exponent = (1.0 + n) / n;

  return Y0 * eps0 / exponent * (std::pow(1.0 + eqps / eps0, exponent) - 1.0);
}

inline double
FlowStrength(Properties const props, double const eqps)
{
  double const& Y0   = props.Y0;
  double const& n    = props.n;
  double const& eps0 = props.eps0;

  if (n <= 0.0) return Y0;

  return Y0 * std::pow(1.0 + eqps / eps0, 1.0 / n);
}

inline double
HardeningRate(Properties const props, double const eqps)
{
  double const& Y0   = props.Y0;
  double const& n    = props.n;
  double const& eps0 = props.eps0;

  if (n <= 0.0) return 0.0;

  return Y0 / (eps0 * n) * std::pow(1.0 + eqps / eps0, (1.0 - n) / n);
}

inline double
ViscoplasticDualKineticPotential(Properties const props, double const delta_eqps, double const dt)
{
  double const& Svis0    = props.Svis0;
  double const& m        = props.m;
  double const& eps_dot0 = props.eps_dot0;

  if (Svis0 <= 0.0) return 0.0;

  double const exponent = (1.0 + m) / m;
  double       psi_star = 0.0;
  if (delta_eqps > 0) {
    psi_star = dt * Svis0 * eps_dot0 / exponent * std::pow(delta_eqps / dt / eps_dot0, exponent);
  }
  return psi_star;
}

inline double
ViscoplasticStress(Properties const props, double const delta_eqps, double const dt)
{
  double const& Svis0    = props.Svis0;
  double const& m        = props.m;
  double const& eps_dot0 = props.eps_dot0;

  if (Svis0 <= 0.0 || m <= 0.0) return 0.0;

  double Svis = 0;
  if (delta_eqps > 0) {
    Svis = Svis0 * std::pow(delta_eqps / dt / eps_dot0, 1.0 / m);
  }

  return Svis;
}

inline double
ViscoplasticHardeningRate(Properties const props, double const delta_eqps, double const dt)
{
  double const& Svis0    = props.Svis0;
  double const& m        = props.m;
  double const& eps_dot0 = props.eps_dot0;

  if (Svis0 <= 0.0) return 0.0;

  double Hvis = 0;
  if (delta_eqps > 0) {
    Hvis = Svis0 / (eps_dot0 * m * dt) * std::pow(delta_eqps / dt / eps_dot0, (1.0 - m) / m);
  }

  return Hvis;
}

}  // namespace j2

class J2PlasticityMaterial : public Material
{
 public:
  static void
  register_supported_material_parameters(MaterialFactoryBase& factory);

  NIMBLE_FUNCTION
  J2PlasticityMaterial(const J2PlasticityMaterial& mat) = default;

  NIMBLE_FUNCTION
  explicit J2PlasticityMaterial(MaterialParameters const& material_parameters);

  NIMBLE_FUNCTION
  int
  NumStateVariables() const override
  {
    return num_state_variables_;
  };

  NIMBLE_FUNCTION
  virtual std::pair<std::string, nimble::Length>
  GetStateVariableLabelAndType(int index) const override
  {
    switch (index) {
    case 0:
      return std::make_pair("plastic_deformation_gradient_xx", nimble::SCALAR);
    case 1:
      return std::make_pair("plastic_deformation_gradient_yy", nimble::SCALAR);
    case 2:
      return std::make_pair("plastic_deformation_gradient_zz", nimble::SCALAR);
    case 3:
      return std::make_pair("plastic_deformation_gradient_xy", nimble::SCALAR);
    case 4:
      return std::make_pair("plastic_deformation_gradient_yz", nimble::SCALAR);
    case 5:
      return std::make_pair("plastic_deformation_gradient_zx", nimble::SCALAR);
    case 6:
      return std::make_pair("plastic_deformation_gradient_yx", nimble::SCALAR);
    case 7:
      return std::make_pair("plastic_deformation_gradient_zy", nimble::SCALAR);
    case 8:
      return std::make_pair("plastic_deformation_gradient_xz", nimble::SCALAR);
    case 9:
      return std::make_pair("equivalent_plastic_strain", nimble::SCALAR);
    default:
      std::cerr << std::endl << "**** Error, bad index in " << __PRETTY_FUNCTION__ << std::endl;
      exit(1);
    }
    return std::make_pair("", nimble::LENGTH_0);
  }
  
  NIMBLE_FUNCTION
  double
  GetStateVariableInitialValue(int index) const override
  {
    switch (index) {
    case 0:
      return 1.0;
    case 1:
      return 1.0;
    case 2:
      return 1.0;
    case 3:
      return 0.0;
    case 4:
      return 0.0;
    case 5:
      return 0.0;
    case 6:
      return 0.0;
    case 7:
      return 0.0;
    case 8:
      return 0.0;
    case 9:
      return 0.0;
    default:
      std::cerr << std::endl << "**** Error, bad index in " << __PRETTY_FUNCTION__ << std::endl;
      exit(1);
    }
    return 0.0;
  }

  NIMBLE_FUNCTION
  double
  GetDensity() const override
  {
    return props_.rho0;
  }

  NIMBLE_FUNCTION
  double
  GetBulkModulus() const override
  {
    return props_.kappa;
  }

  NIMBLE_FUNCTION
  double
  GetShearModulus() const override
  {
    return props_.mu;
  }

  NIMBLE_FUNCTION
  void
  GetStress(
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
      bool          is_output_step) override;

  NIMBLE_FUNCTION
  void
  GetTangent(int num_pts, double* material_tangent) const override;

#ifdef NIMBLE_HAVE_UQ
  NIMBLE_FUNCTION
  void
  GetOffNominalStress(
      const double& bulk_mod,
      const double& shear_mod,
      int           num_pts,
      const double* deformation_gradient_np1,
      double*       stress_np1) override;
#endif

 protected:
  NIMBLE_FUNCTION
  void
  GetStress(
      double                            time_previous,
      double                            time_current,
      nimble::Viewify<1, const double>& deformation_gradient_n,
      nimble::Viewify<1, const double>& deformation_gradient_np1,
      nimble::Viewify<1, const double>& stress_n,
      nimble::Viewify<1>                stress_np1) const override;

 private:
  int    num_state_variables_{10};
  int    dim_{3};
  j2::Properties props_;
};

}  // namespace nimble

#endif  // NIMBLE_MATERIAL_H
