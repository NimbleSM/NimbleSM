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
      throw std::runtime_error(errMsg);
    }
  }

  inline const std::string&
  GetStringParameterValue(const char* parameter_name) const
  {
    try {
      return material_string_parameters_.at(parameter_name);
    } catch (...) {
      std::string errMsg = "String parameter '" + std::string(parameter_name) + "' does not exist";
      throw std::runtime_error(errMsg);
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
  virtual void
  GetStateVariableLabel(int index, char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN]) const = 0;

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
  NIMBLE_FUNCTION
  virtual void
  GetStress(
      double                                                     time_previous,
      double                                                     time_current,
      const nimble_kokkos::DeviceFullTensorIntPtSingleEntryView& deformation_gradient_n,
      const nimble_kokkos::DeviceFullTensorIntPtSingleEntryView& deformation_gradient_np1,
      const nimble_kokkos::DeviceSymTensorIntPtSingleEntryView&  stress_n,
      nimble_kokkos::DeviceSymTensorIntPtSingleEntryView         stress_np1) const
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
  void
  GetStateVariableLabel(int index, char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN]) const override
  {
    printf(
        "\n**** Error, bad index in "
        "ElasticMaterial::GetStateVariableLabel().\n");
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
  void
  GetStateVariableLabel(int index, char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN]) const override
  {
    printf(
        "\n**** Error, bad index in "
        "NeohookeanMaterial::GetStateVariableLabel().\n");
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
  int    num_state_variables_;
  int    dim_;
  double density_;
  double bulk_modulus_;
  double shear_modulus_;
};

}  // namespace nimble

#endif  // NIMBLE_MATERIAL_H
