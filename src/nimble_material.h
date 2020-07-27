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

#include "nimble_kokkos_defs.h"
#include "nimble_data_utils.h"
#include <string.h>

namespace nimble {

  class DataManager;
  class MaterialFactoryBase;
  class MaterialFactory;

  NIMBLE_INLINE_FUNCTION
  int StringLength(const char* str) {
    int len = 0;
    while (str[len] != '\0') {
      len++;
    }
    return len;
  }

  NIMBLE_INLINE_FUNCTION
  bool StringsAreEqual(const char* str1,
                       const char* str2) {
    int len1 = StringLength(str1);
    int len2 = StringLength(str2);
    int len = len1 < len2 ? len1 : len2;
    bool are_equal = true;
    for (int i=0 ; i<len ; ++i) {
      if (str1[i] != str2[i]) {
        are_equal = false;
        break;
      }
    }
    return are_equal;
  }

  class MaterialParameters {

  public:

    static const int MAX_NUM_MAT_PARAM = 64;
    static const int MAX_MAT_MODEL_STR_LEN = 64;

    inline
    MaterialParameters() : num_material_points_(0) {
    }

  inline
  MaterialParameters(const std::string& material_name,
                     const std::map<std::string, std::string>& string_params,
                     const std::map<std::string, double>& double_params,
                     int num_material_points = 0)
      :
      material_name_(material_name),
      material_string_parameters_(string_params),
      material_double_parameters_(double_params),
      num_material_points_(num_material_points) {
  }

    inline
    ~MaterialParameters() {}

    inline
    void AddParameter(const char* parameter_name, double parameter_value) {
      material_double_parameters_.insert(std::make_pair(std::string(parameter_name), parameter_value));
    }

    inline
    void AddStringParameter(const char* parameter_name, const char* parameter_value) {
      material_string_parameters_.insert(std::make_pair(std::string(parameter_name), std::string(parameter_value)));
    }

    inline
    bool IsParameter(const char* parameter_name) const {
      return material_double_parameters_.find(parameter_name) != material_double_parameters_.end();
    }

    inline
    bool IsStringParameter(const char* parameter_name) const {
      return material_string_parameters_.find(parameter_name) != material_string_parameters_.end();
    }

  inline std::string GetMaterialName(bool upper_case = false) const {
    std::string name = material_name_;
    if (upper_case) {
      ConvertStringToUpperCase(name);
    }
    return name;
  }

    inline
    int GetNumParameters() const {
      return material_double_parameters_.size();
    }

    inline
    int GetNumStringParameters() const {
      return material_string_parameters_.size();
    }

    inline
    double GetParameterValue(const char* parameter_name) const {
      return material_double_parameters_.at(parameter_name);
    }

    inline
    const std::string& GetStringParameterValue(const char* parameter_name) const {
      return material_string_parameters_.at(parameter_name);
    }

    inline
    int GetNumMaterialPoints() const {
      return num_material_points_;
    }

    inline
    const std::map<std::string, double>&
    GetParameters() const {
      return material_double_parameters_;
    }

    inline
    const std::map<std::string, std::string>&
    GetStringParameters() const {
      return material_string_parameters_;
    }

    inline
    void Print() const {
      printf("\n--MaterialParameters\n");
      printf("  material name %s\n", material_name_.c_str());
      int num_material_double_parameters_ = material_double_parameters_.size();
      printf("  number of material double parameters %d\n", num_material_double_parameters_);
      for (auto p : material_double_parameters_) {
        printf("  %s %f\n", p.first.c_str(), p.second);
      }
      int num_material_string_parameters_ = material_string_parameters_.size();
      printf("  number of material string parameters %d\n", num_material_string_parameters_);
      for (auto p : material_string_parameters_) {
        printf("  %s %s\n", p.first.c_str(), p.second.c_str());
      }
    }

  private:

    inline void ConvertStringToUpperCase(std::string& s) const {
      for(auto&& c : s) {
        c = std::toupper(c);
      }
    }

    std::string material_name_;

    std::map<std::string, double> material_double_parameters_;
    std::map<std::string, std::string> material_string_parameters_;

    int num_material_points_;
  };

  class Material {

    public:

    NIMBLE_FUNCTION
    Material() {}

    NIMBLE_FUNCTION
    virtual ~Material() {}

    NIMBLE_FUNCTION
    virtual bool IsNGPLAMEModel() const { return false; }

    NIMBLE_FUNCTION
    virtual int NumStateVariables() const = 0;

    NIMBLE_FUNCTION
    virtual void GetStateVariableLabel(int index, char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN]) const = 0;

    NIMBLE_FUNCTION
    virtual double GetStateVariableInitialValue(int index) const = 0;

    NIMBLE_FUNCTION
    virtual double GetDensity() const = 0;

    NIMBLE_FUNCTION
    virtual double GetBulkModulus() const = 0;

    NIMBLE_FUNCTION
    virtual void InitializeRVE(int elem_global_id,
                               int integration_point_id,
                               DataManager& data_manager,
                               bool write_exodus_output,
                               MaterialFactory& factory) {}

    NIMBLE_FUNCTION
    virtual void GetStress(int elem_id,
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
                           bool is_output_step) = 0;

#ifdef NIMBLE_HAVE_KOKKOS
    NIMBLE_FUNCTION
    virtual void GetStress(double time_previous,
                           double time_current,
                           nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_n,
                           nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_np1,
                           nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_n,
                           nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_np1) = 0;
#endif

    NIMBLE_FUNCTION
    virtual void GetTangent(int num_pts,
                            double* material_tangent) const = 0;

#ifdef NIMBLE_HAVE_UQ
    NIMBLE_FUNCTION
    virtual void GetOffNominalStress(const std::vector<double> & params_this_sample,
                                     const int & bulk_mod_idx,
                                     const int & shear_mod_idx,
                                     int num_pts,
                                     const double * const deformation_gradient_np1,
                                     double* stress_np1) = 0;
#endif

  protected:

    const int K_S_XX_ = 0 ;
    const int K_S_YY_ = 1 ;
    const int K_S_ZZ_ = 2 ;
    const int K_S_XY_ = 3 ;
    const int K_S_YZ_ = 4 ;
    const int K_S_ZX_ = 5 ;
    const int K_S_YX_ = 3 ;
    const int K_S_ZY_ = 4 ;
    const int K_S_XZ_ = 5 ;

    const int K_F_XX_ = 0 ;
    const int K_F_YY_ = 1 ;
    const int K_F_ZZ_ = 2 ;
    const int K_F_XY_ = 3 ;
    const int K_F_YZ_ = 4 ;
    const int K_F_ZX_ = 5 ;
    const int K_F_YX_ = 6 ;
    const int K_F_ZY_ = 7 ;
    const int K_F_XZ_ = 8 ;
  };

  class ElasticMaterial : public Material {

  public:
    static void register_supported_material_parameters(MaterialFactoryBase& factory);

    ElasticMaterial(MaterialParameters const & material_parameters);

    NIMBLE_FUNCTION
    virtual ~ElasticMaterial() {}

    NIMBLE_FUNCTION
    int NumStateVariables() const { return num_state_variables_; };

    NIMBLE_FUNCTION
    void GetStateVariableLabel(int index, char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN]) const {
      printf("\n**** Error, bad index in ElasticMaterial::GetStateVariableLabel().\n");
    }

    NIMBLE_FUNCTION
    double GetStateVariableInitialValue(int index) const  {
      printf("\n**** Error, bad index in ElasticMaterial::GetStateVariableInitialValue().\n");
      return 0.0;
    }

    NIMBLE_FUNCTION
    double GetDensity() const { return density_; }

    NIMBLE_FUNCTION
    double GetBulkModulus() const { return bulk_modulus_; }

    NIMBLE_FUNCTION
    void GetStress(int elem_id,
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
                   bool is_output_step);

#ifdef NIMBLE_HAVE_KOKKOS
    NIMBLE_FUNCTION
    void GetStress(double time_previous,
                   double time_current,
                   nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_n,
                   nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_np1,
                   nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_n,
                   nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_np1);
#endif

    NIMBLE_FUNCTION
    void GetTangent(int num_pts,
                    double* material_tangent) const ;

#ifdef NIMBLE_HAVE_UQ
    NIMBLE_FUNCTION
    void GetOffNominalStress(const std::vector<double> & params_this_sample,
                             const int & bulk_mod_idx,
                             const int & shear_mod_idx,
                             int num_pts,
                             const double * const deformation_gradient_np1,
                             double* stress_np1);
#endif

  private:

    int num_state_variables_;
    int dim_;
    double density_;
    double bulk_modulus_;
    double shear_modulus_;
  };

  class NeohookeanMaterial : public Material {

  public:
    static void register_supported_material_parameters(MaterialFactoryBase& factory);

    NeohookeanMaterial(MaterialParameters const & material_parameters);

    NIMBLE_FUNCTION
    virtual ~NeohookeanMaterial() {}

    NIMBLE_FUNCTION
    int NumStateVariables() const { return num_state_variables_; };

    NIMBLE_FUNCTION
    void GetStateVariableLabel(int index, char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN]) const {
      printf("\n**** Error, bad index in NeohookeanMaterial::GetStateVariableLabel().\n");
    }

    NIMBLE_FUNCTION
    double GetStateVariableInitialValue(int index) const  {
      printf("\n**** Error, bad index in NeohookeanMaterial::GetStateVariableInitialValue().\n");
      return 0.0;
    }

    NIMBLE_FUNCTION
    double GetDensity() const { return density_; }

    NIMBLE_FUNCTION
    double GetBulkModulus() const { return bulk_modulus_; }

    NIMBLE_FUNCTION
    void GetStress(int elem_id,
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
                   bool is_output_step);

#ifdef NIMBLE_HAVE_KOKKOS
    NIMBLE_FUNCTION
    void GetStress(double time_previous,
                   double time_current,
                   nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_n,
                   nimble_kokkos::DeviceFullTensorIntPtSingleEntryView deformation_gradient_np1,
                   nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_n,
                   nimble_kokkos::DeviceSymTensorIntPtSingleEntryView stress_np1);
#endif

    NIMBLE_FUNCTION
    void GetTangent(int num_pts,
                    double* material_tangent) const ;

#ifdef NIMBLE_HAVE_UQ
    NIMBLE_FUNCTION
    void GetOffNominalStress(const std::vector<double> & params_this_sample,
                             const int & bulk_mod_idx,
                             const int & shear_mod_idx,
                             int num_pts,
                             const double * const deformation_gradient_np1,
                             double* stress_np1);
#endif

  private:

    int num_state_variables_;
    int dim_;
    double density_;
    double bulk_modulus_;
    double shear_modulus_;
  };

} // namespace nimble

#endif // NIMBLE_MATERIAL_H
