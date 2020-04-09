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

    NIMBLE_INLINE_FUNCTION
    MaterialParameters() : num_material_parameters_(0), num_material_points_(0) {
      for (int i=0 ; i<MAX_NUM_MAT_PARAM ; ++i) {
        for (int j=0 ; j<MAX_MAT_MODEL_STR_LEN ; ++j) {
          material_parameter_names_[i][j] = '\0';
        }
        material_parameter_values_[i] = 0.0;
      }
    }

    NIMBLE_INLINE_FUNCTION
    MaterialParameters(const char material_name[MAX_MAT_MODEL_STR_LEN],
                       int num_material_parameters,
                       const char material_parameter_names[MAX_NUM_MAT_PARAM][MAX_MAT_MODEL_STR_LEN],
                       const double material_parameter_values[MAX_NUM_MAT_PARAM],
                       int num_material_points = 0)
      : num_material_parameters_(num_material_parameters), num_material_points_(num_material_points) {
      for (int i=0 ; i<MAX_MAT_MODEL_STR_LEN ; ++i) {
        material_name_[i] = material_name[i];
      }
      for (int i=0 ; i<MAX_NUM_MAT_PARAM ; ++i) {
        for (int j=0 ; j<MAX_MAT_MODEL_STR_LEN ; ++j) {
          material_parameter_names_[i][j] = material_parameter_names[i][j];
        }
        material_parameter_values_[i] = material_parameter_values[i];
      }
    }

    NIMBLE_INLINE_FUNCTION
    ~MaterialParameters() {}

    NIMBLE_INLINE_FUNCTION
    void AddParameter(const char* parameter_name, double parameter_value) {
      int len = StringLength(parameter_name);
      for (int i=0 ; i<len ; ++i) {
        material_parameter_names_[num_material_parameters_][i] = parameter_name[i];
      }
      material_parameter_names_[num_material_parameters_][len] = '\0';
      material_parameter_values_[num_material_parameters_] = parameter_value;
      num_material_parameters_ += 1;
    }

    NIMBLE_INLINE_FUNCTION
    bool IsParameter(const char* parameter_name) const {
      for (int i=0 ; i<num_material_parameters_ ; ++i) {
        if (StringsAreEqual(parameter_name, material_parameter_names_[i])) {
          return true;
        }
      }
      return false;
    }

    NIMBLE_INLINE_FUNCTION
    void GetMaterialName(char material_name[MAX_MAT_MODEL_STR_LEN],
                         bool upper_case = false) const {
      for (int i=0 ; i<MAX_MAT_MODEL_STR_LEN ; ++i) {
        material_name[i] = material_name_[i];
      }
      if (upper_case) {
        ConvertStringToUpperCase(material_name);
      }
    }

    NIMBLE_INLINE_FUNCTION
    int GetNumParameters() const {
      return num_material_parameters_;
    }

    NIMBLE_INLINE_FUNCTION
    void GetParameterName(int index,
                          char parameter_name[MAX_MAT_MODEL_STR_LEN],
                          bool upper_case = false) const {
      for (int i=0 ; i<MAX_MAT_MODEL_STR_LEN ; ++i) {
        parameter_name[i] = material_parameter_names_[index][i];
      }
      if (upper_case) {
        ConvertStringToUpperCase(parameter_name);
      }
    }

    NIMBLE_INLINE_FUNCTION
    double GetParameterValue(const char* parameter_name) const {
      for (int i=0 ; i<num_material_parameters_ ; ++i) {
        if (StringsAreEqual(parameter_name, material_parameter_names_[i])) {
          return material_parameter_values_[i];
        }
      }
      return 0.0;
    }

    NIMBLE_INLINE_FUNCTION
    double GetParameterValue(int index) const {
      return material_parameter_values_[index];
    }

    NIMBLE_INLINE_FUNCTION
    int GetNumMaterialPoints() const {
      return num_material_points_;
    }

    NIMBLE_INLINE_FUNCTION
    void Print() const {
      printf("\n--MaterialParameters\n");
      printf("  material name %s\n", material_name_);
      printf("  number of material parameters %d\n", num_material_parameters_);
      for (int i=0 ; i<num_material_parameters_ ; ++i) {
        printf("  %s %f\n", material_parameter_names_[i], material_parameter_values_[i]);
      }
    }

  private:

    NIMBLE_INLINE_FUNCTION
    void ConvertStringToUpperCase(char s[MAX_MAT_MODEL_STR_LEN]) const {
      for (int i=0 ; i<MAX_MAT_MODEL_STR_LEN ; ++i) {
        if (s[i] == 'a') s[i] = 'A';
        if (s[i] == 'b') s[i] = 'B';
        if (s[i] == 'c') s[i] = 'C';
        if (s[i] == 'd') s[i] = 'D';
        if (s[i] == 'e') s[i] = 'E';
        if (s[i] == 'f') s[i] = 'F';
        if (s[i] == 'g') s[i] = 'G';
        if (s[i] == 'h') s[i] = 'H';
        if (s[i] == 'i') s[i] = 'I';
        if (s[i] == 'j') s[i] = 'J';
        if (s[i] == 'k') s[i] = 'K';
        if (s[i] == 'l') s[i] = 'L';
        if (s[i] == 'm') s[i] = 'M';
        if (s[i] == 'n') s[i] = 'N';
        if (s[i] == 'o') s[i] = 'O';
        if (s[i] == 'p') s[i] = 'P';
        if (s[i] == 'q') s[i] = 'Q';
        if (s[i] == 'r') s[i] = 'R';
        if (s[i] == 's') s[i] = 'S';
        if (s[i] == 't') s[i] = 'T';
        if (s[i] == 'u') s[i] = 'U';
        if (s[i] == 'v') s[i] = 'V';
        if (s[i] == 'w') s[i] = 'W';
        if (s[i] == 'x') s[i] = 'X';
        if (s[i] == 'y') s[i] = 'Y';
        if (s[i] == 'z') s[i] = 'Z';
      }
    }

    char material_name_[MAX_MAT_MODEL_STR_LEN];
    int num_material_parameters_ = 0;
    char material_parameter_names_[MAX_NUM_MAT_PARAM][MAX_MAT_MODEL_STR_LEN];
    double material_parameter_values_[MAX_NUM_MAT_PARAM];
    int num_material_points_;
  };

  class Material {

    public:

    NIMBLE_FUNCTION
    Material(MaterialParameters const & material_parameters): material_parameters_(material_parameters) {}

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

    NIMBLE_FUNCTION
    MaterialParameters const & GetMaterialParameters() const { return material_parameters_;}

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

    MaterialParameters material_parameters_;
  };

  class ElasticMaterial : public Material {

  public:

    NIMBLE_FUNCTION
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

    NIMBLE_FUNCTION
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
