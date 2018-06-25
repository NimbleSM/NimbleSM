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
//#include "stdlib.h"
#include <string.h>

namespace nimble {

  class DataManager;

  class MaterialParameters {

  public:

    static const int MAX_NUM_MAT_PARAM = 64;
    static const int MAX_MAT_MODEL_STR_LEN = 64;

    NIMBLE_INLINE_FUNCTION
    MaterialParameters() {
      for (int i=0 ; i<MAX_NUM_MAT_PARAM ; ++i) {
        for (int j=0 ; j<MAX_MAT_MODEL_STR_LEN ; ++j) {
          material_parameter_names_[i][j] = '\0';
        }
        material_parameter_values_[i] = 0.0;
      }
    }

    NIMBLE_INLINE_FUNCTION
    MaterialParameters(int num_material_parameters,
                       const char material_parameter_names[MAX_NUM_MAT_PARAM][MAX_MAT_MODEL_STR_LEN],
                       const double material_parameter_values[MAX_NUM_MAT_PARAM]) {
      num_material_parameters_ = num_material_parameters;
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
    double GetParameterValue(const char* parameter_name) {
      for (int i=0 ; i<MAX_NUM_MAT_PARAM ; ++i) {
        if (StringsAreEqual(parameter_name, material_parameter_names_[i])) {
          return material_parameter_values_[i];
        }
      }
      return 0.0;
    }

  private:

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

    int num_material_parameters_ = 0;
    char material_parameter_names_[MAX_NUM_MAT_PARAM][MAX_MAT_MODEL_STR_LEN];
    double material_parameter_values_[MAX_NUM_MAT_PARAM];
  };

  void ParseMaterialParametersString(const char* material_parameters,
                                     int& num_material_parameters,
                                     char material_parameter_names[MaterialParameters::MAX_NUM_MAT_PARAM][MaterialParameters::MAX_MAT_MODEL_STR_LEN],
                                     double material_parameter_values[MaterialParameters::MAX_NUM_MAT_PARAM]);

  class Material {

    public:

    NIMBLE_FUNCTION
    Material(MaterialParameters const & material_parameters): material_parameters_(material_parameters) {}

    NIMBLE_FUNCTION
    virtual ~Material() {}

    NIMBLE_FUNCTION
    virtual int NumStateVariables() const = 0;

    NIMBLE_FUNCTION
    virtual void GetStateVariableLabel(int index, char* label) const = 0;

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
                               bool write_exodus_output) {}

    NIMBLE_FUNCTION
    virtual void GetStress(int elem_id,
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
                           bool is_output_step) = 0;

#ifdef NIMBLE_HAVE_KOKKOS
    NIMBLE_FUNCTION
    virtual void GetStress(double time_previous,
                           double time_current,
                           nimble_kokkos::DeviceFullTensorSingleEntryView deformation_gradient_n,
                           nimble_kokkos::DeviceFullTensorSingleEntryView deformation_gradient_np1,
                           nimble_kokkos::DeviceSymTensorSingleEntryView unrotated_stress_n,
                           nimble_kokkos::DeviceSymTensorSingleEntryView unrotated_stress_np1) = 0;
#endif

    NIMBLE_FUNCTION
    virtual void GetTangent(int num_pts,
                            double* material_tangent) const = 0;

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

  class NeohookeanMaterial : public Material {

  public:

    NIMBLE_FUNCTION
    NeohookeanMaterial(MaterialParameters const & material_parameters);

    NIMBLE_FUNCTION
    virtual ~NeohookeanMaterial() {}

    NIMBLE_FUNCTION
    int NumStateVariables() const { return num_state_variables_; };

    NIMBLE_FUNCTION
    void GetStateVariableLabel(int index, char* label) const {
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
                   const double * const unrotated_stress_n,
                   double* unrotated_stress_np1,
                   const double * const state_data_n,
                   double* state_data_np1,
                   DataManager& data_manager,
                   bool is_output_step);

#ifdef NIMBLE_HAVE_KOKKOS
    NIMBLE_FUNCTION
    void GetStress(double time_previous,
                   double time_current,
                   nimble_kokkos::DeviceFullTensorSingleEntryView deformation_gradient_n,
                   nimble_kokkos::DeviceFullTensorSingleEntryView deformation_gradient_np1,
                   nimble_kokkos::DeviceSymTensorSingleEntryView unrotated_stress_n,
                   nimble_kokkos::DeviceSymTensorSingleEntryView unrotated_stress_np1);
#endif

    NIMBLE_FUNCTION
    void GetTangent(int num_pts,
                    double* material_tangent) const ;

  private:

    int num_state_variables_;
    int dim_;
    double density_;
    double bulk_modulus_;
    double shear_modulus_;
  };

} // namespace nimble

#endif // NIMBLE_MATERIAL_H
