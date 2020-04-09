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

#ifndef NIMBLE_RVE_H
#define NIMBLE_RVE_H

#include "nimble_kokkos_defs.h"
#include "nimble_material.h"
#include "nimble_genesis_mesh.h"
#include "nimble_boundary_condition_manager.h"
#include "nimble_linear_solver.h"

#ifndef NIMBLE_HAVE_DARMA
  #include <vector>
  #include <map>
  #include <string>
#endif

class Block;

namespace nimble {

  class RVE : public Material {

  public:

    enum Boundary_Condition_Strategy {
      IMPOSE_DEFORMATION_GRADIENT = 0,
      PERIODIC_BC = 1,
      UNDEFINED = 2
    };

    NIMBLE_FUNCTION
    RVE(std::map<int, std::string> const & material_parameters_string,
        GenesisMesh const & rve_mesh,
        std::string rve_boundary_condition_strategy);

    NIMBLE_FUNCTION
    virtual ~RVE() {}

    std::map<std::string, double> ParseParametersString(std::string const & material_parameters_string) const;

    NIMBLE_FUNCTION
    int NumStateVariables() const { return 0; }

    NIMBLE_FUNCTION
    void GetStateVariableLabel(int index, char label[MaterialParameters::MAX_MAT_MODEL_STR_LEN]) const {
      printf("\n**** Error, bad index in RVE::GetStateVariableLabel().\n");
    }

    NIMBLE_FUNCTION
    double GetStateVariableInitialValue(int index) const  {
      printf("\n**** Error, bad index in RVE::GetStateVariableInitialValue().\n");
      return 0.0;
    }

    NIMBLE_FUNCTION
    double GetDensity() const ;

    NIMBLE_FUNCTION
    double GetBulkModulus() const ;

    NIMBLE_FUNCTION
    void InitializeRVE(int elem_global_id,
                       int integration_point_id,
                       DataManager& data_manager,
                       bool write_exodus_output,
                       MaterialFactory& factory) override;

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
                   nimble_kokkos::DeviceFullTensorSingleEntryView deformation_gradient_n,
                   nimble_kokkos::DeviceFullTensorSingleEntryView deformation_gradient_np1,
                   nimble_kokkos::DeviceSymTensorSingleEntryView stress_n,
                   nimble_kokkos::DeviceSymTensorSingleEntryView stress_np1) {}
#endif

    NIMBLE_FUNCTION
    void GetTangent(int num_pts,
                    double* material_tangent) const;

#ifdef NIMBLE_HAVE_UQ
    NIMBLE_FUNCTION
    void GetOffNominalStress(const std::vector<double> & params_this_sample,
                             const int & bulk_mod_idx,
                             const int & shear_mod_idx,
                             int num_pts,
                             const double * const deformation_gradient_np1,
                             double* stress_np1) {}
#endif

  private:

    std::map<int, std::string> material_parameters_string_;
    GenesisMesh rve_mesh_;
    Boundary_Condition_Strategy boundary_condition_strategy_;
    BoundaryConditionManager bc_;

    std::vector<int> linear_system_node_ids_;
    std::map<int, std::vector<int>> map_from_linear_system_;
    CGScratchSpace cg_scratch_;

    double origin_x_;
    double origin_y_;
    double origin_z_;

    int num_global_data_;
    std::vector<std::string> global_data_labels_;
    std::vector<double> global_data_;
  };

} // namespace nimble

#endif // NIMBLE_RVE_H
