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

#ifndef NIMBLE_BLOCK_H
#define NIMBLE_BLOCK_H

#include "nimble_element.h"
#include "nimble_material.h"
#include "nimble_genesis_mesh.h"

#ifdef NIMBLE_HAVE_DARMA
  #include "darma.h"
#else
  #include <string>
  #include <vector>
  #include <memory>
#endif

namespace nimble { class MaterialFactory; }

namespace nimble {

  class UqParameters;//Forward declaration to avoid circular header inclusion

  class Block {

  public:

  Block() : macro_material_parameters_("none"), vol_ave_volume_offset_(-1) , rve_boundary_condition_strategy_("none") {}

    virtual ~Block() { }

#ifdef NIMBLE_HAVE_DARMA
    template<typename ArchiveType>
    void serialize(ArchiveType& ar) {

      //ar | macro_material_parameters_ | rve_material_parameters_ | rve_boundary_condition_strategy_ | def_grad_offset_ | stress_offset_ | state_data_offset_ | vol_ave_volume_offset_ | vol_ave_offsets_ | vol_ave_index_to_derived_data_index_;
      //ar | rve_output_global_elem_ids_;


      ar | macro_material_parameters_;
      ar | rve_material_parameters_;
      ar | rve_boundary_condition_strategy_;
      ar | rve_output_global_elem_ids_;
      ar | rve_mesh_;
      // These are set up below
      // ar | element_;
      // ar | material_;
      ar | def_grad_offset_;
      ar | stress_offset_;
      ar | state_data_offset_;
      ar | vol_ave_volume_offset_;
      ar | vol_ave_offsets_;
      ar | vol_ave_index_to_derived_data_index_;
  
      // The purpose for the serialize call can be determined with:
      // ar.is_sizing()
      // ar.is_packing()
      // ar.is_unpacking()
  
      if (ar.is_unpacking()) {
        InstantiateMaterialModel();
        InstantiateElement();
      }
    }
#endif

    void Initialize(std::string const & macro_material_parameters,
                    MaterialFactory& factory);

    void Initialize(std::string const & macro_material_parameters,
                    std::map<int, std::string> const & rve_material_parameters,
                    GenesisMesh const & rve_mesh,
                    std::string rve_boundary_condition_strategy,
                    MaterialFactory& factory);

    void InstantiateMaterialModel(MaterialFactory& factory);

    void InstantiateElement();

    void GetDataLabelsAndLengths(std::vector< std::pair<std::string, Length> >& data_labels_and_lengths);

    double GetDensity() const {
      return material_->GetDensity();
    }

    double GetBulkModulus() const {
      return material_->GetBulkModulus();
    }

    void ComputeLumpedMassMatrix(const double * const reference_coordinates,
                                 int num_elem,
                                 const int * const elem_conn,
                                 double* lumped_mass) const ;

    double ComputeCriticalTimeStep(const double * const node_reference_coordinates,
                                   const double * const node_displacements,
                                   int num_elem,
                                   const int * const elem_conn) const;

    void InitializeElementData(int num_elem_in_block,
                               std::vector<int> const & elem_global_ids_in_block,
                               std::vector<int> const & rve_output_global_elem_ids,
                               std::vector<std::string> const & elem_data_labels,
                               std::vector<std::string> const & derived_elem_data_labels,
                               std::vector<double>& elem_data_n,
                               std::vector<double>& elem_data_np1,
                               MaterialFactory& material_factory,
                               DataManager& data_manager);

    void ComputeInternalForce(const double * const reference_coordinates,
                              const double * const displacement,
                              const double * const velocity,
                              const double * const rve_macroscale_deformation_gradient,
                              double * const internal_force,
                              double time_previous,
                              double time_current,
                              int num_elem,
                              const int * const elem_conn,
                              const int * const elem_global_ids,
                              std::vector<std::string> const & elem_data_labels,
                              std::vector<double> const & elem_data_n,
                              std::vector<double> & elem_data_np1,
                              DataManager& data_manager,
                              bool is_output_step,
#ifdef NIMBLE_HAVE_UQ
                              const bool & uq_enabled,
                              const UqParameters* uq_params,
                              const int & num_uq_samples,
                              const std::vector<double*> offnominal_displacements,
                              const std::vector<double*> offnominal_internal_forces,
                              const std::vector<double*> displacement_sensitivities,
#endif
                              bool compute_stress_only = false) const ;

    template <typename MatT>
    void ComputeTangentStiffnessMatrix(int num_global_unknowns,
                                       const double * const reference_coordinates,
                                       const double * const displacement,
                                       int num_elem,
                                       const int * const elem_conn,
                                       const int * const global_node_ids,
                                       MatT & tangent_stiffness) const ;

    void ComputeDerivedElementData(const double * const reference_coordinates,
                                   const double * const displacement,
                                   int num_elem,
                                   const int * const elem_conn,
                                   int num_elem_data,
                                   std::vector<double> const & elem_data_np1,
                                   int num_derived_elem_data,
                                   std::vector< std::vector<double> >& derived_elem_data);

    std::shared_ptr<Material> const GetMaterialPointer() const { return material_; }

#ifdef NIMBLE_HAVE_UQ
    void SetUqParamsIndexRange( const int & param_low,
                                const int & param_hi,
                                const std::vector<std::string> & param_names) {
      range_of_uq_params_ = std::make_pair(param_low, param_hi);
      bulk_modulus_uq_index_ = -1;
      shear_modulus_uq_index_ = -1;
      for(int p = param_low; p <= param_hi; p++) {
        if(param_names[p]=="bulk_modulus") { bulk_modulus_uq_index_ = p;}
        if(param_names[p]=="shear_modulus") { shear_modulus_uq_index_ = p;}
      }
    }
#endif

  private:

    void DetermineDataOffsets(std::vector<std::string> const & elem_data_labels,
                              std::vector<std::string> const & derived_elem_data_labels);

    std::string macro_material_parameters_;
    std::map<int, std::string> rve_material_parameters_;
    std::string rve_boundary_condition_strategy_;
    std::vector<int> rve_output_global_elem_ids_;
    // todo: can we avoid carrying the rve_mesh around?
    GenesisMesh rve_mesh_;
    std::shared_ptr<Element> element_ = nullptr;
    std::shared_ptr<Material> material_ = nullptr;
    std::vector<int> def_grad_offset_;
    std::vector<int> stress_offset_;
    std::vector<int> state_data_offset_;
    int vol_ave_volume_offset_;
    std::vector<int> vol_ave_offsets_;
    std::map<int, int> vol_ave_index_to_derived_data_index_;
#ifdef NIMBLE_HAVE_UQ
    std::pair<int, int> range_of_uq_params_;
    int bulk_modulus_uq_index_;
    int shear_modulus_uq_index_;
#endif
  };

} // namespace nimble

#endif // NIMBLE_BLOCK_H
