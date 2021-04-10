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

#include "nimble_block_base.h"
#include "nimble_element.h"

#ifdef NIMBLE_HAVE_DARMA
#include "darma.h"
#else
#include <memory>
#include <string>
#include <vector>
#endif

namespace nimble {
class MaterialFactory;
}

namespace nimble {

class Block : public nimble::BlockBase
{
 public:
  Block()
      : BlockBase(),
        vol_ave_volume_offset_(-1)
  {}

  ~Block() override = default;

#ifdef NIMBLE_HAVE_DARMA
  template <typename ArchiveType>
  void
  serialize(ArchiveType& ar)
  {
    // ar | macro_material_parameters_ | rve_material_parameters_ |
    // rve_boundary_condition_strategy_ | def_grad_offset_ | stress_offset_ |
    // state_data_offset_ | vol_ave_volume_offset_ | vol_ave_offsets_ |
    // vol_ave_index_to_derived_data_index_; ar | rve_output_global_elem_ids_;

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

  void
  Initialize(std::string const& macro_material_parameters, MaterialFactory& factory);

  void
  Initialize(
      std::string const&                macro_material_parameters,
      std::map<int, std::string> const& rve_material_parameters,
      GenesisMesh const&                rve_mesh,
      std::string                       rve_boundary_condition_strategy,
      MaterialFactory&                  factory);

  void
  InstantiateMaterialModel(MaterialFactory& factory);

  void
  InstantiateElement() override;

  void
  GetDataLabelsAndLengths(std::vector<std::pair<std::string, Length>>& data_labels_and_lengths);

  void
  ComputeLumpedMassMatrix(
      const double* const reference_coordinates,
      int                 num_elem,
      const int* const    elem_conn,
      double*             lumped_mass) const;

  void
  InitializeElementData(
      int                             num_elem_in_block,
      std::vector<int> const&         elem_global_ids_in_block,
      std::vector<int> const&         rve_output_global_elem_ids,
      std::vector<std::string> const& elem_data_labels,
      std::vector<std::string> const& derived_elem_data_labels,
      std::vector<double>&            elem_data_n,
      std::vector<double>&            elem_data_np1,
      MaterialFactory&                material_factory,
      DataManager&                    data_manager);

  void
  ComputeInternalForce(
      const double*                   reference_coordinates,
      const double*                   displacement,
      const double*                   velocity,
      const double*                   rve_macroscale_deformation_gradient,
      double*                         internal_force,
      double                          time_previous,
      double                          time_current,
      int                             num_elem,
      const int*                      elem_conn,
      const int*                      elem_global_ids,
      std::vector<std::string> const& elem_data_labels,
      std::vector<double> const&      elem_data_n,
      std::vector<double>&            elem_data_np1,
      DataManager&                    data_manager,
      bool                            is_output_step,
      bool compute_stress_only = false) const;

  void
  ComputeDerivedElementData(
      const double* const               reference_coordinates,
      const double* const               displacement,
      int                               num_elem,
      const int* const                  elem_conn,
      int                               num_elem_data,
      std::vector<double> const&        elem_data_np1,
      int                               num_derived_elem_data,
      std::vector<std::vector<double>>& derived_elem_data);

 protected:
  void
  DetermineDataOffsets(
      std::vector<std::string> const& elem_data_labels,
      std::vector<std::string> const& derived_elem_data_labels);

  std::vector<int>   def_grad_offset_;
  std::vector<int>   stress_offset_;
  std::vector<int>   state_data_offset_;
  int                vol_ave_volume_offset_;
  std::vector<int>   vol_ave_offsets_;
  std::map<int, int> vol_ave_index_to_derived_data_index_;

};

}  // namespace nimble

#endif  // NIMBLE_BLOCK_H
