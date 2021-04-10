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

#ifndef NIMBLESM_NIMBLE_UQ_BLOCK_H
#define NIMBLESM_NIMBLE_UQ_BLOCK_H

#ifdef NIMBLE_HAVE_UQ

#include <memory>
#include <string>
#include <vector>

#include "nimble_block.h"

namespace nimble {
class DataManager;
}

namespace nimble_uq {

class Block : public nimble::Block
{
 public:
  Block() : nimble::Block(), bulk_modulus_uq_index_(-1), shear_modulus_uq_index_(-1) {}

  ~Block() override = default;

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
      nimble::DataManager&            data_manager,
      bool                            is_output_step,
      const bool&                     is_off_nominal,
      std::vector<double> const&      uq_params_this_sample,
      bool                            compute_stress_only = false) const;

  // HACK specific to elastic models
  void
  SetUqParameters(const std::map<std::string, int>& param_indices)
  {
    bulk_modulus_uq_index_  = -1;
    shear_modulus_uq_index_ = -1;
    for (auto const& it : param_indices) {
      if (it.first == "bulk_modulus") { bulk_modulus_uq_index_ = it.second; }
      if (it.first == "shear_modulus") { shear_modulus_uq_index_ = it.second; }
    }
  }

 protected:
  // HACK specific to elastic models
  //  std::pair<int, int> range_of_uq_params_;
  int bulk_modulus_uq_index_ = -1;
  int shear_modulus_uq_index_ = -1;
};

}  // namespace nimble_uq

#endif

#endif  // NIMBLESM_NIMBLE_UQ_BLOCK_H
