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

#ifndef NIMBLE_DATA_MANAGER_H
#define NIMBLE_DATA_MANAGER_H

#include "nimble_block.h"
#include "nimble_data_utils.h"
#include "nimble_exodus_output.h"
#include "nimble_genesis_mesh.h"
#include "nimble_linear_solver.h"
#include "nimble_model_data_base.h"
#include "nimble_parser.h"

#ifdef NIMBLE_HAVE_DARMA
  #include "darma.h"
#else
  #include <vector>
  #include <string>
  #include <map>
#endif

namespace nimble {

class ModelDataBase;

struct RVEData {
  std::shared_ptr<nimble::ModelDataBase> model_data_;
  std::vector<double> residual_vector_;
  std::vector<double> linear_solver_solution_;
  CRSMatrixContainer tangent_stiffness_;
  bool write_exodus_output_;
  ExodusOutput exodus_output_;
  std::map<int, std::vector< std::vector<double> > > derived_elem_data_;

#ifdef NIMBLE_HAVE_DARMA
    template<typename ArchiveType>
    void serialize(ArchiveType& ar) {
      ar | model_data_ | residual_vector_ | linear_solver_solution_ | tangent_stiffness_ | write_exodus_output_ | exodus_output_ | derived_elem_data_;
    }
#endif
};

class DataManager {

public:

  explicit DataManager(const nimble::Parser &parser);

  virtual ~DataManager() = default;

#ifdef NIMBLE_HAVE_DARMA
  template<typename ArchiveType>
  void serialize(ArchiveType& ar) {
    ar | macroscale_data_ | rve_data_;
  }
#endif

  std::shared_ptr<nimble::ModelDataBase> GetMacroScaleData()
  { return macroscale_data_; }

  RVEData& AllocateRVEData(int global_element_id,
                           int integration_point_id);

  RVEData& GetRVEData(int global_element_id,
                      int integration_point_id);

 protected:

  const nimble::Parser &parser_;
  std::shared_ptr<nimble::ModelDataBase> macroscale_data_;
  std::map<std::pair<int,int>, RVEData> rve_data_;

 };

} // namespace nimble

#endif
