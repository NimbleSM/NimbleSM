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

#include "nimble_data_manager.h"
#include <set>
#include <sstream>
#include <algorithm>
#include <ctype.h>
#include <stdexcept>

namespace nimble {

  ModelData& DataManager::GetMacroScaleData() {
    return macroscale_data_;
  }

  RVEData& DataManager::AllocateRVEData(int global_element_id,
                                        int integration_point_id) {
    std::pair<int, int> id_pair(global_element_id, integration_point_id);
    if (rve_data_.find(id_pair) != rve_data_.end()) {
      throw std::logic_error("\n****Error in DataManager::AllocateRVEData, requested global_element_id and integration_point_id have already been allocated.\n");
    }
    rve_data_[id_pair] = RVEData();
    RVEData& rve_data = rve_data_[id_pair];
    return rve_data;
  }

  RVEData& DataManager::GetRVEData(int global_element_id,
                                     int integration_point_id) {
    std::pair<int, int> id_pair(global_element_id, integration_point_id);
    RVEData& rve_data = rve_data_.at(id_pair);
    return rve_data;
  }

} // namespace nimble
