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

#ifndef NIMBLE_MPI_REDUCTION_UTILS
#define NIMBLE_MPI_REDUCTION_UTILS

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <exception>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "nimble.mpi.mpicontext.h"
#include "nimble.mpi.serialization.h"

namespace nimble {

std::vector<int>
PackIDSpace(
    const std::vector<int>& raw_node_ids,
    int                     max_nodes_assigned_to_a_rank,
    const mpicontext&       context);

int
GetMaximumNodeId(
    const std::vector<int>& global_node_ids,
    const mpicontext&       context);

std::unique_ptr<int[]>
GetNumberOfNodesAssignedToEachRank(
    const std::vector<int>& global_node_ids,
    const mpicontext&       context);

void
EnsureCheckpoint(const mpicontext& context, const std::string& message);

template <class Key, class Val>
std::vector<std::pair<Key, std::vector<Val>>>
GroupConsecutive(const std::vector<std::pair<Key, Val>>& lst)
{
  if (lst.size() == 0) return {};
  std::vector<int> key_ids;
  key_ids.reserve(lst.size() + 1);
  const Key* group_key  = &lst[0].first;
  int        max_key_id = 0;
  int        lst_size   = lst.size();
  for (const std::pair<Key, Val>& p : lst) {
    if (p.first != *group_key) {
      ++max_key_id;
      group_key = &p.first;
    }
    key_ids.push_back(max_key_id);
  }
  key_ids.push_back(max_key_id + 1);
  std::vector<std::pair<Key, std::vector<Val>>> groups;
  groups.reserve(max_key_id + 1);
  {
    int prev_key_id = 0;
    int group_size  = 1;
    for (int i = 1; i <= lst_size; ++i) {
      if (key_ids[i] != prev_key_id) {
        groups.emplace_back(lst[i - 1].first, std::vector<Val>(group_size));
        auto&                      new_elem = groups.back();
        const std::pair<Key, Val>* src_ptr  = &lst[i - group_size];
        for (Val& q : new_elem.second) {
          q = src_ptr->second;
          ++src_ptr;
        }
        group_size  = 1;
        prev_key_id = key_ids[i];
      } else {
        ++group_size;
      }
    }
  }
  return groups;
}

}  // namespace nimble

#endif  // NIMBLE_MPI_REDUCTION_UTILS
