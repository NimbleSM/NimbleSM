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

#pragma once
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

namespace nimble
{
std::vector<id_pair> GetIDPairs(const std::vector<int>& global_node_ids)
{
  std::vector<id_pair> pairs{global_node_ids.size()};
  int len = global_node_ids.size();
  for (int i = 0; i < len; ++i)
  {
    pairs[i] = {global_node_ids[i], i};
  }
  return pairs;
}
std::vector<id_pair> GetSortedIDPairs(const std::vector<int>& global_node_ids)
{
  std::vector<id_pair> pairs{global_node_ids.size()};
  int len = global_node_ids.size();
  for (int i = 0; i < len; ++i)
  {
    pairs[i] = {global_node_ids[i], i};
  }
  id_pair::SortByGlobalID(pairs);
  return pairs;
}
std::vector<int> PackIDSpace(const std::vector<int>& raw_node_ids,
                             int max_nodes_assigned_to_a_rank,
                             const mpicontext& context)
{
  if (context.is_root())
  {
    std::map<int, int> lookup;
    std::vector<int> buffer(max_nodes_assigned_to_a_rank);
    int sane_node_id = 0;
    for (int node : raw_node_ids)
      if (lookup.count(node) == 0)
      {
        lookup[node] = sane_node_id;
        ++sane_node_id;
      }

    for (int rank = 1, context_size = context.get_size(); rank < context_size; ++rank)
    {
      MPI_Status recv_status;
      int recv_size;
      int err = MPI_Recv(
          buffer.data(), buffer.size(), MPI_INT, rank, 0xfae, context.get_comm(), &recv_status);
      if (err != MPI_SUCCESS)
        throw std::logic_error("Bad recv in PackIDSpace");
      MPI_Get_count(&recv_status, MPI_INT, &recv_size);
      for (int i = 0; i < recv_size; ++i)
      {
        int node = buffer[i];
        if (lookup.count(node) == 0)
        {
          lookup[node] = sane_node_id;
          ++sane_node_id;
        }
        buffer[i] = lookup[node];
      }
      MPI_Send(buffer.data(), recv_size, MPI_INT, rank, 0xfae, context.get_comm());
    }
    buffer.resize(raw_node_ids.size());
    buffer.shrink_to_fit();
    for (int i = 0, max = buffer.size(); i < max; ++i)
    {
      buffer[i] = lookup[raw_node_ids[i]];
    }
    // std::cout << "rank " << context.get_rank() << " new_ids: " << buffer << std::endl;
    return buffer;
  }
  else
  {
    std::vector<int> new_ids(raw_node_ids.size());
    MPI_Send(raw_node_ids.data(), raw_node_ids.size(), MPI_INT, 0, 0xfae, context.get_comm());
    MPI_Recv(
        new_ids.data(), new_ids.size(), MPI_INT, 0, 0xfae, context.get_comm(), MPI_STATUS_IGNORE);
    // std::cout << "rank " << context.get_rank() << " new_ids: " << new_ids << std::endl;
    return new_ids;
  }
}
int GetMaximumNodeId(const std::vector<int>& global_node_ids, const mpicontext& context)
{
  // std::cout << "GetMaximumNodeId recieved input " << global_node_ids << std::endl;
  int my_max_node_id = *std::max_element(global_node_ids.begin(), global_node_ids.end());
  int global_maximum = 0;
  MPI_Allreduce(&my_max_node_id, &global_maximum, 1, MPI_INT, MPI_MAX, context.get_comm());
  // std::cout << "Calculated global maximum to be " << global_maximum << std::endl;
  return global_maximum;
}
std::unique_ptr<int[]> GetNumberOfNodesAssignedToEachRank(const std::vector<int>& global_node_ids,
                                                          const mpicontext& context)
{
  int n_nodes_i_own{(int)global_node_ids.size()};
  std::unique_ptr<int[]> buffer{new int[context.get_size()]};
  MPI_Allgather(&n_nodes_i_own, 1, MPI_INT, buffer.get(), 1, MPI_INT, context.get_comm());
  return buffer;
}
void EnsureCheckpoint(const mpicontext& context, const std::string& message)
{
  MPI_Barrier(context.get_comm());
  if (context.is_root())
    std::cout << message << std::endl;
}
template<class Key, class Val>
std::vector<std::pair<Key, std::vector<Val>>> GroupConsecutive(
    const std::vector<std::pair<Key, Val>>& lst)
{
  if (lst.size() == 0)
    return {};
  std::vector<int> key_ids;
  key_ids.reserve(lst.size() + 1);
  const Key* group_key = &lst[0].first;
  int max_key_id       = 0;
  int lst_size         = lst.size();
  for (const std::pair<Key, Val>& p : lst)
  {
    if (p.first != *group_key)
    {
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
    for (int i = 1; i <= lst_size; ++i)
    {
      if (key_ids[i] != prev_key_id)
      {
        groups.emplace_back(lst[i - 1].first, std::vector<Val>(group_size));
        auto& new_elem                     = groups.back();
        const std::pair<Key, Val>* src_ptr = &lst[i - group_size];
        for (Val& q : new_elem.second)
        {
          q = src_ptr->second;
          ++src_ptr;
        }
        group_size  = 1;
        prev_key_id = key_ids[i];
      }
      else
      {
        ++group_size;
      }
    }
  }
  return groups;
}

}   // namespace nimble
