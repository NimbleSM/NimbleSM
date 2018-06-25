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

#ifndef NIMBLE_MPI_REDUCTION_V4_H
#define NIMBLE_MPI_REDUCTION_V4_H

#ifdef NIMBLE_HAVE_MPI
#include <mpi.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>

#include "nimble.mpi.reduction.reduction_base.h"
#include "nimble_id_pair.h"

#ifndef NIMBLE_HAVE_DARMA
#include <vector>
#endif

namespace nimble
{
namespace reduction_v4
{
struct ReductionInfo : public ReductionInfoBase
{
  using indexmap = std::vector<std::pair<int, std::vector<int>>>;
  struct RankReductionInfo
  {
    std::unique_ptr<int[]> indicies;
    std::unique_ptr<double[]> sendbuffer;
    std::unique_ptr<double[]> recvbuffer;
    int buffersize;
    int n_indicies;
    int rank;
    RankReductionInfo() {}
    RankReductionInfo(int rank, const std::vector<int>& index_list, int buffersize) :
      indicies{new int[index_list.size()]},
      sendbuffer{new double[buffersize]},
      recvbuffer{new double[buffersize]},
      buffersize{buffersize},
      n_indicies{(int)index_list.size()},
      rank{rank}
    {
      std::copy_n(index_list.data(), index_list.size(), indicies.get());
    }
    RankReductionInfo(RankReductionInfo&& source) = default;
    RankReductionInfo& operator=(RankReductionInfo&& source) = default;
    bool okayfieldsizeQ(int field_size) { return field_size * n_indicies <= buffersize; }
    void fitnewfieldsize(int new_field_size) { resizebuffer(n_indicies * new_field_size); }
    void resizebuffer(int newbuffersize)
    {
      sendbuffer.reset(new double[newbuffersize]);
      recvbuffer.reset(new double[newbuffersize]);
      buffersize = newbuffersize;
    }
    template<int field_size>
    void fillsendbuffer(double* source)
    {
      double* scan   = sendbuffer.get();
      int *index_ptr = indicies.get(), *index_ptr_end = index_ptr + n_indicies;
      for (; index_ptr < index_ptr_end; ++index_ptr)
      {
        double* sourcescan = source + (*index_ptr) * field_size;
        for (int j = 0; j < field_size; ++j)
        {
          *scan = *sourcescan;
          ++scan;
          ++sourcescan;
        }
      }
    }
    void send_async(int field_size, MPI_Request* request, int tag = 0)
    {
      int count = field_size * n_indicies;
      MPI_Isend(sendbuffer.get(), count, MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, request);
    }
    void recv_async(int field_size, MPI_Request* request, int tag = 0)
    {
      int count = field_size * n_indicies;
      MPI_Irecv(recvbuffer.get(), count, MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, request);
    }
    template<int field_size>
    void reducefromrecvbuffer(double* source)
    {
      double* scan   = recvbuffer.get();
      int *index_ptr = indicies.get(), *index_ptr_end = index_ptr + n_indicies;
      for (; index_ptr < index_ptr_end; ++index_ptr)
      {
        double* sourcescan = source + (*index_ptr) * field_size;
        for (int j = 0; j < field_size; ++j)
        {
          *sourcescan += *scan;
          ++sourcescan;
          ++scan;
        }
      }
    }
  };
  std::vector<RankReductionInfo> info;
  std::unique_ptr<MPI_Request[]> recv_requests;
  std::unique_ptr<MPI_Request[]> send_requests;
  std::unique_ptr<int[]> unfinished;
  std::mt19937 gen;
  void SeedGenerator()
  {
    std::random_device rd;
    gen.seed(rd());
  }
  ReductionInfo() { SeedGenerator(); }
  ReductionInfo(const indexmap& m, int maximum_field_size) :
    recv_requests{new MPI_Request[m.size()]},
    send_requests{new MPI_Request[m.size()]},
    unfinished{new int[m.size()]}
  {
    info.reserve(m.size());
    for (auto& im : m)
    {
      info.emplace_back(im.first, im.second, (int)(maximum_field_size * im.second.size()));
    }
    SeedGenerator();
  }
  ReductionInfo(ReductionInfo&& ri) = default;
  ReductionInfo& operator=(ReductionInfo&& ri) = default;
  template<int field_size>
  void Reduce(double* data, int tag = 0)
  {
    int n_ranks = info.size();
    std::vector<int> i_ordering(n_ranks);
    for (int i = 0; i < n_ranks; ++i)
      i_ordering[i] = i;
    std::shuffle(i_ordering.begin(), i_ordering.end(), gen);
    for (int i = 0; i < n_ranks; ++i)
    {
      if (!info[i].okayfieldsizeQ(field_size))
        throw std::invalid_argument(
            "too big of a field size. Declare a larger maximum_field_size.");
    }
    for (int i = 0; i < n_ranks; ++i)
    {
      info[i].recv_async(field_size, &recv_requests[i], tag);
    }
    for (int i : i_ordering)
    {
      info[i].fillsendbuffer<field_size>(data);
      info[i].send_async(field_size, &send_requests[i], tag);
    }
    std::iota(unfinished.get(), unfinished.get() + n_ranks, 0);
    int n_unfinished = n_ranks;
    while (n_unfinished > 0)
    {
      int i = 0;
      while (i < n_unfinished)
      {
        int index = unfinished[i];
        int flag  = 0;
        MPI_Test(&recv_requests[index], &flag, MPI_STATUS_IGNORE);
        if (flag)
        {
          info[index].reducefromrecvbuffer<field_size>(data);
          n_unfinished -= 1;
          unfinished[i] = unfinished[n_unfinished];
        }
        else
          ++i;
      }
    }
  }
  static size_t total_indicies(const indexmap& i)
  {
    size_t size = 0;
    for (auto& rankindexpair : i)
      size += rankindexpair.second.size();
    return size;
  }
  void PerformReduction(double* data, int field_size)
  {
    switch (field_size)
    {
      case 1: Reduce<1>(data); break;
      case 2: Reduce<2>(data); break;
      case 3: Reduce<3>(data); break;
      default:
        std::string fs = std::to_string(field_size);
        std::cout << "Bad field size" << std::endl;
        throw std::invalid_argument("Bad field size of " + fs);
    }
  }
};
int GetMaximumNodeId(const std::vector<int>& global_node_ids)
{
  int my_max_node_id = *std::max(global_node_ids.begin(), global_node_ids.end());
  int global_maximum{};
  MPI_Allreduce(&my_max_node_id, &global_maximum, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return global_maximum;
}
std::unique_ptr<int[]> GetNumberOfNodesAssignedToEachRank(const std::vector<int>& global_node_ids,
                                                          int my_rank,
                                                          int total_ranks)
{
  int n_nodes_i_own{(int)global_node_ids.size()};
  std::unique_ptr<int[]> buffer{new int[total_ranks]};
  MPI_Allgather(&n_nodes_i_own, 1, MPI_INT, buffer.get(), 1, MPI_INT, MPI_COMM_WORLD);
  return buffer;
}
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

ReductionInfo::indexmap IdentifySharedNodes(const std::vector<int>& global_node_ids,
                                            int my_rank,
                                            int total_ranks)
{
  std::unique_ptr<int[]> node_counts{
      GetNumberOfNodesAssignedToEachRank(global_node_ids, my_rank, total_ranks)};
  int requisite_buffer_size = *std::max_element(node_counts.get(), node_counts.get() + total_ranks);
  std::vector<id_pair> my_ids = GetSortedIDPairs(global_node_ids);
  std::unique_ptr<id_pair[]> buffer{new id_pair[requisite_buffer_size]};
  std::unique_ptr<id_pair[]> intersection_buffer{new id_pair[requisite_buffer_size]};
  ReductionInfo::indexmap overall_mapping;
  for (int root = 0; root < total_ranks; ++root)
  {
    int n_nodes = node_counts[root];
    if (root == my_rank)
    {
      MPI_Bcast(my_ids.data(), 2 * n_nodes, MPI_INT, root, MPI_COMM_WORLD);
    }
    else
    {
      MPI_Bcast(buffer.get(), 2 * node_counts[root], MPI_INT, root, MPI_COMM_WORLD);
      id_pair* intersection_end = set_intersection(my_ids.begin(),
                                                   my_ids.end(),
                                                   buffer.get(),
                                                   buffer.get() + n_nodes,
                                                   intersection_buffer.get());
      int intersection_size     = std::distance(intersection_buffer.get(), intersection_end);
      if (intersection_size < 0)
        throw std::logic_error{
            "something went wrong when computing the intersection between different node ranges. "
            "Negative intersection size attained."};
      if (intersection_size > 0)
        // The vector passed to overall_mapping describes how this rank
        // Should package and unpackage data it's sending to the rank currently
        // Acting as the "root".
        overall_mapping.emplace_back(
            root, id_pair::GetMinorIDs(intersection_buffer.get(), intersection_size));
    }
  }
  return overall_mapping;
}
ReductionInfo* GenerateReductionInfo(const std::vector<int>& global_node_ids,
                                     int my_rank,
                                     int total_ranks)
{
  return new ReductionInfo{IdentifySharedNodes(global_node_ids, my_rank, total_ranks), 3};
}
void PerformReduction(double* data, int field_size, ReductionInfo& info)
{
  switch (field_size)
  {
    case 1: info.Reduce<1>(data); break;
    case 2: info.Reduce<2>(data); break;
    case 3: info.Reduce<3>(data); break;
    default:
      std::string fs = std::to_string(field_size);
      std::cout << "Bad field size" << std::endl;
      throw std::invalid_argument("Bad field size of " + fs);
  }
}
}   // namespace reduction_v4
}   // namespace nimble
#endif

#endif   // NIMBLE_MPI_REDUCTION_V4_H
