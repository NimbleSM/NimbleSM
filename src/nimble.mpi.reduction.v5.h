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

#ifndef NIMBLE_MPI_REDUCTION_V5_H
#define NIMBLE_MPI_REDUCTION_V5_H

#ifdef NIMBLE_HAVE_MPI
#include <mpi.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>

#include "nimble.mpi.rank_clique_reducer.h"
#include "nimble.mpi.reduction.reduction_base.h"
#include "nimble.mpi.reduction_utils.h"
#include "nimble_id_pair.h"

#ifndef NIMBLE_HAVE_DARMA
#include <vector>
#endif

namespace nimble
{
namespace reduction_v5
{
// globalid is a global node id
// localid is a local node id
typedef int globalid;
typedef int localid;
typedef int rankid;
/*To do: find and replace color_t with colorid */
// color_t represents the color identifier used to create a MPI_Comm type
// When MPI_Comm_split is called, all processors calling MPI_Comm_split with
// The same value passed for color will be grouped into the same comm.
typedef int color_t;   // HAS NOTHING TO DO WITH GRAPHICS OR IMAGES
using length_v   = std::vector<int>;
using rank_v     = std::vector<int>;
using globalid_v = std::vector<globalid>;
using localid_v  = std::vector<localid>;
using color_v    = std::vector<color_t>;
using std::pair;
using std::vector;
struct ReductionInfo : public ReductionInfoBase
{
  struct ReductionClique_tOld
  {
    std::unique_ptr<int[]> indicies;
    std::unique_ptr<double[]> sendbuffer;
    std::unique_ptr<double[]> recvbuffer;
    int buffersize;
    int n_indicies;
    MPI_Comm clique_comm;
    MPI_Request Iallreduce_request;
    bool exists_active_asyncreduce_request = false;
    ReductionClique_tOld() {}
    /* README
    All processors calling MPI_Comm_split with the same comm_color will
    be grouped into the same communicator group. This means that this ReductionClique_t
    will share a clique_comm with all the other ReductionClique_ts created with the
    given comm_color. The way this works is that if a set of ranks R all share
    a set of nodes N (so that every rank in R has every node in N), all ranks
    in R will be grouped into a single MPI_Comm (identified by the comm_color).
    The bookkeeping for this new MPI_Comm is stored in a ReductionClique_t struct,
    along with the bookkeeping necessary to call Allreduce to reduce the data
    for all the nodes in N.
    See http://mpitutorial.com/tutorials/introduction-to-groups-and-communicators/
    for a reference
    */
    ReductionClique_tOld(const std::vector<int>& index_list,
                         int buffersize,
                         int comm_color,
                         const mpicontext& context) :
      indicies{new int[index_list.size()]},
      sendbuffer{new double[buffersize]},
      recvbuffer{new double[buffersize]},
      buffersize{buffersize},
      n_indicies{(int)index_list.size()},
      // Actual instantiation of clique_comm is performed by MPI_Comm_split
      clique_comm{}
    {
      std::copy_n(index_list.data(), index_list.size(), indicies.get());
      // printf("rank %d: Splitting %d with comm color %d\n", context.get_rank(),
      // context.get_comm(), comm_color);
      MPI_Comm_split(context.get_comm(), comm_color, context.get_rank(), &clique_comm);
    }
    ReductionClique_tOld(ReductionClique_tOld&& source) = default;
    ReductionClique_tOld& operator=(ReductionClique_tOld&& source) = default;
    bool okayfieldsizeQ(int field_size) { return field_size * n_indicies <= buffersize; }
    void fitnewfieldsize(int new_field_size) { resizebuffer(n_indicies * new_field_size); }
    void resizebuffer(int newbuffersize)
    {
      sendbuffer.reset(new double[newbuffersize]);
      recvbuffer.reset(new double[newbuffersize]);
      buffersize = newbuffersize;
    }
    template<int field_size>
    void pack(double* source)
    {
      double* destscan = sendbuffer.get();
      int *index_ptr = indicies.get(), *index_ptr_end = index_ptr + n_indicies;
      for (; index_ptr < index_ptr_end; ++index_ptr)
      {
        double* sourcescan = source + (*index_ptr) * field_size;
        for (int j = 0; j < field_size; ++j)
        {
          *destscan = *sourcescan;
          ++destscan;
          ++sourcescan;
        }
      }
    }
    template<int field_size>
    void unpack(double* dest)
    {
      double* sourcescan = recvbuffer.get();
      int *index_ptr = indicies.get(), *index_ptr_end = index_ptr + n_indicies;
      for (; index_ptr < index_ptr_end; ++index_ptr)
      {
        double* destscan = dest + (*index_ptr) * field_size;
        for (int j = 0; j < field_size; ++j)
        {
          *destscan = *sourcescan;
          ++destscan;
          ++sourcescan;
        }
      }
    }
    void EnsureDataSafety(int field_size)
    {
      if (!okayfieldsizeQ(field_size))
      {
        throw std::logic_error("Field size too big");
      }
    }
    // Performs a blocking Allreduce
    template<int field_size>
    void reduce(double* data)
    {
      EnsureDataSafety(field_size);
      pack<field_size>(data);
      MPI_Allreduce(sendbuffer.get(),
                    recvbuffer.get(),
                    n_indicies * field_size,
                    MPI_DOUBLE,
                    MPI_SUM,
                    clique_comm);
      unpack<field_size>(data);
    }
    // Copies data to sendbuffer and starts an asynchronous allreduce operation.
    // asyncreduce_finalize must be called for the data to be copied from the
    // recvbuffer to the databuffer
    template<int field_size>
    void asyncreduce_initialize(double* databuffer)
    {
      if (exists_active_asyncreduce_request)
        throw std::logic_error(
            "asyncreduce_initialize(data) was called "
            "when an active asynchronous reduce already exists");

      EnsureDataSafety(field_size);
      pack<field_size>(databuffer);
      MPI_Iallreduce(sendbuffer.get(),
                     recvbuffer.get(),
                     n_indicies * field_size,
                     MPI_DOUBLE,
                     MPI_SUM,
                     clique_comm,
                     &Iallreduce_request);
      exists_active_asyncreduce_request = true;
    }
    // Returns true if the currently active asynchronous reduce request has completed
    // Returns false if the currently active asynchronous reduce request hasn't completed
    // Throws an exception if there's no currently active asynchronous reduce request
    bool Iallreduce_completedQ()
    {
      if (exists_active_asyncreduce_request)
      {
        int flag = 0;
        MPI_Test(&Iallreduce_request, &flag, MPI_STATUS_IGNORE);
        // flag is set to true if the allreduce has completed
        return flag != 0;
      }
      else
        throw std::logic_error(
            "Iallreduce_completedQ() was called "
            "without an active asynchronous reduce request.");
    }
    // Returns true and unpacks data if the currently active asynchronous reduce request has
    // completed Returns false if the currently active asynchronous reduce request hasn't completed
    // Throws an exception if there's no currently active asynchronous reduce request
    template<int field_size>
    bool asyncreduce_finalize(double* databuffer)
    {
      if (exists_active_asyncreduce_request)
      {
        if (Iallreduce_completedQ())
        {
          unpack<field_size>(databuffer);
          exists_active_asyncreduce_request = false;
          return true;
        }
        else
        {
          return false;
        }
      }
      else
        throw std::logic_error(
            "asyncreduce_finalize(data) was called "
            "without an active asynchronous reduce request.");
    }
  };
  std::vector<ReductionClique_t> cliques;
  std::unique_ptr<int[]> unfinished;

  ReductionInfo() {}
  ReductionInfo(
      const std::pair<int, std::vector<std::pair<int, std::vector<int>>>>& MaxColorAndRawCliqueData,
      const mpicontext& context,
      int maximum_field_size = 3) :
    cliques{},
    unfinished{new int[MaxColorAndRawCliqueData.second.size()]}
  {
    int max_color       = MaxColorAndRawCliqueData.first;
    auto& RawCliqueData = MaxColorAndRawCliqueData.second;
    MPI_Comm discard_comm;
    cliques.reserve(RawCliqueData.size());
    {
      auto RawCliqueDataIterator = RawCliqueData.begin();
      for (int color = 0; color < max_color; ++color)
      {
        if (RawCliqueDataIterator == RawCliqueData.end() || (*RawCliqueDataIterator).first > color)
          MPI_Comm_split(context.get_comm(), MPI_UNDEFINED, 0, &discard_comm);
        else if ((*RawCliqueDataIterator).first == color)
        {
          auto& color_nodelst_pair = *RawCliqueDataIterator;
          ++RawCliqueDataIterator;
          auto& nodelist = color_nodelst_pair.second;
          int color      = color_nodelst_pair.first;

          int buffersize = maximum_field_size * nodelist.size();
          cliques.emplace_back(nodelist, buffersize, color, context);
        }
        else
          throw std::invalid_argument("Something went horribly wrong in ReductionInfo");
        // MPI_Barrier(context.get_comm());
      }
    }
  }
  ReductionInfo(ReductionInfo&& ri) = default;
  ReductionInfo& operator=(ReductionInfo&& ri) = default;
  template<int field_size>
  void Reduce(double* data)
  {
    int n_cliques    = cliques.size();
    int n_unfinished = n_cliques;
    for (int i = 0; i < n_cliques; ++i)
    {
      cliques[i].asyncreduce_initialize<field_size>(data);
      unfinished[i] = i;
    }
    while (n_unfinished > 0)
    {
      for (int i = 0; i < n_unfinished; ++i)
      {
        int clique_id = unfinished[i];
        // If the clique finalizes successfully, remove it from the list of cliques to finalize.
        if (cliques[clique_id].asyncreduce_finalize<field_size>(data))
        {
          unfinished[i] = unfinished[n_unfinished - 1];
          n_unfinished--;
          i--;
        }
      }
    }
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
        // std::cout << "Bad field size" << std::endl;
        throw std::invalid_argument("Bad field size of " + fs);
    }
  }
};

constexpr int CliqueSizesTag                 = 1;
constexpr int DistributeCliqueInformationTag = 2;
constexpr int FindMyCliquesTag               = 3;
length_v CliqueSizes(const globalid_v& global_node_ids,
                     globalid max_node_id,
                     const mpicontext& context)
{
  // EnsureCheckpoint(context, "Entered CliqueSizes with max_node_id " +
  // std::to_string(max_node_id)); EnsureCheckpoint(context, "Calling Creating send buffer");
  length_v sizes_send(max_node_id + 2, 0);
  for (globalid i : global_node_ids)
  {
    if (i <= max_node_id && i >= 0)
      sizes_send[i] = 1;
    else
    {
      // printf("DED");
      throw new std::logic_error("Something went horribly wrong in CliqueSizes");
    }
  }
  // EnsureCheckpoint(context, "Creating recv buffer");
  length_v sizes_recv(max_node_id + 2, 0);

  MPI_Allreduce(
      sizes_send.data(), sizes_recv.data(), max_node_id + 1, MPI_INT, MPI_SUM, context.get_comm());
  // EnsureCheckpoint(context, "Exiting CliqueSizes");
  return sizes_recv;
}

pair<color_t, vector<pair<color_t, globalid_v>>> DistributeCliqueInformation(
    const vector<pair<rank_v, globalid_v>>& clique_info,
    const mpicontext& context)
{
  pair<color_t, vector<pair<color_t, globalid_v>>> output;
  vector<pair<color_t, globalid_v>>& my_clique_info = output.second;
  if (context.get_rank() == 0)
  {
    vector<color_v> clique_assignments(context.get_size());
    int n_cliques = clique_info.size();
    for (color_t color = 0; color < n_cliques; ++color)
    {
      const rank_v& ranks_in_comm = clique_info[color].first;
      for (rankid r : ranks_in_comm)
      {
        clique_assignments[r].push_back(color);
      }
    }
    // puts("Calculated clique_assignments");
    vector<int> send_buffer;
    pair<color_t, vector<pair<color_t, const globalid_v*>>> thing_to_send;
    thing_to_send.first                                        = clique_info.size();
    vector<pair<color_t, const globalid_v*>>& rank_clique_info = thing_to_send.second;
    for (rankid r = 1; r < context.get_size(); ++r)
    {
      color_v& assignment = clique_assignments[r];
      rank_clique_info.clear();
      rank_clique_info.reserve(assignment.size());
      for (color_t color : assignment)
        rank_clique_info.emplace_back(color, &(clique_info[color].second));

      serialization::pack(thing_to_send, send_buffer);
      // printf("DistributeCliqueInformation: Sending %d ints to rank %d\n", send_buffer.size(), r);
      MPI_Send(send_buffer.data(),
               send_buffer.size(),
               MPI_INT,
               r,
               DistributeCliqueInformationTag,
               context.get_comm());
      // printf("DistributeCliqueInformation: Sent to rank %d\n", r);
    }
    {
      output.first        = clique_info.size();
      color_v& assignment = clique_assignments[0];
      my_clique_info.reserve(assignment.size());
      for (color_t color : assignment)
        my_clique_info.emplace_back(color, clique_info[color].second);
    }
  }
  else
  {
    vector<int> recv_buffer;
    MPI_Status recv_status;
    int recv_size;
    // printf("DistributeCliqueInformation rank %d: probing size\n", context.get_rank());
    MPI_Probe(0, DistributeCliqueInformationTag, context.get_comm(), &recv_status);
    MPI_Get_count(&recv_status, MPI_INT, &recv_size);
    // printf("DistributeCliqueInformation rank %d: preparing to recieve %d ints\n",
    // context.get_rank(), recv_size);
    recv_buffer.resize(recv_size);
    MPI_Recv(recv_buffer.data(),
             recv_size,
             MPI_INT,
             0,
             DistributeCliqueInformationTag,
             context.get_comm(),
             MPI_STATUS_IGNORE);
    serialization::unpack(output, recv_buffer);
  }
  // printf("rank %d clique info: ", context.get_rank());
  // std::cout << output << std::endl;
  // printf("DistributeCliqueInformation rank %d: returning my_clique_info\n", context.get_rank());
  return output;
}
pair<color_t, vector<pair<color_t, globalid_v>>> FindMyCliques(
    const std::vector<int>& global_node_ids,
    std::unique_ptr<int[]>& node_counts,
    const mpicontext& context)
{
  using std::pair;
  using std::vector;
  typedef std::vector<int> ranklst;
  typedef std::vector<int> nodelst;

  // EnsureCheckpoint(context, "Entered FindMyCliques");

  int requisite_buffer_size =
      *std::max_element(node_counts.get(), node_counts.get() + context.get_size());

  int max_node_id = GetMaximumNodeId(global_node_ids, context);
  // EnsureCheckpoint(context, "Attained Maximum node id of " + std::to_string(max_node_id));

  // The ith element of ranks_per_node is the number of ranks that share node i
  auto ranks_per_node = CliqueSizes(global_node_ids, max_node_id, context);
  // EnsureCheckpoint(context, "Attained clique sizes");
  if (context.is_root())
  {
    vector<pair<ranklst, int>> ranks_sharing_each_node(max_node_id);
    for (int i = 0; i < max_node_id; ++i)
    {
      if (ranks_per_node[i] > 1)
        ranks_sharing_each_node[i].first.reserve(ranks_per_node[i]);
      ranks_sharing_each_node[i].second = i;
    }
    for (int i : global_node_ids)
      if (ranks_per_node[i] > 1)
        ranks_sharing_each_node[i].first.push_back(0);

    // puts("Allocated ranks_sharing_each_node");
    std::unique_ptr<int[]> gni_other{new int[requisite_buffer_size]};
    for (int rank = 1; rank < context.get_size(); ++rank)
    {
      int size_of_recv = node_counts[rank];
      MPI_Recv(gni_other.get(),
               size_of_recv,
               MPI_INT,
               rank,
               FindMyCliquesTag,
               context.get_comm(),
               MPI_STATUS_IGNORE);
      for (int i = 0; i < size_of_recv; ++i)
      {
        int node_id = gni_other[i];
        if (ranks_per_node[node_id] > 1)
          ranks_sharing_each_node[node_id].first.push_back(rank);
      }
      // MPI_Barrier(context.get_comm());
    }
    auto check_if_empty = [](const pair<rank_v, int>& p) { return p.first.size() == 0; };
    // Removes any ranklst/node_id pairs where the ranklst is empty
    auto new_end = std::remove_if(
        ranks_sharing_each_node.begin(), ranks_sharing_each_node.end(), check_if_empty);

    ranks_sharing_each_node.resize(std::distance(ranks_sharing_each_node.begin(), new_end));
    // Sorts so that identical ranklsts are next to each other.
    // We're gonna perform a gather operation to organize things into cliques.
    std::sort(ranks_sharing_each_node.begin(), ranks_sharing_each_node.end());
    vector<pair<rank_v, globalid_v>> all_cliques = GroupConsecutive(ranks_sharing_each_node);
    // std::cout << "all cliques: " << all_cliques << std::endl;
    // puts("Distributing clique info");
    return DistributeCliqueInformation(all_cliques, context);
  }
  else
  {
    for (int rank = 1; rank < context.get_size(); ++rank)
    {
      if (rank == context.get_rank())
      {
        int size_of_send = global_node_ids.size();
        MPI_Send(
            global_node_ids.data(), size_of_send, MPI_INT, 0, FindMyCliquesTag, context.get_comm());
      }
      // MPI_Barrier(context.get_comm());
    }
    // Non-root ranks don't need to fill in this datastructure
    vector<pair<rank_v, globalid_v>> all_cliques_empty = {};
    return DistributeCliqueInformation(all_cliques_empty, context);
  }
}

/*std::unique_ptr<int[]> GetNumberOfNodesAssignedToEachRank(const std::vector<int>& global_node_ids,
const mpicontext& context) { int n_nodes_i_own { (int)global_node_ids.size() };
  std::unique_ptr<int[]> buffer { new int[context.get_size()] };
  MPI_Allgather(&n_nodes_i_own, 1, MPI_INT, buffer.get(), 1, MPI_INT, context.get_comm());
  return buffer;
}*/
ReductionInfoBase* GenerateReductionInfo(const std::vector<int>& insane_global_node_ids,
                                         const mpicontext& context)
{
  // EnsureCheckpoint(context, "Entering GenerateReductionInfo");
  std::unique_ptr<int[]> node_counts =
      GetNumberOfNodesAssignedToEachRank(insane_global_node_ids, context);
  vector<int> global_node_ids =
      PackIDSpace(insane_global_node_ids,
                  *std::max_element(node_counts.get(), node_counts.get() + context.get_size()),
                  context);
  // EnsureCheckpoint(context, "New sane global_node_ids: " + ss(global_node_ids));
  // EnsureCheckpoint(context, "Attained sane global_node_ids");

  std::map<globalid, localid> localid_lookup;
  for (int i = 0, max = global_node_ids.size(); i < max; ++i)
  {
    localid_lookup[global_node_ids[i]] = i;
  }
  auto my_clique_info = FindMyCliques(global_node_ids, node_counts, context);
  // EnsureCheckpoint(context, "Exited FindMyCliques");
  // Convert from global ids to local ids
  for (pair<color_t, globalid_v>& p : my_clique_info.second)
  {
    for (globalid& i : p.second)
    {
      i = localid_lookup[i];
    }
  }
  // EnsureCheckpoint(context, "reating pointer to new ReductionInfo");
  auto* reduction_info = new ReductionInfo{my_clique_info, context};
  // EnsureCheckpoint(context, "Generated ReductionInfo");
  return reduction_info;
}
}   // namespace reduction_v5
}   // namespace nimble
#endif

#endif   // NIMBLE_MPI_REDUCTION_V5_H
