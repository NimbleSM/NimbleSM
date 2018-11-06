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

#ifdef NIMBLE_HAVE_MPI
#pragma once
#include <mpi.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "nimble.mpi.rank_clique_reducer.h"
#include "nimble.mpi.reduction.reduction_base.h"
#include "nimble.mpi.reduction_utils.h"
#include "nimble.quanta.arrayview.h"
#include "nimble.quanta.cc"
#include "nimble.quanta.functional.cc"
#include "nimble.quanta.out.h"
#include "nimble.quanta.stopwatch.h"
#include "nimble_id_pair.h"

#ifndef NIMBLE_HAVE_DARMA
#include <vector>
#endif

namespace nimble
{
namespace reduction_v6
{
struct ReductionInfo : public ReductionInfoBase
{
  template<class... Args>
  using vector = std::vector<Args...>;
  template<class... Args>
  using hashmap = std::unordered_map<Args...>;
  std::vector<ReductionClique_t> cliques;
  std::vector<int> unfinished;

  ReductionInfo(const vector<int>& clique_colors,
                const vector<int>& clique_ids,
                const vector<int>& clique_assignment,
                const vector<int>& global_ids,
                const mpicontext& context)
  {
    hashmap<int, vector<int>> clique_index_assignment;
    clique_index_assignment.reserve(clique_ids.size());
    int num_ranks = context.get_size();

    hashmap<int, int> counts;
    counts.reserve(clique_ids.size());

    for (int clique_id : clique_assignment)
      if (clique_id > num_ranks)
        counts[clique_id] += 1;

    for (int clique_id : clique_ids)
      clique_index_assignment[clique_id].reserve(counts[clique_id]);

    for (int index = 0, max = clique_assignment.size(); index < max; ++index)
    {
      int clique_id = clique_assignment[index];
      if (clique_id > num_ranks)
      {
        clique_index_assignment[clique_id].emplace_back(index);
      }
    }

    std::vector<MPI_Comm> comms{};
    comms.reserve(clique_ids.size());

    {
      quanta::stopwatch s;
      for (int color : clique_colors)
      {
        MPI_Comm comm = context.split_by_color(color);
        if (comm != MPI_COMM_NULL)
        {
          comms.push_back(comm);
        }
      }
      double _time = s.age();
      /* context.print_formatted(std::to_string(_time), */
      /*                         "\"reduction_v6 time to create comms\": [", */
      /*                         ", ", */
      /*                         "],\n", */
      /*                         std::cerr); */
    }

    if (comms.size() != clique_ids.size())
      throw std::logic_error("**** Error, comms.size() != clique_ids.size().");

    for (size_t i = 0; i < clique_ids.size(); ++i)
    {
      int clique_id    = clique_ids[i];
      MPI_Comm comm    = comms[i];
      auto& index_list = clique_index_assignment[clique_id];
      std::sort(index_list.begin(), index_list.end(), [&](int a, int b) {
        return global_ids[a] < global_ids[b];
      });
      this->cliques.emplace_back(std::move(index_list), index_list.size() * 3, comm);

      // DEBUGGING
      int my_mpi_comm_world_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_comm_world_rank);
      std::cout << "DEBUGGING rank " << my_mpi_comm_world_rank << " is adding a clique with index_list.size() " << index_list.size() << std::endl;
      // DEBUGGING
    }
  }

  ReductionInfo() {}
  ReductionInfo(ReductionInfo&& ri) = default;
  ReductionInfo& operator=(ReductionInfo&& ri) = default;
  template<int field_size, class Lookup>
  void Reduce(Lookup&& data)
  {
    for (auto& clique : cliques)
      clique.asyncreduce_initialize<field_size>(data);

    unfinished.resize(cliques.size());
    std::iota(unfinished.begin(), unfinished.end(), 0);

    while (!unfinished.empty())
    {
      int increment = 0;
      for (size_t i = 0; i < unfinished.size(); i += increment)
      {
        bool reduceFinished = cliques[unfinished[i]].asyncreduce_finalize<field_size>(data);

        increment = !reduceFinished;

        if (reduceFinished)
        {
          unfinished[i] = unfinished.back();
          unfinished.pop_back();
        }
      }
    }
  }
  std::vector<int> GetAllIndices()
  {
    std::set<int> index_set;
    for (auto& clique : cliques)
    {
      // DJL todo, expand this function to also return the mpi ranks associated with each shared node
      // std::vector<int> temp = clique.GetMPICommWorldRanks();
      // int min_mpi_comm_world_rank = clique.GetMPICommWorldRanks()[0];

      int num_indices = clique.GetNumIndices();
      int const * clique_indices = clique.GetIndices();
      for (int i=0 ; i<num_indices ; i++)
      {
        index_set.insert(clique_indices[i]);
      }
    }
    std::vector<int> all_indices(index_set.begin(), index_set.end());
    return all_indices;
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
  template<class Lookup>
  void PerformReduction(Lookup& lookup, int field_size) {
    switch (field_size)
    {
      case 1: Reduce<1>(lookup); break;
      case 2: Reduce<2>(lookup); break;
      case 3: Reduce<3>(lookup); break;
      default:
        std::string fs = std::to_string(field_size);
        // std::cout << "Bad field size" << std::endl;
        throw std::invalid_argument("Bad field size of " + fs);
    }
  }
};

template<class list_of_lists_t, class F>
auto fill_clique_lookup(list_of_lists_t& ids_by_rank, F&& clique_lookup) ->
    typename std::decay<quanta::transformed_iterated_t<F, quanta::elem_t<list_of_lists_t>>>::type
{
  typedef decltype(clique_lookup) lookup_t;
  typedef quanta::elem_t<list_of_lists_t> inner_list_t;
  typedef
      typename std::decay<quanta::transformed_iterated_t<lookup_t, inner_list_t>>::type clique_t;

  std::unordered_map<clique_t, clique_t> remapped_cliques{};

  const clique_t zero{};
  clique_t rank_plus_one{};

  auto generate_clique_id = quanta::make_counter<clique_t>(quanta::len(ids_by_rank) + 1);

  for (auto& id_list : ids_by_rank)
  {
    remapped_cliques.clear();
    rank_plus_one++;

    for (auto& id : id_list)
    {
      auto& current_clique = clique_lookup(id);

      if (current_clique == zero)
      {
        current_clique = rank_plus_one;
      }
      else
      {
        auto remapped_clique_iter = remapped_cliques.find(current_clique);
        if (remapped_clique_iter != remapped_cliques.end())
        {
          current_clique = remapped_clique_iter->second;
        }
        else
        {
          auto new_clique                  = generate_clique_id();
          remapped_cliques[current_clique] = new_clique;
          current_clique                   = new_clique;
        }
      }
    }
  }
  return generate_clique_id.get_count();
}
ReductionInfoBase* GenerateReductionInfo(const std::vector<int>& raw_global_ids,
                                         const mpicontext& context)
{
  using namespace quanta;
  using std::vector;
  int rank           = context.get_rank();
  int numRanks       = context.get_size();
  size_t nodeIDCount = raw_global_ids.size();

  std::vector<int> IDCounts = context.allgather_to_vector((int)nodeIDCount);

  int maxIDCount         = 0;
  int IDCountAccumulator = 0;
  std::vector<int> displacements;
  displacements.reserve(IDCounts.size());
  for (int i : IDCounts)
  {
    displacements.emplace_back(IDCountAccumulator);
    IDCountAccumulator += i;
    if (i > maxIDCount)
      maxIDCount = i;
  }

  vector<int> clique_assignment(raw_global_ids.size());
  int num_cliques;

  if (context.is_root())
  {
    std::vector<int> all_global_ids(IDCountAccumulator);
    // See: https://www.open-mpi.org/doc/v2.1/man3/MPI_Gatherv.3.php
    context.gatherv_recieve(raw_global_ids, all_global_ids, IDCounts, displacements);

    auto ids_by_rank = quanta::partition_into_arrayviews(all_global_ids, IDCounts);

    auto find_cliques_with_unordered_map = [&] {
      std::unordered_map<int, int> clique_lookup;
      clique_lookup.reserve(all_global_ids.size());
      num_cliques = fill_clique_lookup(ids_by_rank, make_indexer(clique_lookup));
      remap(all_global_ids, make_indexer(clique_lookup));
    };

    auto find_cliques_with_vector = [&](int min_id, int max_id) {
      std::vector<int> clique_lookup(max_id - min_id + 1, 0);
      num_cliques = fill_clique_lookup(ids_by_rank, make_indexer(clique_lookup.data() - min_id));
      remap(all_global_ids, make_indexer(clique_lookup.data() - min_id));
    };

    // Use the unordered map version for now
    find_cliques_with_unordered_map();
    /*
    At this point, every global id has been replaced with a clique id. If the clique id is smaller
    than or equal to the number of ranks, then that clique id correspons to a clique with only one
    rank as a member (that rank being one less than the clique id). Otherwise, the clique id
    must correspond to a clique that belongs to multiple ranks.
    */

    context.scatterv_send(all_global_ids, IDCounts, displacements, clique_assignment);
  }
  else
  {
    context.gatherv_send(raw_global_ids);

    context.scatterv_recieve(clique_assignment);
  }

  context.bcast(num_cliques);

  ReductionInfo* reduction_info;
  {
    std::vector<int> clique_ids = clique_assignment;
    std::sort(clique_ids.begin(), clique_ids.end());
    auto new_end_iter = std::unique(clique_ids.begin(), clique_ids.end());
    clique_ids.resize(std::distance(clique_ids.begin(), new_end_iter));
    std::vector<int> clique_colors(num_cliques, MPI_UNDEFINED);
    for (int i : clique_ids)
      clique_colors[i] = 1;
    reduction_info =
        new ReductionInfo(clique_colors, clique_ids, clique_assignment, raw_global_ids, context);
  }
  return reduction_info;
}
}   // namespace reduction_v6
}   // namespace nimble
#endif
