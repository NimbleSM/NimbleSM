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

#ifndef NIMBLE_MPI_REDUCTION_H
#define NIMBLE_MPI_REDUCTION_H

#ifdef NIMBLE_HAVE_MPI
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
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "nimble.mpi.rank_clique_reducer.h"
#include "nimble.mpi.reduction_utils.h"
#include "nimble.quanta.h"

namespace nimble
{
namespace reduction
{
struct ReductionInfo
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
      for (int color : clique_colors)
      {
        MPI_Comm comm = context.split_by_color(color);
        if (comm != MPI_COMM_NULL)
        {
          comms.push_back(comm);
        }
      }
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
  void GetAllIndices(std::vector<int>& indices,
                     std::vector<int>& min_rank_containing_index)
  {
    for (auto& clique : cliques)
    {
      int min_mpi_comm_world_rank = clique.GetMPICommWorldRanks()[0];
      int num_indices = clique.GetNumIndices();
      int const * clique_indices = clique.GetIndices();
      for (int i=0 ; i<num_indices ; i++)
      {
        indices.push_back(clique_indices[i]);
        min_rank_containing_index.push_back(min_mpi_comm_world_rank);
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

ReductionInfo* GenerateReductionInfo(const std::vector<int>& raw_global_ids,
                                     const mpicontext& context);

}   // namespace reduction
}   // namespace nimble
#endif

#endif // NIMBLE_MPI_REDUCTION_H
