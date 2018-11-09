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

#include "nimble.mpi.reduction.h"

nimble::reduction::ReductionInfo* nimble::reduction::GenerateReductionInfo(const std::vector<int>& raw_global_ids,
                                                                           const mpicontext& context)
{
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
      num_cliques = fill_clique_lookup(ids_by_rank, quanta::make_indexer(clique_lookup));
      quanta::remap(all_global_ids, quanta::make_indexer(clique_lookup));
    };

    auto find_cliques_with_vector = [&](int min_id, int max_id) {
      std::vector<int> clique_lookup(max_id - min_id + 1, 0);
      num_cliques = fill_clique_lookup(ids_by_rank, quanta::make_indexer(clique_lookup.data() - min_id));
      quanta::remap(all_global_ids, quanta::make_indexer(clique_lookup.data() - min_id));
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
