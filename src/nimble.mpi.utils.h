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

#ifndef NIMBLE_MPI_UTILS_H
#define NIMBLE_MPI_UTILS_H

#ifdef NIMBLE_HAVE_MPI

#include <mpi.h>
#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "nimble.mpi.reduction.h"

namespace nimble
{

class MPIContainer
{
  std::unique_ptr<reduction::ReductionInfo> MeshReductionInfo;

 public:
  MPIContainer() {}

  // The goal of this function is to set up all the bookkeeping so that the VectorReduction()
  // function can be performed as quickly as possible
  // We'll need, for example, a list of the nodes that are shared and the rank(s) that they're
  // shared with Then, in VectorReduction, we'll send arrays of data to/from ranks that share nodes
  // To start with, each rank has a list of its global nodes (this is passed in as global_node_ids)
  void Initialize(std::vector<int> const& global_node_ids)
  {
    MPI_Comm duplicate_of_world;
    MPI_Comm_dup(MPI_COMM_WORLD, &duplicate_of_world);
    mpicontext context{duplicate_of_world};
    reduction::ReductionInfo* reduction_info = reduction::GenerateReductionInfo(global_node_ids, context);
    MeshReductionInfo.reset(reduction_info);
  }
  void GetPartitionBoundaryNodeLocalIds(std::vector<int>& node_local_ids,
                                        std::vector<int>& min_rank_containing_node)
  {
    MeshReductionInfo->GetAllIndices(node_local_ids, min_rank_containing_node);
  }
  void VectorReduction(int data_dimension, double* data)
  {
    MeshReductionInfo->PerformReduction(data, data_dimension);
  }
  template<class Lookup>
  void VectorReduction(int data_dimension, Lookup&& lookup)
  {
    MeshReductionInfo->PerformReduction(lookup, data_dimension);
  }
};
}   // namespace nimble

#else // NIMBLE_HAVE_MPI

namespace nimble
{
/// \class MPIContainer
/// \brief Empty class of communication container when NimbleSM is compiled without MPI
/// \note This empty class allows to write one code (while being a minor misnomer).
class MPIContainer
{
 public:
  MPIContainer() {}

  void Initialize(std::vector<int> const& global_node_ids) {}

  void VectorReduction(int data_dimension, double* data) {}

  template<class Lookup>
  void VectorReduction(int data_dimension, Lookup& lookup) {}
};
}   // namespace nimble

#endif // NIMBLE_HAVE_MPI
#endif // NIMBLE_MPI_UTILS_H
