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

#ifndef NIMBLE_MPI_RANK_CLIQUE_REDUCER
#define NIMBLE_MPI_RANK_CLIQUE_REDUCER

#include <mpi.h>

#include <algorithm>
#include <memory>
#include <vector>
#include <stdexcept>

#include "nimble.mpi.mpicontext.h"
#include "nimble.mpi.reduction_utils.h"
#include "nimble.quanta.arrayview.h"

namespace nimble
{
struct ReductionClique_t
{
  std::unique_ptr<int[]> indices;
  std::unique_ptr<double[]> sendbuffer;
  std::unique_ptr<double[]> recvbuffer;
  int buffersize;
  int n_indices;
  MPI_Comm clique_comm;
  MPI_Request Iallreduce_request;
  bool exists_active_asyncreduce_request = false;

  ReductionClique_t() {}
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
  ReductionClique_t(const std::vector<int>& index_list,
                    int buffersize,
                    int comm_color,
                    const mpicontext& context) :
    indices{new int[index_list.size()]},
    sendbuffer{new double[buffersize]},
    recvbuffer{new double[buffersize]},
    buffersize{buffersize},
    n_indices{(int)index_list.size()},
    // Actual instantiation of clique_comm is performed by MPI_Comm_split
    clique_comm{}
  {
    std::copy_n(index_list.data(), index_list.size(), indices.get());
    // printf("rank %d: Splitting %d with comm color %d\n", context.get_rank(), context.get_comm(),
    // comm_color);
    MPI_Comm_split(context.get_comm(), comm_color, context.get_rank(), &clique_comm);
  }
  ReductionClique_t(const std::vector<int>& index_list, int buffersize, MPI_Comm clique_comm) :
    indices{new int[index_list.size()]},
    sendbuffer{new double[buffersize]},
    recvbuffer{new double[buffersize]},
    buffersize{buffersize},
    n_indices{(int)index_list.size()},
    clique_comm{clique_comm}
  {
    std::copy_n(index_list.data(), index_list.size(), indices.get());
  }
  ReductionClique_t(ReductionClique_t&& source) = default;
  ReductionClique_t& operator=(ReductionClique_t&& source) = default;
  bool okayfieldsizeQ(int field_size) { return field_size * n_indices <= buffersize; }
  void fitnewfieldsize(int new_field_size) { resizebuffer(n_indices * new_field_size); }
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
    int *index_ptr = indices.get(), *index_ptr_end = index_ptr + n_indices;
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
    int *index_ptr = indices.get(), *index_ptr_end = index_ptr + n_indices;
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
  template<int field_size, class Lookup>
  void pack(Lookup& source)
  {
    double* destscan = sendbuffer.get();

    for (int index : quanta::arrayview_t<int>(indices.get(), n_indices))
    {
      for (int j = 0; j < field_size; ++j)
      {
        *destscan++ = source(index, j);
      }
    }
  }
  template<int field_size, class Lookup>
  void unpack(Lookup& dest)
  {
    double* sourcescan = recvbuffer.get();
    for (int index : quanta::arrayview_t<int>(indices.get(), n_indices))
    {
      for (int j = 0; j < field_size; ++j)
      {
        dest(index, j) = *sourcescan++;
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
  template<int field_size, class Lookup>
  void reduce(Lookup&& data)
  {
    EnsureDataSafety(field_size);
    pack<field_size>(data);
    MPI_Allreduce(sendbuffer.get(),
                  recvbuffer.get(),
                  n_indices * field_size,
                  MPI_DOUBLE,
                  MPI_SUM,
                  clique_comm);
    unpack<field_size>(data);
  }
  // Copies data to sendbuffer and starts an asynchronous allreduce operation.
  // asyncreduce_finalize must be called for the data to be copied from the
  // recvbuffer to the databuffer
  template<int field_size, class Lookup>
  void asyncreduce_initialize(Lookup&& databuffer)
  {
    if (exists_active_asyncreduce_request)
      throw std::logic_error(
          "asyncreduce_initialize(data) was called "
          "when an active asynchronous reduce already exists");

    EnsureDataSafety(field_size);
    pack<field_size>(databuffer);
    MPI_Iallreduce(sendbuffer.get(),
                   recvbuffer.get(),
                   n_indices * field_size,
                   MPI_DOUBLE,
                   MPI_SUM,
                   clique_comm,
                   &Iallreduce_request);
    exists_active_asyncreduce_request = true;
  }
  // Returns true if the currently active asynchronous reduce request has
  // completed Returns false if the currently active asynchronous reduce
  // request hasn't completed Throws an exception if there's no currently
  // active asynchronous reduce request
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
  // Returns true and unpacks data if the currently active asynchronous reduce
  // request has completed Returns false if the currently active asynchronous
  // reduce request hasn't completed Throws an exception if there's no
  // currently active asynchronous reduce request
  template<int field_size, class Lookup>
  bool asyncreduce_finalize(Lookup&& databuffer)
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
  int GetNumIndices()
  {
    return n_indices;
  }
  int const * GetIndices()
  {
    return indices.get();
  }
  std::vector<int> GetMPICommWorldRanks()
  {
    int my_mpi_comm_world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_comm_world_rank);
    int clique_comm_size;
    MPI_Comm_size(clique_comm, &clique_comm_size);
    std::vector<int> mpi_comm_world_ranks(clique_comm_size);
    MPI_Allgather(&my_mpi_comm_world_rank, 1, MPI_INT, mpi_comm_world_ranks.data(), 1, MPI_INT, clique_comm);
    std::sort(mpi_comm_world_ranks.begin(), mpi_comm_world_ranks.end());
    return mpi_comm_world_ranks;
  }
};
}   // namespace nimble

#endif // NIMBLE_MPI_RANK_CLIQUE_REDUCER
