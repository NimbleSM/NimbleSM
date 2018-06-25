#pragma once
#include <mpi.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>

#include "nimble.mpi.mpicontext.h"
#include "nimble.mpi.reduction.reduction_base.h"
#include "nimble.mpi.reduction_utils.h"
#include "nimble.quanta.arrayview.h"
#include "nimble.quanta.cc"
#include "nimble.quanta.functional.cc"
#include "nimble_id_pair.h"
namespace nimble
{
struct ReductionClique_t
{
  std::unique_ptr<int[]> indicies;
  std::unique_ptr<double[]> sendbuffer;
  std::unique_ptr<double[]> recvbuffer;
  int buffersize;
  int n_indicies;
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
    indicies{new int[index_list.size()]},
    sendbuffer{new double[buffersize]},
    recvbuffer{new double[buffersize]},
    buffersize{buffersize},
    n_indicies{(int)index_list.size()},
    // Actual instantiation of clique_comm is performed by MPI_Comm_split
    clique_comm{}
  {
    std::copy_n(index_list.data(), index_list.size(), indicies.get());
    // printf("rank %d: Splitting %d with comm color %d\n", context.get_rank(), context.get_comm(),
    // comm_color);
    MPI_Comm_split(context.get_comm(), comm_color, context.get_rank(), &clique_comm);
  }
  ReductionClique_t(const std::vector<int>& index_list, int buffersize, MPI_Comm clique_comm) :
    indicies{new int[index_list.size()]},
    sendbuffer{new double[buffersize]},
    recvbuffer{new double[buffersize]},
    buffersize{buffersize},
    n_indicies{(int)index_list.size()},
    clique_comm{clique_comm}
  {
    std::copy_n(index_list.data(), index_list.size(), indicies.get());
  }
  ReductionClique_t(ReductionClique_t&& source) = default;
  ReductionClique_t& operator=(ReductionClique_t&& source) = default;
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
  template<int field_size, class Lookup>
  void pack(Lookup& source)
  {
    double* destscan = sendbuffer.get();

    for (int index : quanta::arrayview_t<int>(indicies.get(), n_indicies))
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
    for (int index : quanta::arrayview_t<int>(indicies.get(), n_indicies))
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
                  n_indicies * field_size,
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
                   n_indicies * field_size,
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
};
}   // namespace nimble