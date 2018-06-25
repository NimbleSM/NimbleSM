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

#ifndef NIMBLE_MPI_REDUCTION_V1_H
#define NIMBLE_MPI_REDUCTION_V1_H

#ifdef NIMBLE_HAVE_MPI

#include <mpi.h>

#include <algorithm>
#include <exception>
#include <iostream>

#include "nimble.mpi.reduction.reduction_base.h"

#ifndef NIMBLE_HAVE_DARMA
#include <vector>
#endif

namespace nimble
{
namespace reduction_v1
{
struct ReductionInfo : public ReductionInfoBase
{
  std::vector<int> global_node_ids;
  double* send_buffer;
  double* recv_buffer;
  int total_nodes;
  bool deleteQ;

  ReductionInfo() : deleteQ(false), send_buffer(0), recv_buffer(0) {}
  ReductionInfo(ReductionInfo&& ri) :
    send_buffer(ri.send_buffer),
    recv_buffer(ri.recv_buffer),
    total_nodes(ri.total_nodes),
    global_node_ids(std::move(ri.global_node_ids)),
    deleteQ(ri.deleteQ)
  {
    ri.deleteQ = false;
  }

  ReductionInfo(const std::vector<int>& global_node_ids, int total_nodes) :
    global_node_ids(global_node_ids),
    total_nodes(total_nodes),
    send_buffer(new double[total_nodes]{}),
    recv_buffer(new double[total_nodes]{}),
    deleteQ(true)
  {
  }

  void operator=(ReductionInfo&& ri)
  {
    if (deleteQ)
    {
      if (send_buffer)
        delete[] send_buffer;
      if (recv_buffer)
        delete[] recv_buffer;
    }
    global_node_ids = std::move(ri.global_node_ids);
    send_buffer     = ri.send_buffer;
    recv_buffer     = ri.recv_buffer;
    total_nodes     = ri.total_nodes;
    deleteQ         = ri.deleteQ;
    ri.deleteQ      = false;
  }
  ~ReductionInfo()
  {
    if (deleteQ)
    {
      if (send_buffer)
        delete[] send_buffer;
      if (recv_buffer)
        delete[] recv_buffer;
    }
  }
  void PerformReduction(double* data, int field_size)
  {
    for (int shift = 0; shift < field_size; ++shift)
    {
      auto& global_node_ids = this->global_node_ids;
      int total_nodes       = this->total_nodes;
      double *send_buf = this->send_buffer, *recv_buf = this->recv_buffer;

      for (int i = 0; i < total_nodes; ++i)
        send_buf[i] = 0;
      int n_my_nodes = global_node_ids.size();

      for (int i = 0; i < n_my_nodes; ++i)
      {
        double& data_value           = data[i * field_size + shift];
        send_buf[global_node_ids[i]] = data_value;
      }
      MPI_Allreduce(send_buf, recv_buf, total_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      for (int i = 0; i < n_my_nodes; ++i)
      {
        data[i * field_size + shift] = recv_buf[global_node_ids[i]];
      }
    }
  }
};
ReductionInfo* GenerateReductionInfo(const std::vector<int>& global_node_ids,
                                     int my_rank,
                                     int total_ranks)
{
  int num_local_ids  = static_cast<int>(global_node_ids.size());
  int num_global_ids = 0;
  MPI_Allreduce(&num_local_ids, &num_global_ids, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return new ReductionInfo(global_node_ids, num_global_ids);
}
void PerformReduction(double* data, int field_size, ReductionInfo& info)
{
  for (int shift = 0; shift < field_size; ++shift)
  {
    auto& global_node_ids = info.global_node_ids;
    int total_nodes       = info.total_nodes;
    double *send_buf = info.send_buffer, *recv_buf = info.recv_buffer;

    for (int i = 0; i < total_nodes; ++i)
      send_buf[i] = 0;
    int n_my_nodes = global_node_ids.size();

    for (int i = 0; i < n_my_nodes; ++i)
    {
      double& data_value           = data[i * field_size + shift];
      send_buf[global_node_ids[i]] = data_value;
    }
    MPI_Allreduce(send_buf, recv_buf, total_nodes, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i < n_my_nodes; ++i)
    {
      data[i * field_size + shift] = recv_buf[global_node_ids[i]];
    }
  }
}
}   // namespace reduction_v1
}   // namespace nimble
#endif

#endif   // NIMBLE_MPI_REDUCTION_V1_H
