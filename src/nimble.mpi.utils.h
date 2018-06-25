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

#include "nimble_data_manager.h"
#include "nimble_genesis_mesh.h"

#ifdef NIMBLE_HAVE_MPI
#include <mpi.h>
#endif

namespace nimble
{
void PrintReductionTimingInfo();

}   // namespace nimble
#endif

#ifdef NIMBLE_HAVE_MPI
#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "nimble.mpi.reduction.reduction_base.h"
#include "nimble.mpi.reduction.v1.h"
#include "nimble.mpi.reduction.v2.h"
#include "nimble.mpi.reduction.v3.h"
#include "nimble.mpi.reduction.v4.h"
#include "nimble.mpi.reduction.v5.h"
#include "nimble.mpi.reduction.v6.h"
#include "nimble.quanta.stopwatch.h"
namespace nimble
{
std::pair<int, double> ReductionTimingInfo;
void PrintReductionTimingInfo()
{
  std::pair<int, double> timinginfo = ReductionTimingInfo;
  {
    int flag{};
    MPI_Initialized(&flag);
    if (flag == 0)
      return;
  }
  int reduction_calls            = timinginfo.first;
  double time_spend_on_reduction = timinginfo.second;
  int world_size{}, my_mpi_rank{};
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi_rank);
  if (my_mpi_rank == 0)
  {
    std::cout << "Benchmarking information: " << std::endl;
    std::cout << "A is rank," << std::endl;
    std::cout << "B is the number of reduction calls," << std::endl;
    std::cout << "C is the total cpu time for reduction (seconds)," << std::endl;
    std::cout << "D is the mean cpu time for reduction (seconds)" << std::endl;
    std::cout << "Col A\tCol B\tCol C    \tCol D" << std::endl;
  }
  for (int i = 0; i < world_size; ++i)
  {
    if (i == my_mpi_rank)
    {
      std::cout << i << "\t" << reduction_calls << "\t" << time_spend_on_reduction << "\t"
                << time_spend_on_reduction / reduction_calls << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  double total_time_for_all_nodes{};
  MPI_Reduce(&time_spend_on_reduction,
             &total_time_for_all_nodes,
             1,
             MPI_DOUBLE,
             MPI_SUM,
             0,
             MPI_COMM_WORLD);
  if (my_mpi_rank == 0)
  {
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Info over all ranks:" << std::endl;
    std::cout << "Total cpu time spent on reduction over all nodes:\t" << total_time_for_all_nodes
              << std::endl;
    std::cout << "Average cpu time spent on reduction:\t" << total_time_for_all_nodes / world_size
              << std::endl;
    std::cout << "Average cpu time per reduction:\t"
              << total_time_for_all_nodes / (world_size * reduction_calls) << std::endl;
    std::cout << std::endl << std::endl;
  }
}
std::pair<ReductionInfoBase*, int> GetReductionInfoTypeByVersion(
    int version,
    std::vector<int> const& global_node_ids,
    const mpicontext& context)
{
  if (context.is_root())
    std::cout << "     mpi reduction version " << version << std::endl;
  if (version == 0)
  {
    // 0 means "use default version"
    // The current default version is 5
    version = 6;
  }
  context.print_if_root("\"reduction version\": ", std::cerr);
  context.print_formatted(std::to_string(version), "[", ", ", "],\n", std::cerr);
  ReductionInfoBase* reduction_method;
  quanta::stopwatch s{};
  switch (version)
  {
    case 1:
      reduction_method = reduction_v1::GenerateReductionInfo(
          global_node_ids, context.get_rank(), context.get_size());
      break;
    case 2:
      reduction_method = reduction_v2::GenerateReductionInfo(
          global_node_ids, context.get_rank(), context.get_size());
      break;
    case 3:
      reduction_method = reduction_v3::GenerateReductionInfo(
          global_node_ids, context.get_rank(), context.get_size());
      break;
    case 4:
      reduction_method = reduction_v4::GenerateReductionInfo(
          global_node_ids, context.get_rank(), context.get_size());
      break;
    case 5: reduction_method = reduction_v5::GenerateReductionInfo(global_node_ids, context); break;
    case 6: reduction_method = reduction_v6::GenerateReductionInfo(global_node_ids, context); break;
    default:
      throw new std::invalid_argument(
          "From: MPIContainer::Initialize in mimimultiscale_mpi_utils.cc:\n"
          "Error: Bad reduction version given.\n"
          "Reduction version should be an integer in the range [1,5] "
          "or it should be 0, which indicates that the reduction version "
          "is to be determined automatically.\n"
          "Recieved: " +
          std::to_string(version) + "\n");
  }
  double _elapsed = s.age();
  context.print_if_root("\"GenerateReductionInfo benchmark\": ", std::cerr);
  context.print_formatted(std::to_string(_elapsed), "[", ", ", "]\n", std::cerr);
  return {reduction_method, version};
}

class MPIContainer
{
  std::unique_ptr<ReductionInfoBase> MeshReductionInfo;

  int _reduction_version;

 public:
  MPIContainer() {}

  // The goal of this function is to set up all the bookkeeping so that the VectorReduction()
  // function can be performed as quickly as possible
  // We'll need, for example, a list of the nodes that are shared and the rank(s) that they're
  // shared with Then, in VectorReduction, we'll send arrays of data to/from ranks that share nodes
  // To start with, each rank has a list of its global nodes (this is passed in as global_node_ids)
  void Initialize(std::vector<int> const& global_node_ids, int reduction_version)
  {
    using namespace std;

    MPI_Comm duplicate_of_world;
    MPI_Comm_dup(MPI_COMM_WORLD, &duplicate_of_world);
    mpicontext context{duplicate_of_world};
    std::pair<ReductionInfoBase*, int> reduction_info =
        GetReductionInfoTypeByVersion(reduction_version, global_node_ids, context);
    MeshReductionInfo.reset(reduction_info.first);
    this->_reduction_version = reduction_info.second;
  }
  void VectorReduction(int data_dimension, double* data)
  {
    double& total_reduction_time = ReductionTimingInfo.second;
    quanta::stopwatch s;
    s.reset();
    //The below commented-out code was used to test reduction version 6 using a lookup function
    //Instead of reading directly from a pointer. 
    /*
    if(_reduction_version == 6) {
      auto lookup = [data, field_size = data_dimension](int node_id, int index) -> double& {
        return data[node_id * field_size + index];
      };
      auto& reduction_info = (reduction_v6::ReductionInfo&)*MeshReductionInfo.get();
      reduction_info.PerformReduction(lookup, data_dimension);
    }
    else MeshReductionInfo->PerformReduction(data, data_dimension);
    */
    MeshReductionInfo->PerformReduction(data, data_dimension);
    total_reduction_time += s.age();
    int& total_reduction_calls = ReductionTimingInfo.first;
    total_reduction_calls += 1;
  }
  template<class Lookup>
  void VectorReduction(int data_dimension, Lookup&& lookup)
  {
    double& total_reduction_time = ReductionTimingInfo.second;
    quanta::stopwatch s;
    s.reset();
    switch (_reduction_version)
    {
      case 6:
      {
        auto& reduction_info = (reduction_v6::ReductionInfo&)*MeshReductionInfo.get();
        reduction_info.PerformReduction(lookup, data_dimension);
        break;
      }
      default:
      {
        throw std::runtime_error(
            "Use of kokkos lookup function not supported for reduction versions earlier than "
            "version 6.");
      }
    }
    total_reduction_time += s.age();
    int& total_reduction_calls = ReductionTimingInfo.first;
    total_reduction_calls += 1;
    // To do: write vector reduction code based on lookup rather than based
    // on data pointer!!!
  }
};
}   // namespace nimble

#else
// Throws runtime error if any of these are called as
// They're not implemented when compiling without mpi.
namespace nimble
{
class MPIContainer
{
 public:
  MPIContainer() {}

  void Initialize(std::vector<int> const& global_node_ids, int reduction_version)
  {
    throw std::runtime_error(
        "[Called MPIContainer.Initialize()]: Calling MPI functions in version compiled without "
        "MPI; try compiling with "
        "-DNIMBLE_HAVE_MPI to compile with MPI");
  }

  void VectorReduction(int data_dimension, double* data)
  {
    throw std::runtime_error(
        "[Called MPIContainer.VectorReduction]: Calling MPI functions in version compiled without "
        "MPI; try compiling with "
        "-DNIMBLE_HAVE_MPI to compile with MPI");
  }

  template<class Lookup>
  void VectorReduction(int data_dimension, Lookup& lookup)
  {
    throw std::runtime_error(
        "[Called MPIContainer.VectorReduction]: Calling MPI functions in version compiled without "
        "MPI; try compiling with "
        "-DNIMBLE_HAVE_MPI to compile with MPI");
  }
};
}   // namespace nimble
#endif
