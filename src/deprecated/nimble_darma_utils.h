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

#ifndef NIMBLE_DARMA_UTILS_H
#define NIMBLE_DARMA_UTILS_H

#include <darma.h>
#include <darma/impl/array/index_range.h>
#include <darma/impl/task_collection/access_handle_collection.h>
#include <darma/impl/task_collection/task_collection.h>
#include <darma/interface/app/backend_hint.h>
#include <darma/serialization/serializers/all.h>

#include "nimble_boundary_condition_manager.h"
#include "nimble_data_manager.h"
#include "nimble_exodus_output.h"
#include "nimble_genesis_mesh.h"
#include "nimble_parser.h"

namespace nimble {

enum ProgressBarFlag
{
  NONE = 0,
  FIRST_STEP,
  PRINT_PROGRESS,
  LAST_STEP
};

// todo: this is needed only because I can't pass a std::string as a functor
// argument
enum DistributedVectorReductionQuantity
{
  NO_QUANTITY = 0,
  LUMPED_MASS,
  INTERNAL_FORCE
};

void
PrintProgress(ProgressBarFlag progress_bar_flag, int step, int num_steps);

typedef std::map<int, std::vector<std::vector<double>>> DerivedElementDataType;

struct ReductionMin
{
  static constexpr auto is_indexed = true;
  template <typename ArgT>
  void
  reduce(ArgT const& src, ArgT& dest) const
  {
    dest = src < dest ? src : dest;
  }
};

struct ReductionMax
{
  static constexpr auto is_indexed = true;
  template <typename ArgT>
  void
  reduce(ArgT const& src, ArgT& dest) const
  {
    dest = src > dest ? src : dest;
  }
};

struct ReadGenesisFiles
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> macroscale_mesh_collection,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> rve_mesh_collection,
      Parser                                                                          parser) const;
};

struct InitializeDataManager
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> macroscale_mesh_collection,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> rve_mesh_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection,
      Parser                                                                          parser) const;
};

struct IdentifyGloballySharedNodes
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> mesh_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>>
          boundary_node_global_ids_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>> boundary_ranks_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>>
          my_global_node_ids_collection) const;
};

struct InitializeBoundaryConditionManager
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> mesh_collection,
      darma_runtime::AccessHandleCollection<BoundaryConditionManager, darma_runtime::Range1D<int>>
             boundary_condition_manager_collection,
      Parser parser) const;
};

struct InitializeExodusOutput
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                      index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>>  mesh_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>>  data_manager_collection,
      darma_runtime::AccessHandleCollection<ExodusOutput, darma_runtime::Range1D<int>> exodus_output_collection,
      Parser                                                                           parser) const;
};

struct ComputeLumpedMass
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> macroscale_mesh_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection) const;
};

struct DistributedVectorReduction
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>>
          boundary_node_global_ids_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>> boundary_ranks_collection,
      darma_runtime::AccessHandleCollection<std::vector<double>, darma_runtime::Range1D<int>>
                                         global_reduction_buffer_collection,
      DistributedVectorReductionQuantity distributed_vector_reduction_quantity,
      int                                step) const;
};

struct ComputeCriticalTimeStep
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> macroscale_mesh_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection,
      Parser                                                                          parser) const;
};

struct ApplyInitialConditions
{
  void
  operator()(
      darma_runtime::Index1D<int> index,
      darma_runtime::AccessHandleCollection<BoundaryConditionManager, darma_runtime::Range1D<int>>
          boundary_condition_manager_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection) const;
};

struct ApplyKinematicBC
{
  void
  operator()(
      darma_runtime::Index1D<int> index,
      darma_runtime::AccessHandleCollection<BoundaryConditionManager, darma_runtime::Range1D<int>>
          boundary_condition_manager_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection,
      double                                                                          time_current,
      double                                                                          time_previous) const;
};

struct ComputeDerivedElementData
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> macroscale_mesh_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection,
      darma_runtime::AccessHandleCollection<DerivedElementDataType, darma_runtime::Range1D<int>>
          derived_element_data_collection) const;
};

struct ExodusWriteStep
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                      index,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>>  data_manager_collection,
      darma_runtime::AccessHandleCollection<ExodusOutput, darma_runtime::Range1D<int>> exodus_output_collection,
      darma_runtime::AccessHandleCollection<DerivedElementDataType, darma_runtime::Range1D<int>>
             derived_element_data_collection,
      double time_current) const;
};

struct Initialize
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> macroscale_mesh_collection,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> rve_mesh_collection,
      darma_runtime::AccessHandleCollection<BoundaryConditionManager, darma_runtime::Range1D<int>>
          boundary_condition_manager_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>>
          boundary_node_global_ids_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>> boundary_ranks_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>>
                                                                                       my_global_node_ids_collection,
      darma_runtime::AccessHandleCollection<ExodusOutput, darma_runtime::Range1D<int>> exodus_output_collection,
      Parser                                                                           parser) const;
};

struct ExplicitTimeStep
{
  void
  operator()(
      darma_runtime::Index1D<int>                                                     index,
      darma_runtime::AccessHandleCollection<GenesisMesh, darma_runtime::Range1D<int>> macroscale_mesh_collection,
      darma_runtime::AccessHandleCollection<BoundaryConditionManager, darma_runtime::Range1D<int>>
          boundary_condition_manager_collection,
      darma_runtime::AccessHandleCollection<DataManager, darma_runtime::Range1D<int>> data_manager_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>>
          boundary_node_global_ids_collection,
      darma_runtime::AccessHandleCollection<std::vector<int>, darma_runtime::Range1D<int>> boundary_ranks_collection,
      darma_runtime::AccessHandleCollection<std::vector<double>, darma_runtime::Range1D<int>>
                      global_reduction_buffer_collection,
      int             step,
      int             num_steps,
      double          time_current,
      double          time_previous,
      ProgressBarFlag progress_bar_flag,
      bool            is_output_step) const;
};

}  // namespace nimble

#endif
