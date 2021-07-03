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

#if defined(NIMBLE_HAVE_ARBORX) && defined(NIMBLE_HAVE_MPI)

#include "arborx_parallel_contact_manager.h"

#include "nimble_data_manager.h"
#include "nimble_defs.h"
#include "nimble_model_data_base.h"
#include "nimble_vector_communicator.h"

#include "contact/arborx_utils.h"

#ifdef NIMBLE_HAVE_KOKKOS
#include "nimble_kokkos_model_data.h"
#endif

#include <mpi.h>

#include "ArborX.hpp"
#include "ArborX_Config.hpp"

#include <Kokkos_Core.hpp>
#include <iostream>
#include <random>
#include <set>
#include <vector>

#if defined(KOKKOS_ENABLE_OPENMP)
#include <omp.h>
#endif

namespace nimble {

namespace details {

constexpr int dim = 3;

struct PredicateTypeNodesRank
{
  nimble_kokkos::DeviceContactEntityArrayView nodes_;
  int                                         rank_;
};

struct OutputData
{
  int    index_;
  int    rank_;
  double coord_[dim];
  bool   has_force_;
  double force_node_[dim];
};

struct PairData
{
  int pred_rank_  = 0;
  int pred_index_ = -1;
  int prim_rank_  = 0;
  int prim_index_ = -1;
  //
  PairData(int r_, int i_, int f_r_, int f_i_) : pred_rank_(r_), pred_index_(i_), prim_rank_(f_r_), prim_index_(f_i_) {}
  //
  PairData() = default;
  //
  ~PairData() = default;
  //
  bool
  operator<(const PairData& rhs) const
  {
    return (
        std::tie(pred_rank_, pred_index_, prim_rank_, prim_index_) <
        std::tie(rhs.pred_rank_, rhs.pred_index_, rhs.prim_rank_, rhs.prim_index_));
  }
};

struct ContactCallback
{
  const int                                    rank_;
  nimble_kokkos::DeviceContactEntityArrayView& faces_;
  const double                                 penalty_;
  std::set<PairData>&                          list_;
  //
  template <typename Predicate, typename OutputFunctor>
  KOKKOS_FUNCTION void
  operator()(Predicate const& pred, int f_primitive, OutputFunctor const& out) const
  {
    auto const& p_data     = getData(pred);      // <- type OutputData
    auto const& p_geometry = getGeometry(pred);  // <- type Box
    //
    //--- Define copy contact entity
    //--- For the node, we only store the coordinates.
    //
    ContactEntity myNode(ContactEntity::ContactEntityType::NODE, 0, p_data.coord_, 0.0, 0);
    ContactEntity myFace;
    faces_(f_primitive).ExportGeometryInto(myFace);
    //
    double gap                    = 0.0;
    double normal[dim]            = {0., 0., 0.};
    bool   inside                 = false;
    double facet_coordinates[dim] = {0., 0., 0.};
    //
    ContactManager::Projection(myNode, myFace, inside, gap, &normal[0], &facet_coordinates[0]);
    //
    double force[dim] = {0., 0., 0.};
    if (inside) {
      //
      details::getContactForce(penalty_, gap, normal, force);
      //
      bool noSkip = false;
#pragma omp critical
      {
        auto result = list_.insert(PairData{p_data.rank_, p_data.index_, rank_, f_primitive});
        noSkip      = result.second;
      }
      //
      myFace.SetNodalContactForces(force, &facet_coordinates[0]);
      //
      if (noSkip) {
        Kokkos::atomic_assign(&faces_(f_primitive).contact_status_, true);
        Kokkos::atomic_add(&faces_(f_primitive).force_1_x_, myFace.force_1_x_);
        Kokkos::atomic_add(&faces_(f_primitive).force_1_y_, myFace.force_1_y_);
        Kokkos::atomic_add(&faces_(f_primitive).force_1_z_, myFace.force_1_z_);
        Kokkos::atomic_add(&faces_(f_primitive).force_2_x_, myFace.force_2_x_);
        Kokkos::atomic_add(&faces_(f_primitive).force_2_y_, myFace.force_2_y_);
        Kokkos::atomic_add(&faces_(f_primitive).force_2_z_, myFace.force_2_z_);
        Kokkos::atomic_add(&faces_(f_primitive).force_3_x_, myFace.force_3_x_);
        Kokkos::atomic_add(&faces_(f_primitive).force_3_y_, myFace.force_3_y_);
        Kokkos::atomic_add(&faces_(f_primitive).force_3_z_, myFace.force_3_z_);
      }
    }
    //
    out(
        {f_primitive,
         rank_,
         {p_data.coord_[0], p_data.coord_[1], p_data.coord_[2]},
         inside,
         {-force[0], -force[1], -force[2]}});
  }
};

}  // namespace details
}  // namespace nimble

namespace ArborX {

template <>
struct AccessTraits<nimble::details::PredicateTypeNodesRank, PredicatesTag>
{
  static std::size_t
  size(nimble::details::PredicateTypeNodesRank const& v)
  {
    return v.nodes_.extent(0);
  }

  KOKKOS_FUNCTION static auto
  get(nimble::details::PredicateTypeNodesRank const& v, std::size_t i)
  {
    nimble::ContactEntity& e = v.nodes_(i);
    ArborX::Point          point1(e.bounding_box_x_min_, e.bounding_box_y_min_, e.bounding_box_z_min_);
    ArborX::Point          point2(e.bounding_box_x_max_, e.bounding_box_y_max_, e.bounding_box_z_max_);
    ArborX::Box            box(point1, point2);
    //
    return attach(
        intersects(box),
        nimble::details::OutputData{
            static_cast<int>(i), v.rank_, {e.coord_1_x_, e.coord_1_y_, e.coord_1_z_}, false, {0, 0, 0}});
  }
  using memory_space = nimble_kokkos::kokkos_device_memory_space;
};
}  // namespace ArborX

namespace nimble {

using memory_space = nimble_kokkos::kokkos_device_memory_space;
using kokkos_device =
    Kokkos::Device<nimble_kokkos::kokkos_device_execution_space, nimble_kokkos::kokkos_device_memory_space>;

/*!
 * Contact Manager specific to ArborX library
 *
 * @param interface Contact Interface
 */
ArborXParallelContactManager::ArborXParallelContactManager(
    std::shared_ptr<ContactInterface> interface,
    nimble::DataManager&              data_manager)
    : ParallelContactManager(interface, data_manager)
{
}

void
ArborXParallelContactManager::ComputeParallelContactForce(int step, bool debug_output, nimble::Viewify<2> contact_force)
{
  if (model_data == nullptr) {
    auto model_ptr = this->data_manager_.GetMacroScaleData().get();
    model_data     = dynamic_cast<nimble_kokkos::ModelData*>(model_ptr);
  }

  auto field_ids      = this->data_manager_.GetFieldIDs();
  auto displacement_d = model_data->GetDeviceVectorNodeData(field_ids.displacement);

  auto contact_force_h = model_data->GetHostVectorNodeData(field_ids.contact_force);
  auto contact_force_d = model_data->GetDeviceVectorNodeData(field_ids.contact_force);
  Kokkos::deep_copy(contact_force_d, (double)(0.0));

  this->ApplyDisplacements(displacement_d);

  //--- Set vector to store force
  this->startTimer("Contact:ResetData");
  ContactManager::ZeroContactForce();

  //--- Reset the contact_status flags
  for (size_t jj = 0; jj < contact_faces_d_.extent(0); ++jj) contact_faces_d_(jj).ResetContactData();

  for (size_t jj = 0; jj < contact_nodes_d_.extent(0); ++jj) contact_nodes_d_(jj).ResetContactData();
  this->stopTimer("Contact:ResetData");

  //--- Constraint per ContactManager::ComputeContactForce
  if (penalty_parameter_ <= 0.0) {
    throw std::logic_error("\nError in ComputeContactForce(), invalid penalty_parameter.\n");
  }

  Kokkos::View<details::OutputData*, kokkos_device> results("results", 0);
  Kokkos::View<int*, kokkos_device>                 offset("offset", 0);
  std::set<details::PairData>                       list_;
  auto                                              comm = MPI_COMM_WORLD;

  this->startTimer("ArborX::Search::Def");
  ArborX::DistributedTree<memory_space> dtree(comm, kokkos_device::execution_space{}, contact_faces_d_);
  this->stopTimer("ArborX::Search::Def");

  this->startTimer("ArborX::Search::Query");
  dtree.query(
      kokkos_device::execution_space{},
      details::PredicateTypeNodesRank{contact_nodes_d_, m_rank},
      details::ContactCallback{m_rank, contact_faces_d_, penalty_parameter_, list_},
      results,
      offset);
  this->stopTimer("ArborX::Search::Query");

  this->startTimer("Contact::EnforceInteraction");
  nimble_kokkos::DeviceContactEntityArrayView contact_nodes = contact_nodes_d_;
  nimble_kokkos::DeviceScalarNodeView         force         = force_d_;
  auto                                        numNodes      = contact_nodes_d_.extent(0);
  Kokkos::parallel_for(
      "Update Node Force", numNodes, KOKKOS_LAMBDA(const int i_node) {
        auto& myNode = contact_nodes(i_node);
        for (int j = offset(i_node); j < offset(i_node + 1); ++j) {
          auto tmpOutput = results(j);
          if (!tmpOutput.has_force_) continue;
          //
          myNode.set_contact_status(true);
          myNode.force_1_x_ = tmpOutput.force_node_[0];
          myNode.force_1_y_ = tmpOutput.force_node_[1];
          myNode.force_1_z_ = tmpOutput.force_node_[2];
          //
          myNode.ScatterForceToContactManagerForceVector(force);
        }
      });
  //
  for (size_t iface = 0; iface < contact_faces_d_.extent(0); ++iface) {
    auto& myFace = contact_faces_d_(iface);
    if (myFace.contact_status()) myFace.ScatterForceToContactManagerForceVector(force_d_);
  }
  this->stopTimer("Contact::EnforceInteraction");

  this->GetForces(contact_force_d);
  Kokkos::deep_copy(contact_force_h, contact_force_d);

  // Perform a reduction to obtain correct values on MPI boundaries
  constexpr int vector_dim           = 3;
  auto          myVectorCommunicator = this->data_manager_.GetVectorCommunicator();
  myVectorCommunicator->VectorReduction(vector_dim, contact_force_h);
}

}  // namespace nimble

#endif
