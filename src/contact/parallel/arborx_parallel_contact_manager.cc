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
#include "nimble_kokkos_defs.h"

#include <ArborX.hpp>
#include <Kokkos_Core.hpp>

#include <iostream>
#include <random>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace nimble {

namespace details {

constexpr int dim = 3;

struct PredicateTypeNodesRank
{
    nimble_kokkos::DeviceContactEntityArrayView nodes_;
    int rank_;
};

struct OutputData
{
  int index_;
  int rank_;
  double coord_[dim];
  bool has_force_;
  double force_node_[dim];
};

struct PairData
{
  int pred_rank_ = 0;
  int pred_index_ = -1;
  int prim_rank_ = 0;
  int prim_index_ = -1;
  //
  PairData(int r_, int i_, int f_r_, int f_i_)
      : pred_rank_(r_), pred_index_(i_), prim_rank_(f_r_), prim_index_(f_i_)
  {}
  //
  bool operator<(const PairData &rhs) const {
    return (std::tie(pred_rank_, pred_index_, prim_rank_, prim_index_)
            < std::tie(rhs.pred_rank_, rhs.pred_index_, rhs.prim_rank_, rhs.prim_index_));
  }
};

struct ContactCallback
{
  const int rank_;
  nimble_kokkos::DeviceContactEntityArrayView &faces_;
  const double penalty_;
  std::set< PairData > &list_;
  //
  template <typename Predicate, typename OutputFunctor>
  KOKKOS_FUNCTION void operator()(Predicate const &pred, int f_primitive,
                                  OutputFunctor const &out) const
  {
    auto const& p_data = getData(pred);  // <- type OutputData
    auto const& p_geometry = getGeometry(pred); // <- type Box
    //
    //--- Define copy contact entity
    //--- For the node, we only store the coordinates.
    //
    ContactEntity myNode(ContactEntity::ContactEntityType::NODE, 0, p_data.coord_, 0.0, 0);
    ContactEntity myFace;
    faces_(f_primitive).ExportGeometryInto(myFace);
    //
    double gap = 0.0;
    double normal[dim] = {0., 0., 0.};
    bool inside = false;
    double facet_coordinates[dim] = {0., 0., 0.};
    //
    ContactManager::Projection(myNode, myFace, inside, gap, &normal[0],
                               &facet_coordinates[0]);
    //
    double force[dim] = {0., 0., 0.};
    if (inside) {
      ///
      /// TODO Check the factor 3 with RJ and NM
      ///
      details::getContactForce(penalty_ / static_cast<double>(3), gap, normal, force);
      //
      auto result = list_.insert(PairData{p_data.rank_, p_data.index_, rank_, f_primitive});
      if (result.second) { // if the insertion took place
        myFace.SetNodalContactForces(force, &facet_coordinates[0]);
        faces_(f_primitive).set_contact_status(true);
        faces_(f_primitive).force_1_x_ += myFace.force_1_x_;
        faces_(f_primitive).force_1_y_ += myFace.force_1_y_;
        faces_(f_primitive).force_1_z_ += myFace.force_1_z_;
        faces_(f_primitive).force_2_x_ += myFace.force_2_x_;
        faces_(f_primitive).force_2_y_ += myFace.force_2_y_;
        faces_(f_primitive).force_2_z_ += myFace.force_2_z_;
        faces_(f_primitive).force_3_x_ += myFace.force_3_x_;
        faces_(f_primitive).force_3_y_ += myFace.force_3_y_;
        faces_(f_primitive).force_3_z_ += myFace.force_3_z_;
      }
    }
    //
    out({f_primitive, rank_, {p_data.coord_[0], p_data.coord_[1], p_data.coord_[2]},
         inside, {-force[0], -force[1], -force[2]}});
  }
};

} }

namespace ArborX
{
    template <>
    struct AccessTraits<nimble_kokkos::DeviceContactEntityArrayView, PrimitivesTag>
    {
        // size returns the number of elements in the View
        static std::size_t size(nimble_kokkos::DeviceContactEntityArrayView const &v) { return v.size(); }

        /// Returns an ArborX::Box for each contact entity within the nimble view
        ///
        /// \param v
        /// \param i
        /// \return ArborX::Box
        KOKKOS_FUNCTION static ArborX::Box get(nimble_kokkos::DeviceContactEntityArrayView const &v, std::size_t i)
        {
          nimble::ContactEntity &e = v(i);
          ArborX::Point point1(e.bounding_box_x_min_, e.bounding_box_y_min_, e.bounding_box_z_min_);
          ArborX::Point point2(e.bounding_box_x_max_, e.bounding_box_y_max_, e.bounding_box_z_max_);
          ArborX::Box box(point1, point2);
          return box;
        }
        using memory_space = nimble_kokkos::kokkos_device_memory_space;
    };

    template <>
    struct AccessTraits<nimble::details::PredicateTypeNodesRank, PredicatesTag>
    {
        static std::size_t size(nimble::details::PredicateTypeNodesRank const &v) { return v.nodes_.extent(0); }

        KOKKOS_FUNCTION static auto get(nimble::details::PredicateTypeNodesRank const &v, std::size_t i)
        {
          nimble::ContactEntity &e = v.nodes_(i);
          ArborX::Point point1(e.bounding_box_x_min_, e.bounding_box_y_min_, e.bounding_box_z_min_);
          ArborX::Point point2(e.bounding_box_x_max_, e.bounding_box_y_max_, e.bounding_box_z_max_);
          ArborX::Box box(point1, point2);
          //
          return attach(intersects(box),
                        nimble::details::OutputData{static_cast<int>(i), v.rank_,
                                                    {e.coord_1_x_, e.coord_1_y_, e.coord_1_z_},
                                                    false, {0, 0, 0}});
        }
        using memory_space = nimble_kokkos::kokkos_device_memory_space;
    };
} // namespace ArborX

namespace nimble {

  using memory_space = nimble_kokkos::kokkos_device_memory_space;
  using kokkos_device = Kokkos::Device< nimble_kokkos::kokkos_device_execution_space, nimble_kokkos::kokkos_device_memory_space>;

  /*!
   * Contact Manager specific to ArborX library
   *
   * @param interface Contact Interface
   */
  ArborXParallelContactManager::ArborXParallelContactManager(std::shared_ptr<ContactInterface> interface)
        : ParallelContactManager(interface)
  { }

  void ArborXParallelContactManager::ComputeParallelContactForce(int step, bool debug_output) {

    //--- Set vector to store force
    ContactManager::zeroContactForce();

    //--- Reset the contact_status flags
    for (size_t jj = 0; jj < contact_faces_d_.extent(0); ++jj)
      contact_faces_d_(jj).ResetContactData();

    for (size_t jj = 0; jj < contact_nodes_d_.extent(0); ++jj)
      contact_nodes_d_(jj).ResetContactData();

    //--- Constraint per ContactManager::ComputeContactForce
    if (penalty_parameter_ <= 0.0) {
      throw std::logic_error("\nError in ComputeContactForce(), invalid penalty_parameter.\n");
    }

    Kokkos::View<details::OutputData *, kokkos_device> results("results", 0);
    Kokkos::View<int *, kokkos_device> offset("offset", 0);
    std::set< details::PairData > list_;

    auto comm = MPI_COMM_WORLD;
    ArborX::DistributedSearchTree<memory_space> dtree(comm,
                                                      kokkos_device::execution_space{},
                                                      contact_faces_d_);

    dtree.query(kokkos_device::execution_space{},
                details::PredicateTypeNodesRank{contact_nodes_d_, m_rank},
                details::ContactCallback{m_rank, contact_faces_d_, penalty_parameter_, list_},
                results, offset);

    for (size_t inode = 0; inode < contact_nodes_d_.extent(0); ++inode) {
      auto &myNode = contact_nodes_d_(inode);
      for (int j = offset(inode); j < offset(inode+1); ++j) {
        auto tmpOutput = results(j);
        if (!tmpOutput.has_force_)
          continue;
        //
        contact_nodes_d_(inode).set_contact_status(true);
        myNode.force_1_x_ = tmpOutput.force_node_[0];
        myNode.force_1_y_ = tmpOutput.force_node_[1];
        myNode.force_1_z_ = tmpOutput.force_node_[2];
        //
        myNode.ScatterForceToContactManagerForceVector(force_d_);
        //
      }
    }

    for (size_t iface = 0; iface < contact_faces_d_.extent(0); ++iface) {
      auto &myFace = contact_faces_d_(iface);
      if (myFace.contact_status() > 0.0)
        myFace.ScatterForceToContactManagerForceVector(force_d_);
    }

  }

}

#endif
