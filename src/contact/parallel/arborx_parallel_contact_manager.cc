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

#include <random>
#include <vector>
#include <iostream>

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
          // TODO ACCESS TRAIT, TO RETURN POINT OR BOX
          nimble::ContactEntity &e = v(i);
          ArborX::Point point1(e.bounding_box_x_min_, e.bounding_box_y_min_, e.bounding_box_z_min_);
          ArborX::Point point2(e.bounding_box_x_max_, e.bounding_box_y_max_, e.bounding_box_z_max_);
          ArborX::Box box(point1, point2);

          return box;
        }
        using memory_space = nimble_kokkos::kokkos_device_memory_space;
    };

    template <>
    struct AccessTraits<nimble_kokkos::DeviceContactEntityArrayView, PredicatesTag>
    {
        static std::size_t size(nimble_kokkos::DeviceContactEntityArrayView const &v) { return v.size(); }

        KOKKOS_FUNCTION static auto get(nimble_kokkos::DeviceContactEntityArrayView const &v, std::size_t i)
        {
          nimble::ContactEntity &e = v(i);
          ArborX::Point point1(e.bounding_box_x_min_, e.bounding_box_y_min_, e.bounding_box_z_min_);
          ArborX::Point point2(e.bounding_box_x_max_, e.bounding_box_y_max_, e.bounding_box_z_max_);
          ArborX::Box box(point1, point2);


          //
          // What does Intersects returns, how is it used afterwards?
          // The intent with the "unspecified" return type in the doc
          // (https://github.com/arborx/ArborX/wiki/ArborX%3A%3Aintersects)
          // is to consider the return type as an implementation detail.
          //
          // If needed, `decltype(ArborX::intersects(ArborX::Box{}))` spells out the type.
          //
          return intersects(box);
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
  {
    Kokkos::View<int *, nimble_kokkos::kokkos_device> indices("indices", 0);
    Kokkos::View<int *, nimble_kokkos::kokkos_device> offset("offset", 0);
    Kokkos::View<int *, nimble_kokkos::kokkos_device> ranks("ranks", 0);
    updateCollisionData(indices, offset, ranks);
  }

  void ArborXParallelContactManager::updateCollisionData(
    Kokkos::View<int *, nimble_kokkos::kokkos_device> &indices,
    Kokkos::View<int *, nimble_kokkos::kokkos_device> &offset,
    Kokkos::View<int *, nimble_kokkos::kokkos_device> &ranks
  ) {

    auto comm = MPI_COMM_WORLD;
    ArborX::DistributedSearchTree<memory_space> dtree(comm, kokkos_device::execution_space{}, contact_faces_d_);

    dtree.query(kokkos_device::execution_space{}, contact_nodes_d_,
                indices, offset, ranks);

  }

  void ArborXParallelContactManager::ComputeParallelContactForce(int step, bool debug_output) {

    //--- Constraint per ContactManager::ComputeContactForce
    if (penalty_parameter_ <= 0.0) {
        throw std::logic_error("\nError in ComputeContactForce(), invalid penalty_parameter.\n");
    }

    Kokkos::View<int *, kokkos_device> indices("indices", 0);
    Kokkos::View<int *, kokkos_device> offset("offset", 0);
    Kokkos::View<int *, kokkos_device> ranks("ranks", 0);

    //--- Update the geometric collision information
    updateCollisionData(indices, offset, ranks);

    //--- Set vector to store force
    ContactManager::zeroContactForce();
    contact_interface->SetContactForce(force_d_);

    // Reset the contact_status flags
    for (size_t jj = 0; jj < contact_faces_d_.extent(0); ++jj)
      contact_faces_d_(jj).set_contact_status(false);

    for (size_t jj = 0; jj < contact_nodes_d_.extent(0); ++jj)
      contact_nodes_d_(jj).set_contact_status(false);

    //
    // The next loop does not track which node is in contact with which face
    // In theory, we could simply loop on the entries in indices.
    //
    double gap = 0.0;
    double normal[3] = {0., 0., 0.};
    ContactManager::PROJECTION_TYPE flag = UNKNOWN;
    ContactEntity::vertex projected;
    for (size_t inode = 0; inode < contact_nodes_d_.extent(0); ++inode) {
      auto &myNode = contact_nodes_d_(inode);
      for (int j = offset(inode); j < offset(inode+1); ++j) {
        if (ranks(j) != m_rank)
          continue;
        //
        //--- Determine whether the node is projected inside the triangular face
        //
        auto &myFace = contact_faces_d_(indices(j));
        ContactManager::SimpleClosestPointProjectionSingle(myNode, myFace,
            &flag, &projected, gap, &normal[0]);
        if ((flag != UNKNOWN) && (gap < 0.0)) {
          contact_faces_d_(indices(j)).set_contact_status(true);
          contact_nodes_d_(inode).set_contact_status(true);
          contact_interface->EnforceNodeFaceInteraction(myNode, myFace, 3, gap, normal, projected.coords_);
        }
      }
    }

  }

}

#endif
