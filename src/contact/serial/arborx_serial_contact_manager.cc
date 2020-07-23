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

#ifdef NIMBLE_HAVE_ARBORX
#include "arborx_serial_contact_manager.h"

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
  using arborx_bvh = ArborX::BVH<memory_space>;

  /*!
   * Contact Manager specific to ArborX library
   *
   * Global variables contact_faces_d_
   *
   * @param interface Contact Interface
   */
  ArborXSerialContactManager::ArborXSerialContactManager(std::shared_ptr<ContactInterface> interface)
        : SerialContactManager(interface)
  {
    /// TODO Asj NM & RJ whether this is needed at constructor time
    updateCollisionData();
  }

  void ArborXSerialContactManager::updateCollisionData() {

    //
    // primitives are faces
    //
    arborx_bvh bvh{nimble_kokkos::kokkos_device_execution_space{},
                   contact_faces_d_};

    //
    // indices : position of the primitives that satisfy the predicates.
    // offsets : predicate offsets in indices.
    //
    // indices stores the indices of the objects that satisfy the predicates.
    // offset stores the locations in the indices view that start a predicate, that is,
    // predicates(q) is satisfied by indices(o) for primitives(q) <= o < primitives(q+1).
    //
    // Following the usual convention, offset(n) = nnz, where n is the number of queries
    // that were performed and nnz is the total number of collisions.
    //
    // (From https://github.com/arborx/ArborX/wiki/ArborX%3A%3ABoundingVolumeHierarchy%3A%3Aquery )
    //
    Kokkos::View<int *, nimble_kokkos::kokkos_device> indices("indices", 0);
    Kokkos::View<int *, nimble_kokkos::kokkos_device> offset("offset", 0);

    // Number of queries, n = size of contact_nodes_d
    // Size of offset = n + 1
    bvh.query(nimble_kokkos::kokkos_device_execution_space{}, contact_nodes_d_,
              indices, offset);

    // Reset the contact_status flags
    for (size_t jj = 0; jj < contact_faces_d_.extent(0); ++jj)
      contact_faces_d_[jj].set_contact_status(false);

    //
    // The next loop does not track which node is in contact with which face
    // In theory, we could simply loop on the entries in indices.
    //
    for (size_t inode = 0; inode < contact_nodes_d_.extent(0); ++inode) {
      for (int j = offset(inode); j < offset(inode+1); ++j) {
        contact_faces_d_[indices(j)].set_contact_status(true);
      }
    }

    ///--- For debugging purposes ---
    std::cout << " offset size " << offset.extent(0)
              << " nnz " << offset(contact_nodes_d_.extent(0))
              << "\n";

    std::cout << " indices size " << indices.extent(0) << "\n";
    //--------------------------------
  }

  void ArborXSerialContactManager::ComputeSerialContactForce(int step, bool debug_output) {
    //--- Update the geometric collision information
    updateCollisionData();
    //---
  }

}

#endif
