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

// EXAMPLE FROM https://github.com/arborx/ArborX/blob/eddb1d2ceacd8d4bd7bd313c9288ccc6c0840c0d/examples/access_traits/example_host_access_traits.cpp#L48
//      std::vector<ArborX::Point> points;
//
//      // Fill vector with random points in [-1, 1]^3
//      std::uniform_real_distribution<float> dis{-1., 1.};
//      std::default_random_engine gen;
//      auto rd = [&]() { return dis(gen); };
//      std::generate_n(std::back_inserter(points), 100, [&]() {
//          return ArborX::Point{rd(), rd(), rd()};
//      });
//
//      // Pass directly the vector of points to use the access traits defined above
//      ArborX::BVH<Kokkos::HostSpace> bvh{Kokkos::DefaultHostExecutionSpace{},
//                                         points};
//
//      // As a supported alternative, wrap the vector in an unmanaged View
//      bvh = ArborX::BVH<Kokkos::HostSpace>{
//              Kokkos::DefaultHostExecutionSpace{},
//              Kokkos::View<ArborX::Point *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>{
//                      points.data(), points.size()}};
// /EXAMPLE

      // Setting-up the query. Example from https://github.com/arborx/ArborX/blob/eddb1d2ceacd8d4bd7bd313c9288ccc6c0840c0d/examples/access_traits/example_cuda_access_traits.cpp

      // ArborX bounded volume hierarchy setup
      // JLP:  /!\ Does contact_nodes_d_ contains all the information required for our contact search? contact_faces_d_ ?
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

      // Define a copy of contact_nodes_d_ to View in ArborX
      // Number of queries, n = size of contact_nodes_d
      // Size of offset = n + 1
      bvh.query(nimble_kokkos::kokkos_device_execution_space{}, contact_nodes_d_,
                indices, offset);

      ///--- For debugging purposes ---
      for (int i = 0; i < offset.extent(0); ++i)
        std::cout << " offset[" << i << "] = " << offset(i) << "\n";

      for (int i = 0; i < indices.extent(0); ++i)
        std::cout << " indices[" << i << "] = " << indices(i) << "\n";

      for (int i = 0; i < offset.extent(0) - 1; ++i) {
        for (int j = offset(i); j < offset(i + 1); ++j) {
          std::cout << " i " << i << " offset " << offset(i) << " x " << offset(i+1)
                    << " indices " << indices(j) 
                    << std::endl;
        }
      }
      //--------------------------------

  }

  void ArborXSerialContactManager::ComputeSerialContactForce(int step, bool debug_output) {

    // Update collision objects, this will build the trees
//    for ( auto &&node : contact_nodes_ )
//      node.RecomputeKdop();
//    for ( auto &&face : contact_faces_ )
//      face.RecomputeKdop();

/*
    //
    // Extracted from BVH version
    //
    m_nodes->set_entity_data(contact_nodes_);
    m_faces->set_entity_data(contact_faces_);

    m_nodes->broadphase(*m_faces);

    m_last_results.clear();
    m_nodes->for_each_result< NarrowphaseResult >( [this]( const NarrowphaseResult &_res ){
      m_last_results.emplace_back( _res );
    } );
*/

  }

}

#endif
