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
    template <typename T, typename Tag>
    struct AccessTraits<std::vector<T>, Tag>
    {
        static std::size_t size(std::vector<T> const &v) { return v.size(); }
        KOKKOS_FUNCTION static T const &get(std::vector<T> const &v, std::size_t i)
        {
          // TODO ACCESS TRAIT, TO RETURN POINT OR BOX
          return v[i];
        }
        using memory_space = Kokkos::HostSpace;
    };
} // namespace ArborX

namespace nimble {

  using arborx_bvh_type = ArborX::BVH<nimble_kokkos::kokkos_device_memory_space>;



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

      arborx_bvh_type bvh2{nimble_kokkos::kokkos_device_execution_space{},
                           contact_nodes_d_};

  }

  void ArborXSerialContactManager::ComputeSerialContactForce(int step, bool debug_output) {
  }

//  ArborXSerialContactManager::ArborXSerialContactManager(ArborXSerialContactManager &&) noexcept = default;
//  ArborXSerialContactManager & ArborXSerialContactManager::operator=( ArborXSerialContactManager &&) noexcept = default;
//  ArborXSerialContactManager::~ArborXSerialContactManager() = default;

}

#endif