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

#ifndef NIMBLESM_ARBORX_UTILS_H
#define NIMBLESM_ARBORX_UTILS_H

#include "nimble_contact_manager.h"
#include "nimble_defs.h"

#ifdef NIMBLE_HAVE_ARBORX

#include <ArborX.hpp>

//
// Need to specialize AccessTrait< ..., {PrimitivesTag, PredicatesTag} >
// Here the first template parameter is a container of nimble::ContactEntity
//
// See https://github.com/arborx/ArborX/wiki/ArborX%3A%3AAccessTraits
//
namespace ArborX {
template <>
struct AccessTraits<nimble_kokkos::DeviceContactEntityArrayView, PrimitivesTag>
{
  // size returns the number of elements in the View
  static std::size_t
  size(nimble_kokkos::DeviceContactEntityArrayView const& v)
  {
    return v.size();
  }

  /// Returns an ArborX::Box for each contact entity within the nimble view
  ///
  /// \param v
  /// \param i
  /// \return ArborX::Box
  KOKKOS_FUNCTION static ArborX::Box
  get(nimble_kokkos::DeviceContactEntityArrayView const& v, std::size_t i)
  {
    nimble::ContactEntity& e = v(i);
    ArborX::Point          point1(e.bounding_box_x_min_, e.bounding_box_y_min_, e.bounding_box_z_min_);
    ArborX::Point          point2(e.bounding_box_x_max_, e.bounding_box_y_max_, e.bounding_box_z_max_);
    ArborX::Box            box(point1, point2);

    return box;
  }
  using memory_space = nimble_kokkos::kokkos_device_memory_space;
};

template <>
struct AccessTraits<nimble_kokkos::DeviceContactEntityArrayView, PredicatesTag>
{
  static std::size_t
  size(nimble_kokkos::DeviceContactEntityArrayView const& v)
  {
    return v.size();
  }

  KOKKOS_FUNCTION static auto
  get(nimble_kokkos::DeviceContactEntityArrayView const& v, std::size_t i)
  {
    nimble::ContactEntity& e = v(i);
    ArborX::Point          point1(e.bounding_box_x_min_, e.bounding_box_y_min_, e.bounding_box_z_min_);
    ArborX::Point          point2(e.bounding_box_x_max_, e.bounding_box_y_max_, e.bounding_box_z_max_);
    ArborX::Box            box(point1, point2);

    //
    // What does Intersects returns, how is it used afterwards?
    // The intent with the "unspecified" return type in the doc
    // (https://github.com/arborx/ArborX/wiki/ArborX%3A%3Aintersects)
    // is to consider the return type as an implementation detail.
    //
    // If needed, `decltype(ArborX::intersects(ArborX::Box{}))` spells out the type.
    //
    return ArborX::attach(intersects(box), (int)i);
  }
  using memory_space = nimble_kokkos::kokkos_device_memory_space;
};
}  // namespace ArborX

#endif  // #ifdef NIMBLE_HAVE_ARBORX

#endif  // NIMBLESM_ARBORX_UTILS_H
