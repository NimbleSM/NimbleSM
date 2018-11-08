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

#include "nimble.quanta.h"

void
nimble::quanta::radixsort256(uint32_t* array, uint32_t* buffer, uint32_t count)
{
  uint32_t counts[4][256]{};   // counts are default-initialized to 0
  uint32_t* scan[256];
  for (decltype(count) i = 0; i < count; ++i)
  {
    auto value = array[i];
    ++counts[0][value & 0xff];
    ++counts[1][(value >> 8) & 0xff];
    ++counts[2][(value >> 16) & 0xff];
    ++counts[3][(value >> 24) & 0xff];
  }

  radixsort_partial<0, 8>(array, buffer, counts[0], scan, count);
  radixsort_partial<8, 8>(buffer, array, counts[1], scan, count);
  radixsort_partial<16, 8>(array, buffer, counts[2], scan, count);
  radixsort_partial<24, 8>(buffer, array, counts[3], scan, count);
}

void
nimble::quanta::radix_order256(uint32_t* array, uint32_t* ordering, uint32_t* buffer, uint32_t count)
{
  uint32_t counts[4][256]{};
  uint32_t* scan[256];
  for (uint32_t i = 0; i < count; ++i)
  {
    ordering[i]    = i;
    uint32_t value = array[i];
    ++counts[0][value & 0xff];
    ++counts[1][(value >> 8) & 0xff];
    ++counts[2][(value >> 16) & 0xff];
    ++counts[3][(value >> 24) & 0xff];
  }

  radix_order_partial<0, 8>(array, ordering, buffer, counts[0], scan, count);
  radix_order_partial<8, 8>(array, buffer, ordering, counts[1], scan, count);
  radix_order_partial<16, 8>(array, ordering, buffer, counts[2], scan, count);
  radix_order_partial<24, 8>(array, buffer, ordering, counts[3], scan, count);
}

void
nimble::quanta::PackIDs(std::vector<int>& id_array)
{
  size_t size = id_array.size();
  std::vector<uint32_t> ordering(size);
  std::vector<uint32_t> buffer(size);
  radix_order256((uint32_t*)id_array.data(), ordering.data(), buffer.data(), size);
  remap_ids_with_ordering(id_array, ordering);
}

