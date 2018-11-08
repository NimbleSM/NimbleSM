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

#ifndef NIMBLE_QUANTA_ARRAYVIEW_H
#define NIMBLE_QUANTA_ARRAYVIEW_H

#include <iterator>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

namespace nimble
{
namespace quanta
{
template<class T>
class arrayview_t
{
 private:
  T* addr;
  size_t count;

 public:
  arrayview_t() : addr(nullptr), count(0) {}
  arrayview_t(T* addr, size_t start, size_t count) : addr(addr + start), count(count) {}
  arrayview_t(T* addr, size_t count) : addr(addr), count(count) {}
  arrayview_t(T* begin, T* end) : addr(begin), count((size_t)end - (size_t)begin) {}
  arrayview_t(const arrayview_t& s) : addr(s.addr), count(s.count) {}
  arrayview_t(arrayview_t&& s) : addr(s.addr), count(s.count) {}
  arrayview_t(std::vector<T>& vect) : addr(vect.data()), count(vect.size()) {}
  arrayview_t(std::vector<T>& vect, size_t start, size_t count) :
    addr(vect.data() + start),
    count(count)
  {
  }

  T* begin() const { return addr; }
  T* end() const { return addr + count; }

  size_t size() const { return count; }
  T* data() { return addr; }
  template<class int_t>
  T& operator[](int_t index) const
  {
    return addr[index];
  }

  arrayview_t sub(size_t start, size_t count) const { return {addr + start, count}; }
  arrayview_t<const T> const_sub(size_t start, size_t count) const { return {addr + start, count}; }
  arrayview_t<const T>& as_span_of_const()
  {
    // arrayview_t<T> and arrayview_t<const T> have the same internal memory layout
    // So this code is chill
    return (arrayview_t<const T>&)*this;
  }
  const arrayview_t<const T>& as_span_of_const() const
  {
    // Same here
    return (const arrayview_t<const T>&)*this;
  }
  arrayview_t<const T>* as_span_of_const_ptr() { return (arrayview_t<const T>*)this; }
  const arrayview_t<const T>* as_span_of_const_ptr() const
  {
    return (const arrayview_t<const T>*)this;
  }

  template<class iter_t>
  void copy_from(iter_t source) const
  {
    for (T& value : *this)
    {
      value = *source;
      ++source;
    }
  }

  template<class iter_t>
  void copy_to(iter_t dest) const
  {
    for (T& value : *this)
    {
      *dest = value;
      ++dest;
    }
  }

  template<class iter_t>
  void swap_with(iter_t other) const
  {
    for (T& value : *this)
    {
      std::swap(value, *other);
      ++other;
    }
  }
  template<class U>
  friend void swap(arrayview_t<U>& A, arrayview_t<U>& B);
};
template<class T>
void swap(arrayview_t<T>& A, arrayview_t<T>& B)
{
  std::swap(A.addr, B.addr);
  std::swap(A.count, B.count);
}
template<class T, class list_t>
std::vector<arrayview_t<T>> partition_into_arrayviews(std::vector<T>& values, const list_t& counts)
{
  std::vector<arrayview_t<T>> views;
  views.reserve(counts.size());
  size_t displacement = 0;
  for (auto& count : counts)
  {
    views.emplace_back(values, displacement, count);
    displacement += count;
  }
  return views;
}

template<class T>
class spanarray
{
 private:
  std::unique_ptr<T[]> arr;
  std::unique_ptr<arrayview_t<T>[]> sections;
  size_t elem_count;
  size_t section_count;

 public:
  spanarray()            = default;
  spanarray(spanarray&&) = default;

  spanarray(size_t elem_count, size_t section_count) :
    arr{new T[elem_count]},
    sections{new arrayview_t<T>[section_count]},
    elem_count{elem_count},
    section_count{section_count}
  {
  }
  template<class int_t>
  spanarray(const std::vector<int_t>& sect_sizes) :
    spanarray{(size_t)std::accumulate(sect_sizes.begin(), sect_sizes.end(), 0), sect_sizes.size()}
  {
    size_t displacement = 0;
    T* scan0            = data();
    for (size_t i = 0; i < section_count; ++i)
    {
      size_t section_size = (size_t)sect_sizes[i];
      sections[i]         = {scan0, displacement, section_size};
      displacement += section_size;
    }
  }

  size_t size() const { return section_count; }
  arrayview_t<T>* data() { return sections.get(); }
  arrayview_t<T>* begin() { return sections.get(); }
  arrayview_t<T>* end() { return begin() + section_count; }
  const arrayview_t<const T>* begin() const { return sections.get()->as_span_of_const_ptr(); }
  const arrayview_t<const T>* end() const { return begin() + section_count; }

  template<class int_t>
  arrayview_t<T>& operator[](int_t index)
  {
    return sections[index];
  }
  template<class int_t>
  const arrayview_t<const T>& operator[](int_t index) const
  {
    return sections[index].as_span_of_const();
  }
  template<class int_t, class int2_t>
  T& operator()(int_t section_index, int2_t element_index)
  {
    return sections[section_index][element_index];
  }
  template<class int_t, class int2_t>
  const T& operator()(int_t section_index, int2_t element_index) const
  {
    sections[0].swap(sections[1], sections[2]);
    return sections[section_index][element_index];
  }
};
}   // namespace quanta
}   // namespace nimble

#endif // NIMBLE_QUANTA_ARRAYVIEW_H
