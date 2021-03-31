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

#include <algorithm>
#include <array>
#include <cstring>
#include <exception>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
#pragma once

namespace nimble {
namespace serialization {
template <class T>
struct meta
{
  static constexpr bool mem = std::is_pod<T>::value;
  // static constexpr int _size = mem ? sizeof(T) : 0;
};
template <class T>
struct meta<T*>
{
  static constexpr bool mem = meta<T>::mem;
  // static constexpr size_t _size = meta<T>::_size;
};
template <class T, size_t arr_size>
struct meta<std::array<T, arr_size>>
{
  static constexpr bool mem = meta<T>::mem;
  // static constexpr size_t _size = mem ? meta<T>::_size * arr_size : 0;
};
template <class A, class B>
struct meta<std::pair<A, B>>
{
  static constexpr bool mem = meta<A>::mem && meta<B>::mem;
  // static constexpr size_t _size = mem ? meta<A>::_size + meta<B>::_size : 0;
};
template <class T, bool constantsize, bool... constant_size>
struct serial;
template <class T>
using serial_t = serial<typename std::remove_cv<T>::type, meta<T>::mem>;
template <class T>
struct serial<T, false>
{
  // No size function declared
  // Not even the gods can help you now
};
template <class T>
struct serial<T, true>
{
  constexpr static size_t
  size(const T& obj)
  {
    return sizeof(T);
  }
  static void
  pack(const T& src, char*& dest)
  {
    std::memcpy(dest, &src, sizeof(T));
    dest += sizeof(T);
  }
  static void
  unpack(T& dest, char*& src)
  {
    std::memcpy(&dest, src, sizeof(T));
    src += sizeof(T);
  }
};
template <class T>
struct serial<T*, false>
{
  static size_t
  size(const T* obj)
  {
    return serial_t<T>::size(*obj);
  }
  static void
  pack(const T* src, char*& dest)
  {
    serial_t<T>::pack(*src, dest);
  }
  static void
  unpack(T* dest, char*& src)
  {
    serial_t<T>::unpack(*dest, src);
  }
};
template <class T>
struct serial<T*, true>
{
  static size_t
  size(const T* obj)
  {
    return serial_t<T>::size(*obj);
  }
  static void
  pack(const T* src, char*& dest)
  {
    serial_t<T>::pack(*src, dest);
  }
  static void
  unpack(T* dest, char*& src)
  {
    serial_t<T>::unpack(*dest, src);
  }
};
template <class A, class B>
struct serial<std::pair<A, B>, false>
{
  static size_t
  size(const std::pair<A, B>& obj)
  {
    return serial_t<A>::size(obj.first) + serial_t<B>::size(obj.second);
  }
  static void
  pack(const std::pair<A, B>& src, char*& dest)
  {
    serial_t<A>::pack(src.first, dest);
    serial_t<B>::pack(src.second, dest);
  }
  static void
  unpack(std::pair<A, B>& dest, char*& src)
  {
    serial_t<A>::unpack(dest.first, src);
    serial_t<B>::unpack(dest.second, src);
  }
};
template <class A, class B>
struct serial<std::pair<A, B>, true>
{
  constexpr static size_t
  size(const std::pair<A, B>& obj)
  {
    return serial_t<A>::size(obj.first) + serial_t<B>::size(obj.second);
  }
  static void
  pack(const std::pair<A, B>& obj, char*& dest)
  {
    serial_t<A>::pack(obj.first, dest);
    serial_t<B>::pack(obj.second, dest);
  }
  static void
  unpack(std::pair<A, B>& dest, char*& src)
  {
    serial_t<A>::unpack(dest.first, src);
    serial_t<B>::unpack(dest.second, src);
  }
};
template <class T>
struct serial<std::vector<T>, false, false>
{
  static size_t
  size(const std::vector<T>& vect)
  {
    size_t total_size = 0;
    for (auto& elem : vect) total_size += serial<T, false>::size(elem);
    return total_size + sizeof(size_t);
  }
  static void
  pack(const std::vector<T>& vect, char*& dest)
  {
    serial_t<size_t>::pack(vect.size(), dest);
    for (auto& elem : vect) serial_t<T>::pack(elem, dest);
  }
  static void
  unpack(std::vector<T>& dest, char*& src)
  {
    size_t size;
    serial_t<size_t>::unpack(size, src);
    dest.resize(size);
    for (T& elem : dest) serial<T, false>::unpack(elem, src);
  }
};
template <class T>
struct serial<std::vector<T>, false, true>
{
  static size_t
  size(const std::vector<T>& vect)
  {
    return serial_t<T>::size(vect[0]) * vect.size() + sizeof(size_t);
  }
  static void
  pack(const std::vector<T>& src, char*& dest)
  {
    serial_t<size_t>::pack(src.size(), dest);
    for (auto& elem : src) serial_t<T>::pack(elem, dest);
  }
  static void
  unpack(std::vector<T>& dest, char*& src)
  {
    size_t size;
    serial_t<size_t>::unpack(size, src);
    dest.resize(size);
    for (T& elem : dest) serial_t<T>::unpack(elem, src);
  }
};
template <class T>
struct serial<std::vector<T>, false>
{
  static size_t
  size(const std::vector<T>& obj)
  {
    return serial<std::vector<T>, false, meta<T>::mem>::size(obj);
  }
  static void
  pack(const std::vector<T>& src, char*& dest)
  {
    serial<std::vector<T>, false, meta<T>::mem>::pack(src, dest);
  }
  static void
  unpack(std::vector<T>& dest, char*& src)
  {
    serial<std::vector<T>, false, meta<T>::mem>::unpack(dest, src);
  }
};
template <class T, size_t arr_size>
struct serial<std::array<T*, arr_size>, true>
{
  constexpr static size_t
  size(const std::array<T*, arr_size>& arr)
  {
    return sizeof(T) * arr_size;
  }
  static void
  pack(const std::array<T*, arr_size>& src, char*& dest)
  {
    for (auto* elem : src) serial<T, true>::pack(*elem, dest);
  }
  static void
  unpack(std::array<T*, arr_size>& dest, char*& src)
  {
    for (auto* elem : dest) serial<T, true>::unpack(*elem, src);
  }
};
template <class T, size_t arr_size>
struct serial<std::array<T, arr_size>, false>
{
  static size_t
  size(const std::array<T, arr_size>& arr)
  {
    size_t total_size = 0;
    for (auto& elem : arr) total_size += serial<T, false>::size(elem);
    return total_size;
  }
  static void
  pack(const std::array<T, arr_size>& src, char*& dest)
  {
    for (auto& elem : src) serial<T, false>::pack(elem, dest);
  }
  static void
  unpack(std::array<T, arr_size>& dest, char*& src)
  {
    for (T& elem : dest) serial<T, false>::unpack(elem, src);
  }
};
template <class T, class dataT>
size_t
pack(const T& object, std::vector<dataT>& vect)
{
  static_assert(std::is_pod<dataT>::value, "vector must be of POD type");
  size_t size = serial_t<T>::size(object);
  vect.resize(size / sizeof(dataT) + 1);
  char* dataptr = (char*)&vect[0];
  serial_t<T>::pack(object, dataptr);
  return size;
}

template <class T, class dataT>
size_t
pack_avoid_resize(const T& object, std::vector<dataT>& vect)
{
  static_assert(std::is_pod<dataT>::value, "vector must be of POD type");
  size_t size = serial_t<T>::size(object);
  if (size > vect.size()) vect.resize(size / sizeof(dataT) + 1);
  char* dataptr = (char*)&vect[0];
  serial_t<T>::pack(object, dataptr);
  return size;
}

template <class T, class dataT>
void
unpack(T& object, const std::vector<dataT>& vect)
{
  static_assert(std::is_pod<dataT>::value, "vector must be of POD type");
  char* dataptr = (char*)&vect[0];
  serial_t<T>::unpack(object, dataptr);
}
}  // namespace serialization
}  // namespace nimble
