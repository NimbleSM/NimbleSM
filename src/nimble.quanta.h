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

#ifndef NIMBLE_QUANTA_H
#define NIMBLE_QUANTA_H

namespace nimble {
namespace quanta {
template <class T>
T
declval_exact() noexcept;

template <class T>
typename std::remove_reference<T>::type&
declref() noexcept;

// Gets the literal type of a container when it's iterator is dereferenced
// for example, 'iterated_t<std::vector<int>>' is of type 'int&', and
// 'iterated_t<const std::vector<int>>' is of type 'const int&'
template <class list_t>
using iterated_t = decltype(*std::begin(declref<list_t>()));

// Gets the unmodified type of a container when it's iterator is defererenced
// For example, 'elem_t<std::vector<int>>' is of type 'int'
template <class list_t>
using elem_t = typename std::decay<iterated_t<list_t>>::type;

template <class F, class... Args>
using return_t = decltype(declref<F>()(declval_exact<Args>()...));

template <class F, class... list_t>
using transformed_iterated_t = return_t<F, iterated_t<list_t>...>;

template <class T>
auto
len(T&& obj) -> decltype(obj.size())
{
  return obj.size();
}

template <class list_t>
struct indexer_t
{
  list_t& list;
  constexpr indexer_t(list_t& list) : list(list) {}
  constexpr indexer_t(const indexer_t&) = default;
  template <class int_t>
  auto
  operator()(int_t&& index) const -> decltype(list[index])
  {
    return list[index];
  }
};

template <class T>
struct indexer_t<T*>
{
  T* const ptr;
  constexpr indexer_t(T* const ptr) : ptr(ptr) {}
  indexer_t(const indexer_t&) = default;
  template <class int_t>
  T&
  operator()(int_t&& index) const
  {
    return ptr[index];
  }
};

template <class T>
indexer_t<T>
make_indexer(T& list)
{
  return indexer_t<T>(list);
}

template <class T>
indexer_t<T*>
make_indexer(T* ptr)
{
  return indexer_t<T*>(ptr);
}

template <class int_t = int>
class invoke_counter_t
{
  mutable int_t count;

 public:
  invoke_counter_t() : count{} {}
  invoke_counter_t(int_t i) : count{i} {}
  invoke_counter_t(const invoke_counter_t&) = default;
  invoke_counter_t(invoke_counter_t&&)      = default;
  invoke_counter_t&
  operator=(const invoke_counter_t&) = default;
  invoke_counter_t&
  operator=(invoke_counter_t&&) = default;
  void
  reset(const int_t& value = 0)
  {
    count = value;
  }
  int_t&
  get_count() const
  {
    return count;
  }
  auto
  operator()() const -> decltype(count++)
  {
    return count++;
  }
  auto
  increment() const -> decltype(count++)
  {
    return count++;
  }
};
template <class int_t>
invoke_counter_t<int_t>
make_counter(int_t initial = 0)
{
  return invoke_counter_t<int_t>{initial};
}

template <class list_t, class F>
void
remap(list_t& list, F&& func)
{
  for (auto& elem : list) { elem = func(elem); }
}

}  // namespace quanta
}  // namespace nimble

#endif  // NIMBLE_QUANTA_H
