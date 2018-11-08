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

namespace nimble
{
namespace quanta
{
template<class F1, class... F2>
struct composition_t
{
  typedef composition_t<F1> f_t;
  typedef composition_t<F2...> g_t;
  g_t g;
  f_t f;
  composition_t(const composition_t<F1>& f, const composition_t<F2>&... g) : g{g...}, f{f} {}
  composition_t(composition_t<F1>&& f, composition_t<F2>&&... g) :
    g{std::move(g)...},
    f{std::move(f)}
  {
  }
  composition_t(const composition_t& comp) : g(comp.g), f(comp.f) {}
  composition_t(composition_t&& comp) : g(std::move(comp.g)), f(std::move(comp.f)) {}

  template<class... Args>
  auto operator()(Args&&... args) -> decltype(f(g(std::forward<decltype(args)>(args)...)))
  {
    return f(g(std::forward<decltype(args)>(args)...));
  }
  template<class... Args>
  auto operator()(Args&&... args) const -> decltype(f(g(std::forward<decltype(args)>(args)...)))
  {
    return f(g(std::forward<decltype(args)>(args)...));
  }
};
template<class F1>
struct composition_t<F1>
{
  typedef F1 f_t;
  f_t f;
  composition_t(const f_t& f) : f(f) {}
  composition_t(f_t&& f) : f(std::move(f)) {}
  composition_t(const composition_t& comp) : f(comp.f) {}
  composition_t(composition_t&& comp) : f(std::move(comp.f)) {}

  template<class... Args>
  auto operator()(Args&&... args) -> decltype(f(std::forward<decltype(args)>(args)...))
  {
    return f(std::forward<decltype(args)>(args)...);
  }
  template<class... Args>
  auto operator()(Args&&... args) const -> decltype(f(std::forward<decltype(args)>(args)...))
  {
    return f(std::forward<decltype(args)>(args)...);
  }
};
template<class F1>
struct composition_t<F1&>
{
  typedef F1 f_t;
  f_t& f;
  composition_t(f_t& f) : f(f) {}
  composition_t(const composition_t& comp) : f(comp.f) {}

  template<class... Args>
  auto operator()(Args&&... args) -> decltype(f(std::forward<decltype(args)>(args)...))
  {
    return f(std::forward<decltype(args)>(args)...);
  }
  template<class... Args>
  auto operator()(Args&&... args) const -> decltype(f(std::forward<decltype(args)>(args)...))
  {
    return f(std::forward<decltype(args)>(args)...);
  }
};

template<class F1>
auto compose(F1&& f) -> composition_t<F1>
{
  return composition_t<F1>{f};
}
template<class F1, class F2, class... F3>
auto compose(F1&& f1, F2&& f2, F3&&... f3) -> composition_t<F1, F2, F3...>
{
  return composition_t<F1, F2, F3...>{f1, f2, f3...};
}
template<class T>
struct fwd_or_move
{
  typedef T& type;
  typedef T base_type;
  template<class U>
  static U& pass(U& u)
  {
    return u;
  }
};
template<class T>
struct fwd_or_move<T&>
{
  typedef T& type;
  typedef T base_type;
  template<class U>
  static U& pass(U& u)
  {
    return u;
  }
};
template<class T>
struct fwd_or_move<T&&>
{
  typedef T&& type;
  typedef T base_type;
  template<class U>
  static U&& pass(U&& u)
  {
    return std::move(u);
  }
};
template<class F, class List_t>
void apply_discard(List_t&& list, F&& f)
{
  for (auto& elem : list)
  {
    f(fwd_or_move<decltype(list)>::pass(elem));
  }
}
template<class T>
struct assigner_t
{
  T& value;
  assigner_t(T& value) : value(value) {}
  assigner_t(const assigner_t& ass) : value(ass.value) {}
  template<class U>
  T& operator()(U&& thing) const
  {
    value = std::forward<decltype(thing)>(thing);
    return value;
  }
};
struct ostream_printer_t
{
  std::ostream& stream;
  ostream_printer_t(std::ostream& o) : stream(o) {}
  ostream_printer_t(const ostream_printer_t&) = default;
  template<class T>
  std::ostream& operator()(T& value) const
  {
    return (stream << value);
  }
};
template<class T>
assigner_t<T> make_assigner(T& value)
{
  return assigner_t<T>(value);
}
template<class list_t>
struct indexer_t
{
  list_t& list;
  constexpr indexer_t(list_t& list) : list(list) {}
  constexpr indexer_t(const indexer_t&) = default;
  template<class int_t>
  auto operator()(int_t&& index) const -> decltype(list[index])
  {
    return list[index];
  }
};
template<class T>
struct indexer_t<T*>
{
  T* const ptr;
  constexpr indexer_t(T* const ptr) : ptr(ptr) {}
  indexer_t(const indexer_t&) = default;
  template<class int_t>
  T& operator()(int_t&& index) const
  {
    return ptr[index];
  }
};
template<class T>
indexer_t<T> make_indexer(T& list)
{
  return indexer_t<T>(list);
}
template<class T>
indexer_t<T*> make_indexer(T* ptr)
{
  return indexer_t<T*>(ptr);
}

template<class int_t = int>
class invoke_counter_t
{
  mutable int_t count;

 public:
  invoke_counter_t() : count{} {}
  invoke_counter_t(int_t i) : count{i} {}
  invoke_counter_t(const invoke_counter_t&) = default;
  invoke_counter_t(invoke_counter_t&&)      = default;
  invoke_counter_t& operator=(const invoke_counter_t&) = default;
  invoke_counter_t& operator=(invoke_counter_t&&) = default;
  void reset(const int_t& value = 0) { count = value; }
  int_t& get_count() const { return count; }
  auto operator()() const -> decltype(count++) { return count++; }
  auto increment() const -> decltype(count++) { return count++; }
};
template<class int_t>
invoke_counter_t<int_t> make_counter(int_t initial = 0)
{
  return invoke_counter_t<int_t>{initial};
}

template<class list_t, class int_t, class F, class G>
void increment_and_call_on_matching(list_t&& list,
                                    int_t min,
                                    int_t max,
                                    F&& on_match,
                                    G&& on_no_match)
{
  for (auto& val : list)
  {
    for (; min < val && min <= max; ++min)
    {
      on_no_match(min);
    }
    on_match(min);
    ++min;
  }
  for (; min <= max; ++min)
  {
    on_no_match(min);
  }
}

template<class list_t, class F>
void remap(list_t& list, F&& func)
{
  for (auto& elem : list)
  {
    elem = func(elem);
  }
}
}   // namespace quanta
}   // namespace nimble
