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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <random>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#define TESTPRINT(x) \
  std::cout << "Running " #x " ... " << ((x) ? "(success)" : "(faliure)") << std::endl;
namespace nimble
{
namespace quanta
{
template<class T>
T declval_exact() noexcept;

template<class T>
typename std::remove_reference<T>::type& declref() noexcept;

// Gets the literal type of a container when it's iterator is dereferenced
// for example, 'iterated_t<std::vector<int>>' is of type 'int&', and
// 'iterated_t<const std::vector<int>>' is of type 'const int&'
template<class list_t>
using iterated_t = decltype(*std::begin(declref<list_t>()));

// Gets the unmodified type of a container when it's iterator is defererenced
// For example, 'elem_t<std::vector<int>>' is of type 'int'
template<class list_t>
using elem_t = typename std::decay<iterated_t<list_t>>::type;

template<class F, class... Args>
using return_t = decltype(declref<F>()(declval_exact<Args>()...));

template<class F, class... list_t>
using transformed_iterated_t = return_t<F, iterated_t<list_t>...>;

template<class F, class... list_t>
using transformed_elem_t = return_t<F, elem_t<list_t>...>;

template<class T>
auto len(T&& obj) -> decltype(obj.size())
{
  return obj.size();
}

// Does nothing with the inputs it's given
// Used in place of C++17 fold expression with comma operator
template<class... In>
void devour(In&&...)
{
}

// Shuffles the input lists as a group
// So if A[i] moves to A[j], then B[i] also moves to B[j]
template<class gen_t, class... Iterator_t>
void shuffle_together(size_t count, gen_t&& gen, Iterator_t... lists)
{
  typedef std::tuple<decltype(*lists)...> tuple_t;
  std::vector<tuple_t> groups;
  groups.reserve(count);
  for (size_t i = 0; i < count; ++i)
  {
    groups.emplace_back(*lists...);
    devour(++lists...);
  }
  std::shuffle(groups.begin(), groups.end(), gen);
}


//void radix_order256(uint32_t* array, uint32_t* ordering, uint32_t* buffer, uint32_t count);

// Takes a list of ids and remaps it in-place so that they're all consecutive
// Does so based on an ordering
/*Example:
                   0  1  2  3  4  5  6   7
  Suppose array = {1, 3, 2, 9, 3, 9, 10, 2}
  Then the order will be {0, 2, 7, 1, 4, 3, 5, 6}
  After calling this function,
  array becomes   {0, 2, 1, 3, 2, 3, 4, 1}
*/

template<class id_array_t, class order_t>
void remap_ids_with_ordering(id_array_t& array, const order_t& order)
{
  typedef elem_t<id_array_t> id_t;

  auto scan = std::begin(order);
  auto end  = std::end(order);
  if (scan == end)
    return;

  id_t reference = array[*scan];
  id_t counter   = 0;
  array[*scan]   = counter;
  ++scan;
  while (scan != end)
  {
    id_t& id = array[*scan];
    ++scan;
    if (id == reference)
    {
      id = counter;
    }
    else
    {
      reference = id;
      id        = ++counter;
    }
  }
}

//void PackIDs(std::vector<int>& id_array);

// Checks that listA[i] == listA[j] iff listB[i] == listB[j]
template<class list1_t, class list2_t>
bool is_bijected(const list1_t& listA, const list2_t& listB)
{
  std::unordered_map<elem_t<list1_t>, elem_t<list2_t>> A_to_B;
  std::unordered_map<elem_t<list2_t>, elem_t<list1_t>> B_to_A;
  auto itA    = std::begin(listA);
  auto itB    = std::begin(listB);
  auto endOfA = std::end(listA);
  auto endOfB = std::end(listB);
  while (itA != endOfA && itB != endOfB)
  {
    auto& a        = *itA;
    auto& b        = *itB;
    auto A2B_index = A_to_B.find(a);
    auto B2A_index = B_to_A.find(b);
    bool a_foundQ  = A2B_index != A_to_B.end();
    bool b_foundQ  = B2A_index != B_to_A.end();
    // If one was found but the other wasn't, there's a mismatch and it's not
    // bijective
    if (a_foundQ != b_foundQ)
      return false;

    if (a_foundQ && b_foundQ)
    {
      // If both were found, A_to_B[a] must map to b, and B_to_A[b] must map to
      // a.

      // This dereferences the iterators to check that
      if (A2B_index->second != b || B2A_index->second != a)
        return false;
    }
    else
    {
      // If neither of them were found, record the new mapping from both sides
      // So that a maps to b and vice versa
      A_to_B[a] = b;
      B_to_A[b] = a;
    }
    ++itA;
    ++itB;
  }
  if (itA != endOfA || itB != endOfB)
    return false;
  return true;
}
template<class id_array_t>
auto remap_ids_with_map(const id_array_t& ids) -> std::vector<elem_t<id_array_t>>
{
  typedef elem_t<id_array_t> int_t;
  std::vector<int_t> new_ids;
  new_ids.reserve(ids.size());
  int_t num_ids = 0;
  std::unordered_map<elem_t<id_array_t>, int_t> mapping;
  for (auto& elem : ids)
  {
    auto position_in_map = mapping.find(elem);
    int_t id_to_add;
    if (position_in_map == mapping.end())
    {
      mapping[elem] = num_ids;
      new_ids.push_back(num_ids);
      ++num_ids;
    }
    else
    {
      new_ids.push_back(position_in_map->second);
    }
  }
  return new_ids;
}
template<class id_array_t>
id_array_t& remap_ids_with_map_in_place(id_array_t& ids)
{
  typedef elem_t<id_array_t> int_t;
  int_t num_ids = 0;
  std::unordered_map<int_t, int_t> mapping;
  mapping.reserve(ids.size());
  for (auto& elem : ids)
  {
    auto position_in_map = mapping.find(elem);
    if (position_in_map == mapping.end())
    {
      mapping[elem] = num_ids;
      elem          = num_ids;
      ++num_ids;
    }
    else
    {
      elem = position_in_map->second;
    }
  }
  return ids;
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

template<class list_t, class F>
void remap(list_t& list, F&& func)
{
  for (auto& elem : list)
  {
    elem = func(elem);
  }
}

#ifdef TEST_QUANTA
#include <chrono>

namespace tests
{

// If test_shuffle_together worked, then A and B will be shuffled identically
// And they'll have actually been shuffled, so they should be distinct from C
template<class gen_t>
bool test_shuffle_together(size_t count, gen_t&& gen)
{
  if (count < 20)
  {
    std::cerr << "test_shuffle_together(" << count
              << ", [gen]):\n"
                 "  Possible test faliure due to small input size.\n"
                 "  Test assumes that lists are changed after being shuffled,\n"
                 "  Although there is a small probability that the shuffle routine\n"
                 "  Worked correctly but didn't make any changes to the lists. This\n"
                 "  Probability is inversely proportional to the factorial of "
                 "count,\n"
                 "  Which is why it's recommended to use an input size greater than "
                 "20.";
  }
  // Initializes all the lists so that they'll be {0, 1, 2, 3, ...}
  std::vector<size_t> A(count), B(count), C(count);
  for (size_t i = 0; i < count; ++i)
  {
    A[i] = i;
    B[i] = i;
    C[i] = i;
  }
  // Shuffles A and B together
  shuffle_together(count, gen, A.begin(), B.begin());
  // End result should be that A == B but A != C
  // Since A and B should have been shuffled identically,
  // But C wasn't shuffled
  // For vectors, the == operator does an element-wise comparison
  return A == B && A != C;
}

template<class gen_t>
bool test_is_bijected(size_t count, gen_t&& gen)
{
  std::vector<size_t> A(count), B(count);
  // Ensures at least 1 collision
  std::uniform_int_distribution<size_t> dist(1, count - 1);
  // mask cannot be 0
  size_t mask = dist(gen);
  for (size_t i = 0; i < count; ++i)
  {
    size_t val = dist(gen);
    A[i]       = val;
    B[i]       = val ^ mask;
  }
  // test_positive checks that A and B are confirmed to be bijective
  bool test_positive_case = is_bijected(A, B) && is_bijected(B, A);

  size_t index_to_change = 0;
  while (index_to_change == 0)
  {
    index_to_change = dist(gen);
  }
  // Breaks the bijection so that B[0] == B[index_to_change] but A[0] !=
  // A[index_to_change]
  B[0] = B[index_to_change];
  A[0] = A[index_to_change] ^ mask;
  // test_negative checks that A and B are confirmed NOT to be bijective
  bool test_negative_case = !is_bijected(A, B) && !is_bijected(B, A);
  return test_positive_case && test_negative_case;
}

void test_all()
{
   std::random_device gen{};
   TESTPRINT(tests::test_shuffle_together(1000000, gen));
   TESTPRINT(tests::test_is_bijected(1000000, gen));
 }
}   // namespace tests

auto now = []() { return std::chrono::high_resolution_clock::now(); };
int main(int argc, char** argv)
{
  std::vector<int> vect(200000000);
  for (uint64_t i = 0, seed = 0; i < vect.size(); ++i)
  {
    seed = seed * 2039809803534549 + 809380980345342;
    seed >>= 16;
    vect[i] = seed;
  }
  std::cout << vect.back() << std::endl;
  auto t0 = now();
  if (argc == 1)
  {
    std::cout << "PackIDs()" << std::endl;
    nimble::quanta::PackIDs(vect);
  }
  else if (argc == 2)
  {
    std::cout << "RemapIDs()" << std::endl;
    vect = nimble::quanta::remap_ids_with_map(vect);
  }
  else
  {
    std::cout << "RemapInPlace()" << std::endl;
    nimble::quanta::remap_ids_with_map_in_place(vect);
  }
  auto t1 = now();
  std::cout << (t1 - t0).count() * 1e-9 << std::endl;
  std::cout << vect.back() << std::endl;
}
#endif

}   // namespace quanta
}   // namespace nimble

#endif // NIMBLE_QUANTA_H
