#pragma once
#include <cstddef>

namespace ut
{
template <class...>
class TypeList
{
   public:
    constexpr static size_t size = 0;
};
template <class A, class... B>
class TypeList<A, B...>
{
   public:
    using Head                   = A;
    using Tail                   = TypeList<B...>;
    constexpr static size_t size = 1 + sizeof...(B);
};
}  // namespace ut
