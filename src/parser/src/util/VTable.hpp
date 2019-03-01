#pragma once
#include <array>
#include "util/Box.hpp"
#include "util/TypeList.hpp"

namespace ut
{
template <class F, class... Conv>
struct FG
{
    static void app(F func, Box_t<Conv>... value) { func(value.as<Conv>()...); }
    static void app_boxed(Box func, Conv... values)
    {
        func.as<F>()(std::forward<Conv>(values)...);
    }
};
template <class...>
class VisitorDispatchTable;
template <class Visitor, template <class...> class ParamList, class... T>
class VisitorDispatchTable<Visitor, ParamList<T...>>
{
   public:
    using table_t = void (*)(Visitor, Box_t<T>...);
    constexpr static auto table() -> table_t { return FG<Visitor, T...>::app; }
};
template <class Visitor,
          template <class...> class ParamList,
          template <class...> class Variant,
          class... T,
          class... U,
          class... V>
class VisitorDispatchTable<Visitor, ParamList<T...>, Variant<U...> const&, V...>
{
   public:
    using ElemType = typename TypeList<
        typename VisitorDispatchTable<Visitor, ParamList<T..., U>, V...>::table_t...>::Head;
    using table_t = std::array<ElemType, sizeof...(U)>;
    constexpr static auto table() -> table_t
    {
        return {
            VisitorDispatchTable<Visitor, TypeList<T..., U const&>, V...>::table()...};
    }
};
template <class Visitor,
          template <class...> class ParamList,
          template <class...> class Variant,
          class... T,
          class... U,
          class... V>
class VisitorDispatchTable<Visitor, ParamList<T...>, Variant<U...>&, V...>
{
   public:
    using ElemType = typename TypeList<
        typename VisitorDispatchTable<Visitor, ParamList<T..., U>, V...>::table_t...>::Head;
    using table_t = std::array<ElemType, sizeof...(U)>;
    constexpr static auto table() -> table_t
    {
        return {VisitorDispatchTable<Visitor, TypeList<T..., U&>, V...>::table()...};
    }
};
template <class Visitor,
          template <class...> class ParamList,
          template <class...> class Variant,
          class... T,
          class... U,
          class... V>
class VisitorDispatchTable<Visitor, ParamList<T...>, Variant<U...>&&, V...>
{
   public:
    using ElemType = typename TypeList<
        typename VisitorDispatchTable<Visitor, ParamList<T..., U>, V...>::table_t...>::Head;
    using table_t = std::array<ElemType, sizeof...(U)>;
    constexpr static auto table() -> table_t
    {
        return {VisitorDispatchTable<Visitor, TypeList<T..., U&&>, V...>::table()...};
    }
};
}  // namespace ut
