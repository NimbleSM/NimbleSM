#pragma once
#include "ext-includes.hpp"

template <class First, class... Rest>
using FirstOf = First;

template <class T>
T iof();
template <class Range>
using ElemType = typename std::decay<decltype(*std::begin(iof<Range>()))>::type;

template <class Func, class... Input>
using OutputType = decltype(iof<Func>()(iof<Input>()...));

template <class Func, class... Input>
using DecayedOutputType = typename std::decay<OutputType<Func, Input...>>::type;

template <class T>
struct StoreTypeHelper
{
    using type = T;
};
template <class T>
struct StoreTypeHelper<T&&>
{
    using type = T;
};
template <class T>
struct StoreTypeHelper<T&>
{
    using type = T&;
};

template <class T>
using StoreType = typename StoreTypeHelper<T>::type;

template <class T>
auto Preserve(T&& boi) -> StoreType<T>
{
    return boi;
}

template <class T>
constexpr T* NullOf()
{
    return nullptr;
}

template <class Range, class Func, template <class...> class map = std::unordered_map>
using GatherByMap_t
    = map<DecayedOutputType<Func, ElemType<Range>>, std::vector<ElemType<Range>>>;
template <class Range, class Func, class Transform, template <class...> class map = std::unordered_map>
using GatherTransformByMap_t
    = map<DecayedOutputType<Func, ElemType<Range>>,
          std::vector<DecayedOutputType<Transform, ElemType<Range>>>>;


template <class Range, class Func>
auto GatherBy(Range&& range, Func func) -> GatherByMap_t<Range, Func>
{
    GatherByMap_t<Range, Func> dest;
    for (auto&& elem : range)
        dest[func(elem)].push_back(elem);

    return dest;
}
template <class Range, class Func, class Transform>
auto GatherTransformBy(Range&& range, Func func, Transform transform)
    -> GatherTransformByMap_t<Range, Func, Transform>
{
    GatherTransformByMap_t<Range, Func, Transform> dest;
    for (auto&& elem : range)
        dest[func(elem)].push_back(transform(elem));

    return dest;
}

template <class T, class index_t = int>
auto invert(std::vector<T> const& vals) -> std::unordered_map<T, index_t>
{
    std::unordered_map<T, index_t> map;
    index_t                        i = 0;
    for (auto const& val : vals)
    {
        map.emplace(val, i++);
    }
    return map;
}
