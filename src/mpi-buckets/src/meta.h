#pragma once
#include <type_traits>
#include <utility>
template <class First, class... Rest>
using FirstOf = First;

template <class T>
T iof();
template <class Range>
using ElemType = typename std::decay<decltype(std::begin(iof<Range>()))>::type;

template <class Func, class... Input>
using OutputType = decltype(iof<Func>()(iof<Input>()...));

template <class Func, class... Input>
using DecayedOutputType = typename std::decay<OutputType<Func, Input...> >::type;

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

template<class T> constexpr T* NullOf() {
    return nullptr; 
}
