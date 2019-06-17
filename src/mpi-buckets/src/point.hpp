#pragma once
#include "defines.hpp"
#include "meta.hpp"

template <class T>
struct Point3
{
    T x, y, z;

    template <class Func>
    constexpr auto map(Func&& func) const noexcept(noexcept(func(x)))
        -> Point3<OutputType<Func, T const&>>
    {
        return {func(x), func(y), func(z)};
    }
    template <class Func>
    CONSTEXPR_CPP14 auto map(Func&& func) noexcept(noexcept(func(x)))
        -> Point3<OutputType<Func, T&>>
    {
        return {func(x), func(y), func(z)};
    }
    template <class Func>
    constexpr auto apply(Func&& func) const noexcept(noexcept(func(x, y, z)))
        -> OutputType<Func, T const&, T const&, T const&>
    {
        return func(x, y, z);
    }
    
    template <class Func>
    CONSTEXPR_CPP14 auto apply(Func&& func) noexcept(noexcept(func(x, y, z)))
        -> OutputType<Func, T&, T&, T&>
    {
        return func(x, y, z);
    }
};

using Point3d = Point3<double>;
