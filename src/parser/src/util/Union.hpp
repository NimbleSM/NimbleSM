#pragma once
#include <cstddef>
#include <new>
#include <type_traits>

namespace ut
{
template <class...>
union Union;
// Base case
template <class T>
union Union<T>
{
    constexpr static bool trivially_destructible
        = std::is_trivially_destructible<T>::value;
    Union(T const& thing) : value{thing} {}
    Union(T&& thing) : value{std::move(thing)} {}

    T    value;
    void init(T const& new_value) { new (&value) T(new_value); }
    void assign(T const& new_value) { value = new_value; }
    void init(T&& new_value) { new (&value) T(std::move(new_value)); }
    void assign(T&& new_value) { value = std::move(new_value); }
    ~Union(){};
};
// General case
template <class T, class... Rest>
union Union<T, Rest...>
{
    constexpr static bool trivially_destructible
        = std::is_trivially_destructible<T>::value
        && Union<Rest...>::trivially_destructible;

    T              value;
    Union<Rest...> rest;
    Union() {}
    template <class Q>
    Union(Q&& thing) : rest{std::forward<decltype(thing)>(thing)}
    {
    }
    Union(T& thing) : value{thing} {}
    Union(T const& thing) : value{thing} {}
    Union(T&& thing) : value{std::move(thing)} {}
    void init(T&& new_value) { new (&value) T(std::move(new_value)); }
    void init(T const& new_value) { new (&value) T(new_value); }
    void init(T& new_value) { new (&value) T(new_value); }
    void assign(T&& new_value) { value = std::move(new_value); }
    void assign(T const& new_value) { value = new_value; }
    void assign(T& new_value) { value = new_value; }
    template <class Q>
    void init(Q&& new_value)
    {
        rest.init(std::forward<decltype(new_value)>(new_value));
    }
    template <class Q>
    void assign(Q&& new_value)
    {
        rest.assign(std::forward<decltype(new_value)>(new_value));
    }
    ~Union(){};
};
}  // namespace ut
