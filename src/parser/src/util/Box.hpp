#pragma once
#include <type_traits>

namespace ut
{
struct Box
{
    template <class T>
    static Box box(T&& thing)
    {
        using PureT = typename std::decay<T>::type;
        return Box{const_cast<PureT*>(&thing)};
    }
    void* addr;
    template <class T>
    auto as() -> T
    {
        using PureT = typename std::remove_reference<T>::type;
        return const_cast<T>(*reinterpret_cast<PureT*>(addr));
    }
};
template <class T>
using Box_t = Box;
}  // namespace ut
