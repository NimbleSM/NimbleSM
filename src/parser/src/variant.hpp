#pragma once
#include <array>
#include <cstddef>
#include <type_traits>

template <class...>
union Union;
// Base case
template <class T>
union Union<T>
{
    T value;
    ~Union(){};
};
// General case
template <class T, class... Rest>
union Union<T, Rest...>
{
    T              value;
    Union<Rest...> rest;
    ~Union(){};
};

template <class T>
void destroy_as(void* i)
{
    ((T*)i)->~T();
}

template <class... T>
struct RawVariant
{
    Union<T...> value;
    size_t      index;

    ~RawVariant()
    {
        using destructor_t = void (*)(void*);
        constexpr static destructor_t lookup[]{destroy_as<T>...};
        lookup[index](&value);
    }
};

template <class... T>
struct variant : public RawVariant<T...>
{
    using RawVariant<T...>::RawVariant;
};

template <class... T>
using void_t = void;

template <class F, class... Conv>
void convert_apply(F func, void_t<Conv>*... value)
{
    func(*(Conv*)value...);
}

template <class F, class... T>
void visit(F&& f, RawVariant<T...>& v)
{
    using func_t = void (*)(F&, void*);
    constexpr static func_t lookup[]{convert_apply<F&, T>...};
    lookup[v.index](f, &v.value);
}
template <class F, class T, class... Conv>
constexpr auto convert_apply_2()
    -> std::array<void (*)(F&, void*, void*), sizeof...(Conv)>
{
    using func_t = void (*)(F&, void*, void*);
    return std::array<func_t, sizeof...(Conv)>{convert_apply<F, T, Conv>...};
}
template <class F, class... T, class... Q>
void visit(F&& f, RawVariant<T...>& v, RawVariant<Q...>& q)
{
    using func_t              = void (*)(F&, void*, void*);
    constexpr static size_t N = sizeof...(T), M = sizeof...(Q);
    if (v.index > N || q.index > M)
        return;
    constexpr static std::array<std::array<func_t, M>, N> lookup
        = {convert_apply_2<F, T, Q...>()...};
    lookup[v.index][q.index](f, &v.value, &q.value);
}
// template<class F, class A, class B, class C, class D>
// void visit(F&& f, RawVariant<A, B, C, D>& v) {
//    switch(v.index) {
//        case 0: f((A&)v.value); break;
//        case 1: f((B&)v.value); break;
//        case 2: f((C&)v.value); break;
//        case 3: f((D&)v.value); break;
//        default: break;
//    }
//}
