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
void destruct_as(void* i)
{
    ((T*)i)->~T();
}

template <class... T>
struct Variant;

template <class... T>
struct Variant
{
    Union<T...> value;
    size_t      index;

    ~Variant()
    {
        using destructor_t = void (*)(void*);
        constexpr static destructor_t lookup[]{destruct_as<T>...};
        lookup[index](&value);
    }
};

void arrayTest() { int myArray[2][2]{{1, 2}, {3, 4}}; }
template <class... T>
using void_t = void;
template <class F, class... Conv>
void convert_apply(F func, void_t<Conv>*... value)
{
    func(*(Conv*)value...);
}
template <class F, class... T>
void visit(F&& f, Variant<T...>& v)
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
void visit(F&& f, Variant<T...>& v, Variant<Q...>& q)
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
// void visit(F&& f, Variant<A, B, C, D>& v) {
//    switch(v.index) {
//        case 0: f((A&)v.value); break;
//        case 1: f((B&)v.value); break;
//        case 2: f((C&)v.value); break;
//        case 3: f((D&)v.value); break;
//        default: break;
//    }
//}

struct S0;
struct S1;
struct S2;
struct S3;
struct S4;
struct S5;
struct S6;
struct S7;
struct S8;
struct S9;
template <class... T>
void f(T*...);
struct FCaller
{
    template <class... T>
    auto operator()(T&&... thang) const -> void
    {
        f(thang...);
    }
};
constexpr FCaller func;

void foo(FCaller& func, Variant<S0*, S1*, S2*, S3*>& v) { visit(func, v); }
void foo(FCaller&                     func,
         Variant<S0*, S1*, S2*, S3*>& v,
         Variant<S0*, S1*, S2*, S3*>& q)
{
    visit(func, v, q);
}
#include <string>
void foo(int);
void foo(std::string&);
