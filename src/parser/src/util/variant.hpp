#include "util/Box.hpp"
#include "util/TypeList.hpp"
#include "util/Union.hpp"
#include "util/VTable.hpp"

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>


namespace ut
{
template <class T>
constexpr T const& recursive_index(T const& value)
{
    return value;
}
template <class T, class... Int>
constexpr auto recursive_index(T const& arr, size_t index, Int... rest)
    -> decltype(recursive_index(arr[index], rest...))
{
    return recursive_index(arr[index], rest...);
}


template <class Target, class...>
class IndexLookup
{
   public:
    constexpr static size_t index = 0;
};
template <class Target, class... Rest>
class IndexLookup<Target, Target, Rest...>
{
   public:
    constexpr static size_t index = 0;
};
template <class Target, class Alt, class... Rest>
class IndexLookup<Target, Alt, Rest...>
{
   public:
    constexpr static size_t index = 1 + IndexLookup<Target, Rest...>::index;
};

template <class T>
void del(Box x)
{
    x.as<T&>().~T();
}
template <class U, class T>
void init_as(U& u, Box t)
{
    u.init(t.as<T>());
}
template <class U, class T>
void assign_as(U& u, Box t)
{
    u.assign(t.as<T>());
}

template <class... T>
struct variant
{
    using assigner_t = void (*)(Union<T...>&, Box);
    Union<T...>           value;
    constexpr static bool trivially_destructible
        = Union<T...>::trivially_destructible;
    constexpr static bool call_destructor = !trivially_destructible;
    size_t                index;
    template <class Q>
    variant(Q&& thing)
      : value{std::forward<decltype(thing)>(thing)}
      , index{IndexLookup<typename std::remove_reference<Q>::type, T...>::index}
    {
        static_assert(IndexLookup<typename std::remove_reference<Q>::type, T...>::index
                          < sizeof...(T),
                      "Bad type");
    }
    variant(variant const& other)
    {
        constexpr static assigner_t lookup[]{init_as<Union<T...>, T const&>...};
        lookup[other.index](value, Box{&other.value});
        index = other.index;
    }
    variant(variant&& other)
    {
        constexpr static assigner_t lookup[]{init_as<Union<T...>, T&&>...};
        lookup[other.index](value, Box{&other.value});
        index = other.index;
    }
    auto operator=(variant& thing) -> variant&
    {
        constexpr static assigner_t init_lookup[]{init_as<Union<T...>, T&>...};
        constexpr static assigner_t assign_lookup[]{assign_as<Union<T...>, T&>...};

        if (index == thing.index)
        {
            assign_lookup[index](value, Box{&thing.value});
        }
        else
        {
            this->~variant();
            init_lookup[index](value, Box{&thing.value});
        }
        return *this;
    }
    auto operator=(variant const& thing) -> variant&
    {
        constexpr static assigner_t init_lookup[]{init_as<Union<T...>, T const&>...};
        constexpr static assigner_t assign_lookup[]{
            assign_as<Union<T...>, T const&>...};

        if (index == thing.index)
        {
            assign_lookup[index](value, Box{&thing.value});
        }
        else
        {
            this->~variant();
            init_lookup[index](value, Box{&thing.value});
        }
        return *this;
    }
    auto operator=(variant&& thing) -> variant&
    {
        constexpr static assigner_t init_lookup[]{init_as<Union<T...>, T&&>...};
        constexpr static assigner_t assign_lookup[]{
            assign_as<Union<T...>, T const&>...};

        if (index == thing.index)
        {
            assign_lookup[index](value, Box{&thing.value});
        }
        else
        {
            this->~variant();
            init_lookup[index](value, Box{&thing.value});
        }
        return *this;
    }
    template <class Q>
    auto operator=(Q&& thing) -> variant&
    {
        using Qplain = typename std::remove_reference<Q>::type;
        constexpr static size_t Q_index = IndexLookup<Qplain, T...>::index;
        static_assert(Q_index < sizeof...(T), "Item not in variant");
        if (Q_index == index)
        {
            ((Qplain&)value) = std::forward<decltype(thing)>(thing);
        }
        else
        {
            this->~variant();
            value.init(std::forward<decltype(thing)>(thing));
            index = Q_index;
        }
        return *this;
    }
    ~variant()
    {
        if (call_destructor)
        {
            using destructor_t = void (*)(Box);
            constexpr static destructor_t lookup[]{del<T>...};
            lookup[index](Box::box(&value));
        }
    }
    template <class... In>
    void operator()(In&&... inputs)
    {
        using caller_t = void (*)(Box, decltype(inputs)...);
        constexpr static caller_t lookup[]{FG<T, decltype(inputs)...>::app_boxed...};
        lookup[index](std::forward<decltype(inputs)>(inputs)...);
    }
};



template <class F, class... V>
void visit(F&& func, V&&... variants)
{
    constexpr static auto lookup
        = VisitorDispatchTable<F&, TypeList<>, decltype(variants)...>::table();
    const auto func_ptr = recursive_index(lookup, variants.index...);
    func_ptr(func, Box::box(variants.value)...);
}
struct UniversalInvoker
{
    template <class F, class... Args>
    constexpr auto operator()(F&& func, Args&&... args) const
        -> decltype(func(std::forward<decltype(args)>(args)...))
    {
        return func(std::forward<decltype(args)>(args)...);
    }
};
constexpr auto invoke = UniversalInvoker{};
template <class... F, class... V>
void visit(variant<F...>&& func, V&&... variants)
{
    visit(invoke, std::move(func), std::forward<decltype(variants)>(variants)...);
}
template <class... F, class... V>
void visit(variant<F...>& func, V&&... variants)
{
    visit(invoke, func, std::forward<decltype(variants)>(variants)...);
}
template <class... F, class... V>
void visit(variant<F...> const& func, V&&... variants)
{
    visit(invoke, func, std::forward<decltype(variants)>(variants)...);
}
}  // namespace ut
using ut::variant;
using ut::visit;
