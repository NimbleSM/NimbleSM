#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

template <class F, class... V>
void visit(F&& func, V&&... variants);

template <class...>
class TypeSet
{
   public:
    constexpr static size_t size = 0;
};
struct Box
{
    template <class T>
    static Box wrap(T&& thing)
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
using Wrap_t = Box;

template <class F, class... Conv>
struct FG
{
    static void app(F func, Wrap_t<Conv>... value)
    {
        func(value.template as<Conv>()...);
    }
    static void app_boxed(Box func, Conv... values)
    {
        func.as<F>()(std::forward<Conv>(values)...);
    }
};
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

template <class A, class... B>
class TypeSet<A, B...>
{
   public:
    using Head                   = A;
    using Tail                   = TypeSet<B...>;
    constexpr static size_t size = 1 + sizeof...(B);
};

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
    Union(T const& thing) : value{thing} {}
    Union(T&& thing) : value{std::move(thing)} {}
    void init(T&& new_value) { new (&value) T(std::move(new_value)); }
    void assign(T&& new_value) { value = std::move(new_value); }
    void init(T const& new_value) { new (&value) T(new_value); }
    void assign(T const& new_value) { value = new_value; }
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

template <class T>
void del(void* i)
{
    ((T*)i)->~T();
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
            using destructor_t = void (*)(void*);
            constexpr static destructor_t lookup[]{del<T>...};
            lookup[index](&value);
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



template <class...>
class VTable;
template <class Func, template <class...> class ParamList, class... T>
class VTable<Func, ParamList<T...>>
{
   public:
    using table_t = void (*)(Func, Wrap_t<T>...);
    constexpr static auto table() -> table_t { return FG<Func, T...>::app; }
};
template <class Func,
          template <class...> class ParamList,
          template <class...> class Variant,
          class... T,
          class... U,
          class... V>
class VTable<Func, ParamList<T...>, Variant<U...> const&, V...>
{
   public:
    using ElemType =
        typename TypeSet<typename VTable<Func, ParamList<T..., U>, V...>::table_t...>::Head;
    using table_t = std::array<ElemType, sizeof...(U)>;
    constexpr static auto table() -> table_t
    {
        return {VTable<Func, ParamList<T..., U const&>, V...>::table()...};
    }
};
template <class Func,
          template <class...> class ParamList,
          template <class...> class Variant,
          class... T,
          class... U,
          class... V>
class VTable<Func, ParamList<T...>, Variant<U...>&, V...>
{
   public:
    using ElemType =
        typename TypeSet<typename VTable<Func, ParamList<T..., U>, V...>::table_t...>::Head;
    using table_t = std::array<ElemType, sizeof...(U)>;
    constexpr static auto table() -> table_t
    {
        return {VTable<Func, ParamList<T..., U&>, V...>::table()...};
    }
};
template <class Func,
          template <class...> class ParamList,
          template <class...> class Variant,
          class... T,
          class... U,
          class... V>
class VTable<Func, ParamList<T...>, Variant<U...>&&, V...>
{
   public:
    using ElemType =
        typename TypeSet<typename VTable<Func, ParamList<T..., U>, V...>::table_t...>::Head;
    using table_t = std::array<ElemType, sizeof...(U)>;
    constexpr static auto table() -> table_t
    {
        return {VTable<Func, ParamList<T..., U&&>, V...>::table()...};
    }
};
template <class F, class... V>
void visit(F&& func, V&&... variants)
{
    constexpr static auto lookup
        = VTable<F&, TypeSet<>, decltype(variants)...>::table();
    const auto func_ptr = recursive_index(lookup, variants.index...);
    func_ptr(func, Box::wrap(variants.value)...);
}
