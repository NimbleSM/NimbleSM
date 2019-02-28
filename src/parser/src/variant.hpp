#include <array>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>

template <class...>
class TypeSet
{
   public:
    constexpr static size_t size = 0;
};
struct ForwardWrapper
{
    template <class T>
    static ForwardWrapper wrap(T&& thing)
    {
        using PureT = typename std::decay<T>::type;
        return ForwardWrapper{const_cast<PureT*>(&thing)};
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
    Union(T const& thing) : value{thing} {}
    Union(T&& thing) : value{std::move(thing)} {}

    T    value;
    void init(T const& new_value) { new (&value) T(new_value); }
    void assign(T const& new_value) { value = new_value; }
    ~Union(){};
};
// General case
template <class T, class... Rest>
union Union<T, Rest...>
{
    T              value;
    Union<Rest...> rest;
    Union() {}
    template <class Q>
    Union(Q&& thing) : rest{std::forward<decltype(thing)>(thing)}
    {
    }
    Union(T const& thing) : value{thing} {}
    Union(T&& thing) : value{std::move(thing)} {}
    void init(T const& new_value) { new (&value) T(new_value); }
    void assign(T const& new_value) { value = new_value; }
    template <class Q>
    void init(Q const& new_value)
    {
        rest.init(new_value);
    }
    template <class Q>
    void assign(Q const& new_value)
    {
        rest.assign(new_value);
    }
    ~Union(){};
};

template <class T>
void destroy_as(void* i)
{
    ((T*)i)->~T();
}

template <class... T>
struct variant
{
    Union<T...> value;
    size_t      index;
    template <class Q>
    variant(Q&& thing)
      : value{std::forward<decltype(thing)>(thing)}
      , index{IndexLookup<typename std::remove_reference<Q>::type, T...>::index}
    {
        static_assert(IndexLookup<typename std::remove_reference<Q>::type, T...>::index
                          < sizeof...(T),
                      "Bad type");
    }
    variant(variant const& v) = delete;
    variant(variant&& other)  = delete;
    template <class Q>
    auto operator=(Q&& thing) -> variant&
    {
        using Qplain = typename std::remove_reference<Q>::type;
        constexpr static size_t Q_index = IndexLookup<Qplain, T...>::index;
        static_assert(Q_index < sizeof...(T), "Item not in variant");
        if (Q_index == index)
        {
            ((Qplain&)value) = thing;
        }
        else
        {
            ~variant();
            value.init(thing);
        }
    }
    ~variant()
    {
        using destructor_t = void (*)(void*);
        constexpr static destructor_t lookup[]{destroy_as<T>...};
        lookup[index](&value);
    }
};

template <class T>
using Wrap_t = ForwardWrapper;

template <class F, class... Conv>
struct FuncGroup
{
    static void convert_apply(F func, Wrap_t<Conv>... value)
    {
        func(value.template as<Conv>()...);
    }
};


template <class...>
class VTable;
template <class Func, template <class...> class ParamList, class... T>
class VTable<Func, ParamList<T...>>
{
   public:
    using table_t = void (*)(Func, Wrap_t<T>...);
    constexpr static auto table() -> table_t
    {
        return FuncGroup<Func, T...>::convert_apply;
    }
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
    func_ptr(func, ForwardWrapper::wrap(variants.value)...);
}
