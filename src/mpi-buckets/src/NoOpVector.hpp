#include <memory>
#include <vector>

struct EmptyValue
{
    EmptyValue operator=(EmptyValue value) noexcept { return {}; }
    EmptyValue operator=(EmptyValue value) const noexcept { return {}; }

    constexpr bool operator==(EmptyValue) const noexcept { return true; }
    constexpr bool operator!=(EmptyValue) const noexcept { return false; }
    constexpr bool operator<(EmptyValue) const noexcept { return false; }
    constexpr bool operator>(EmptyValue) const noexcept { return false; }
    constexpr bool operator<=(EmptyValue) const noexcept { return true; }
    constexpr bool operator>=(EmptyValue) const noexcept { return true; }
};

struct EmptyValueIterator
{
    size_t              index;
    EmptyValue          operator*() const { return {}; }
    EmptyValueIterator& operator++()
    {
        ++index;
        return *this;
    }
    EmptyValueIterator operator++(int) const { return {index + 1}; }
};
template <class Allocator>
class NoopVector : private Allocator
{
   public:
    using value_type      = EmptyValue;
    using allocator_type  = Allocator;
    using size_type       = size_t;
    using difference_type = std::ptrdiff_t;
    using reference       = EmptyValue;
    using const_reference = const EmptyValue;
    using pointer         = typename std::allocator_traits<Allocator>::pointer;
    using const_pointer = typename std::allocator_traits<Allocator>::const_pointer;
    using iterator      = EmptyValue*;
    using const_iterator         = EmptyValue const&;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    NoopVector()                  = default;
    NoopVector(const NoopVector&) = default;
    NoopVector(NoopVector&&)      = default;
    NoopVector(NoopVector const&, Allocator const& alloc) : Allocator(alloc) {}
    NoopVector(NoopVector&&, Allocator const& alloc) : Allocator(alloc) {}
    NoopVector(std::initializer_list<EmptyValue> init,
               const Allocator&                  alloc = Allocator())
    {
    }
    NoopVector(const Allocator& alloc) : Allocator(alloc) {}

    NoopVector operator=(NoopVector) const { return {}; }
    void       assign(size_type count, EmptyValue const& value) {}
    template <class InputIt>
    void assign(InputIt first, InputIt last)
    {
        for (; first != last; ++first)
        {
            EmptyValue value(*first);
        }
    }
    void           assign(std::initializer_list<EmptyValue> list) {}
    allocator_type get_allocator() const { return Allocator(*this); }

    EmptyValue       at(size_type pos) { return {}; }
    EmptyValue       operator[](size_type) noexcept { return {}; }
    EmptyValue const operator[](size_type) const noexcept { return {}; }
    EmptyValue       front() {}
    EmptyValue       back() {}
    EmptyValue*      data() { return nullptr; }

    iterator       begin() noexcept { return EmptyValue(); }
    const_iterator begin() const noexcept { return EmptyValue(); }
    iterator       end() noexcept { return EmptyValue(); }
    const_iterator end() const noexcept { return EmptyValue(); }


    size_type size() const noexcept { return 0; }
};
