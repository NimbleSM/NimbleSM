#pragma once
#include <utility>

template <class Number = double>
struct NumericRange
{
    Number min, max;
    constexpr NumericRange(Number center) noexcept
      : NumericRange(center, center)
    {
    }
    constexpr NumericRange(Number min, Number max) noexcept : min(min), max(max)
    {
    }
    constexpr NumericRange() : min(), max() {}
    constexpr NumericRange(NumericRange const&) = default;
    constexpr NumericRange(NumericRange&&)      = default;
    auto operator=(NumericRange const&) noexcept -> NumericRange& = default;
    auto operator=(NumericRange&&) noexcept -> NumericRange& = default;

    auto centerAt(Number val) & noexcept -> NumericRange&
    {
        min = val;
        max = val;
        return *this;
    }
    auto centerAt(Number val) && noexcept -> NumericRange&&
    {
        return std::move(centerAt(val));
    }
    auto include(Number val) & noexcept -> NumericRange&
    {
        if (val < min)
            min = val;
        else if (val > max)
            max = val;
        return *this;
    }
    auto include(Number val) && noexcept -> NumericRange&&
    {
        return std::move(include(val));
    }
    auto include(NumericRange const& range) & noexcept -> NumericRange&
    {
        if (range.min < min)
            min = range.min;
        if (range.max > max)
            max = range.max;
        return *this;
    }
    auto include(NumericRange const& range) && noexcept -> NumericRange&&
    {
        return std::move(include(range));
    }
    constexpr auto contains(Number val) const noexcept -> bool
    {
        return min <= val && val <= max;
    }

    constexpr auto intersects(NumericRange const& other) const noexcept -> bool
    {
        return min <= other.max && other.min <= max;
    }
};
