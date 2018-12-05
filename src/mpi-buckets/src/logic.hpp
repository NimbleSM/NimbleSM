#pragma once
constexpr auto all_true(bool first) noexcept -> bool { return first; }
template <class... T>
constexpr auto all_true(bool first, bool second, T... rest) noexcept -> bool {
    return all_true(first && second, rest...);
}

template <class... T>
constexpr static auto any_true(bool first, bool second, T... rest) noexcept -> bool {
    return all_true(first || second, rest...);
}
constexpr static auto any_true(bool first) noexcept -> bool { return first; }

template <class Item, class... In>
constexpr static auto equals_anyof(Item&& item, In&&... values) -> bool {
    return any_true((item == values)...);
}
