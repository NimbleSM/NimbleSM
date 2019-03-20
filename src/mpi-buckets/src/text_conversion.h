#pragma once
#include "BoundingBox.h"
#include "GridCell.h"
#include "GridIndex.h"
#include "meta.h"

#include <iostream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

template <class First, class... Items>
auto printSpaceDelimited(std::ostream& stream, First&& first, Items&&... items)
    -> std::ostream&
{
    {
        stream << first;
        char const _list[]{(stream << ' ' << items, char(0))...};
    }
    return stream;
}

auto operator<<(std::ostream& out, GridIndex index) -> std::ostream&
{
    return printSpaceDelimited(out, index.x, index.y, index.z);
}
template <class T>
auto operator<<(std::ostream& out, NumericRange<T> range) -> std::ostream&
{
    return out << "[" << range.min << ", " << range.max << "]";
}
auto operator<<(std::ostream& out, BoundingBox const& box) -> std::ostream&
{
    out << '(';
    return printSpaceDelimited(out, box.xBounds, box.yBounds, box.zBounds) << ')';
}

auto operator<<(std::ostream&                                     stream,
                std::unordered_map<GridIndex, BoundingBox> const& boxmap)
    -> std::ostream&
{
    for (auto& grid_box_pair : boxmap)
    {
        auto index = grid_box_pair.first;
        auto box   = grid_box_pair.second;

        stream << "Cell " << index << ": " << box << '\n';
    }
    return stream;
}
