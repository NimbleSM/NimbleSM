#pragma once
#include <utility>
#include "GridHash.h"
#include "NumericRange.h"
#include "defines.h"
#include "point.h"

struct BoundingBox
{
    using Range = NumericRange<double>;
    static void fitBounds(double val, double& lower, double& upper) noexcept
    {
        if (val < lower)
            lower = val;
        else if (val > upper)
            upper = val;
    }
    Range xBounds, yBounds, zBounds;
    constexpr BoundingBox(Point3<double> point) noexcept
      : xBounds(point.x) /* xBounds = [x, x] */
      , yBounds(point.y) /* yBounds = [y, y] */
      , zBounds(point.z) /* zBounds = [z, z] */
    {
    }
    constexpr BoundingBox(Range xBounds, Range yBounds, Range zBounds) noexcept
      : xBounds(xBounds) /* Invokes copy constructor; it's cheap */
      , yBounds(yBounds)
      , zBounds(zBounds)
    {
    }

    constexpr BoundingBox()                   = default;
    constexpr BoundingBox(BoundingBox const&) = default;
    constexpr BoundingBox(BoundingBox&&)      = default;

    auto operator=(BoundingBox const&) noexcept -> BoundingBox& = default;
    auto operator=(BoundingBox&&) noexcept -> BoundingBox& = default;

    auto centerAt(Point3<double> point) & noexcept -> BoundingBox&
    {
        xBounds.centerAt(point.x);
        yBounds.centerAt(point.y);
        zBounds.centerAt(point.z);
        return *this;
    }
    auto centerAt(Point3<double> point) && noexcept -> BoundingBox&&
    {
        return std::move(centerAt(point));
    }
    auto include(Point3<double> point) & noexcept -> BoundingBox&
    {
        xBounds.include(point.x);
        yBounds.include(point.y);
        zBounds.include(point.z);
        return *this;
    }
    auto include(Point3<double> point) && noexcept -> BoundingBox&&
    {
        return std::move(include(point));
    }
    auto include(Point3<double> point) const & noexcept -> BoundingBox
    {
        BoundingBox box(*this);
        box.include(point);
        return box;
    }
    auto include(BoundingBox const& box) & noexcept -> BoundingBox&
    {
        xBounds.include(box.xBounds);
        yBounds.include(box.yBounds);
        zBounds.include(box.zBounds);
        return *this;
    }
    auto include(BoundingBox&& box) && noexcept -> BoundingBox&&
    {
        return std::move(include(box));
    }
    template <class Key, class Map>
    constexpr auto includeInMap(Map&& m, Key&& key) const -> void
    {
        auto iter = m.find(key);
        if (iter == m.end())
        {
            m[key] = *this;
        }
        else
        {
            iter->second.include(*this);
        }
    }
    constexpr auto contains(Point3<double> point) const noexcept -> bool
    {
        return xBounds.contains(point.x) && yBounds.contains(point.y)
               && zBounds.contains(point.z);
    }
};
