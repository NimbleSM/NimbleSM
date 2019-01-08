#pragma once
#include "GridHash.h"
#include "defines.h"
#include "NumericRange.h"
#include <utility>

struct BoundingBox {
    using Range = NumericRange<double>;
    static void fitBounds(double val, double& lower, double& upper) noexcept {
        if(val < lower) lower = val;
        else if(val > upper) upper = val; 
    }
    Range xBounds, yBounds, zBounds; 
    constexpr BoundingBox(double x, double y, double z) noexcept
        : xBounds(x)
        , yBounds(y)
        , zBounds(z) {}
    constexpr BoundingBox(Range xBounds, Range yBounds, Range zBounds) noexcept
        : xBounds(xBounds)
        , yBounds(yBounds)
        , zBounds(zBounds) {}
    
    constexpr BoundingBox() = default;
    constexpr BoundingBox(BoundingBox const&) = default;
    constexpr BoundingBox(BoundingBox &&) = default;
    auto operator=(BoundingBox const&) noexcept -> BoundingBox& = default;
    auto operator=(BoundingBox &&) noexcept -> BoundingBox& = default;

    auto centerAt(double x, double y, double z) & noexcept -> BoundingBox& {
        xBounds.centerAt(x);
        yBounds.centerAt(y);
        zBounds.centerAt(z);
        return *this;
    }
    auto centerAt(double x, double y, double z) && noexcept -> BoundingBox&& {
        return std::move(centerAt(x, y, z)); 
    }
    auto include(double x, double y, double z) & noexcept -> BoundingBox& {
        xBounds.include(x); 
        yBounds.include(y); 
        zBounds.include(z);
        return *this; 
    }
    auto include(double x, double y, double z) && noexcept -> BoundingBox&& {
        return std::move(include(x, y, z));
    }
    auto include(double x, double y, double z) const & noexcept -> BoundingBox {
        BoundingBox box(*this);
        box.include(x, y, z);
        return box;
    }
    auto include(BoundingBox const& box) & noexcept -> BoundingBox& {
        xBounds.include(box.xBounds);
        yBounds.include(box.yBounds);
        zBounds.include(box.zBounds);
        return *this;
    }
    auto include(BoundingBox && box) && noexcept -> BoundingBox&& {
        return std::move(include(box)); 
    }
    template<class Key, class Map>
    auto includeInMap(Map&& m, Key&& key) const -> void {
        auto iter = m.find(key); 
        if(iter == m.end()) {
            m[key] = *this;
        } else {
            iter->second.include(*this);
        }
    }
    auto contains(double x, double y, double z) const noexcept -> bool {
        return xBounds.contains(x) && yBounds.contains(y) && zBounds.contains(z); 
    }
};
