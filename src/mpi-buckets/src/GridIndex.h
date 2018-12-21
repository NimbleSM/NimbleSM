#pragma once
#include "HashFunction.h"
#include "GridHash.h"

struct GridIndex {
    int64_t xIndex, yIndex, zIndex; 

    constexpr GridIndex(double x, double y, double z, double scale) noexcept
        : GridIndex(x * scale, y * scale, z * scale) {}
    constexpr GridIndex(double x, double y, double z) noexcept
        : GridIndex(GridHash::ifloor(x), GridHash::ifloor(y), GridHash::ifloor(z)) {}
    constexpr GridIndex(int64_t xIndex, int64_t yIndex, int64_t zIndex) noexcept
        : xIndex(xIndex)
        , yIndex(yIndex)
        , zIndex(zIndex) {}

    //constexpr GridIndex(GridIndex const&) = default;
    //constexpr GridIndex(GridIndex &&) = default; 
    auto operator==(GridIndex index) const noexcept -> bool {
        return xIndex == index.xIndex && yIndex == index.yIndex && zIndex == index.zIndex; 
    }
    auto operator!=(GridIndex index) const noexcept -> bool {
        return xIndex != index.xIndex || yIndex != index.yIndex || zIndex != index.zIndex; 
    }
};


template<>
struct std::hash<GridIndex> {
    using argument_type = GridIndex; 
    using result_type = uint64_t; 

    HashFunction applyHash;
    constexpr hash() = default;
    constexpr hash(uint64_t a, uint64_t b) noexcept : applyHash(a, b) {}
    constexpr hash(HashFunction const& hashFunc) noexcept : applyHash(hashFunc) {}

    constexpr hash(hash const&) = default;
    constexpr hash(hash &&) = default;

    constexpr auto operator()(GridIndex const& g) const noexcept -> uint64_t {
        return applyHash(g.xIndex, g.yIndex, g.zIndex); 
    }
};
