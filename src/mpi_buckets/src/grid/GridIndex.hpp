#pragma once
#include "GridHash.hpp"
#include "HashFunction.hpp"
#include "point.hpp"

struct GridIndex : public Point3<int64_t>
{
    using Base = Point3<int64_t>;

    constexpr GridIndex() : Base{} {}
    constexpr GridIndex(Point3<double> point, double scale) noexcept
      : Base({GridHash::ifloor(point.x * scale),
              GridHash::ifloor(point.y * scale),
              GridHash::ifloor(point.z * scale)})
    {
    }
    constexpr GridIndex(Point3<double> point) noexcept
      : Base({GridHash::ifloor(point.x),
              GridHash::ifloor(point.y),
              GridHash::ifloor(point.z)})
    {
    }
    GridIndex(GridIndex const&) = default;
    GridIndex(GridIndex&&) = default; 
    
    GridIndex& operator=(GridIndex const&) = default; 
    GridIndex& operator=(GridIndex&&) = default; 

    // constexpr GridIndex(GridIndex const&) = default;
    // constexpr GridIndex(GridIndex &&) = default;
    auto operator==(GridIndex index) const noexcept -> bool
    {
        return x == index.x && y == index.y && z == index.z;
    }
    auto operator!=(GridIndex index) const noexcept -> bool
    {
        return x != index.x || y != index.y || z != index.z;
    }
};

template <>
struct std::hash<GridIndex>
{
    using argument_type = GridIndex;
    using result_type   = uint64_t;

    HashFunction applyHash;
    constexpr hash() = default;
    constexpr hash(uint64_t a, uint64_t b) noexcept : applyHash(a, b) {}
    constexpr hash(HashFunction const& hashFunc) noexcept : applyHash(hashFunc)
    {
    }

    constexpr hash(hash const&) = default;
    constexpr hash(hash&&)      = default;

    constexpr auto operator()(GridIndex const& g) const noexcept -> uint64_t
    {
        return applyHash(g.x, g.y, g.z);
    }
};
