#pragma once
#include <utility>
#include "GridHash.h"
#include "GridIndex.h"
#include "point.h"

class GridCellBounds : Point3<double>
{
   public:
    struct FloorAndScale
    {
        double         cell_size;
        double         scale;
        constexpr auto operator()(double val) const noexcept -> double
        {
            return GridHash::ifloor(val * scale) * cell_size;
        }
    };

   private:
    using Base = Point3<double>;
    // This one is private because scale should always be the inverse of
    // cell_size Calling the GridCellBounds(x, y, z, cell_size) constructor
    // Invokes this one with scale=1.0/cell_size
    constexpr GridCellBounds(Point3d const point,
                             double const  cell_size,
                             double const  scale) noexcept
      : Base{GridHash::ifloor(point.x * scale) * cell_size,
             GridHash::ifloor(point.y * scale) * cell_size,
             GridHash::ifloor(point.z * scale) * cell_size}
    {
    }

   public:
    constexpr GridCellBounds(GridIndex index, double cell_size) noexcept
      : Base{index.x * cell_size, index.y * cell_size, index.z * cell_size}
    {
    }
    constexpr GridCellBounds(Base base, double cell_size) noexcept
      : GridCellBounds(base, cell_size, 1.0 / cell_size)
    {
    }
    constexpr auto contains(Base point, double cell_size) const noexcept -> bool
    {
        return ((x <= point.x) && (point.x < x + cell_size))
               && ((y <= point.y) && (point.y < y + cell_size))
               && ((z <= point.z) && (point.z < z + cell_size));
    }

    auto moveTo(GridIndex index, double cell_size) noexcept -> GridCellBounds&
    {
        x = index.x * cell_size;
        x = index.y * cell_size;
        x = index.z * cell_size;
        return *this;
    }
};
