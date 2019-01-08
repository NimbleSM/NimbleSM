#pragma once 
#include "GridHash.h"
#include "GridIndex.h"
#include <utility>

class GridCellBounds {
   private:
    //This one is private because scale should always be the inverse of cell_size
    //Calling the GridCellBounds(x, y, z, cell_size) constructor 
    //Invokes this one with scale=1.0/cell_size
    constexpr GridCellBounds(double x, double y, double z, double cell_size, double scale) noexcept
        : lowerX(GridHash::ifloor(x * scale) * cell_size)
        , lowerY(GridHash::ifloor(y * scale) * cell_size)
        , lowerZ(GridHash::ifloor(y * scale) * cell_size) 
        , cell_size(cell_size) {}
   public:
    double lowerX, lowerY, lowerZ;
    const double cell_size;

    constexpr GridCellBounds(GridIndex index, double cell_size) noexcept 
        : lowerX(index.xIndex * cell_size)
        , lowerY(index.yIndex * cell_size)
        , lowerZ(index.zIndex * cell_size)
        , cell_size(cell_size) {}
    constexpr GridCellBounds(double x, double y, double z, double cell_size) noexcept
        : GridCellBounds(x, y, z, cell_size, 1.0 / cell_size) {}
    constexpr auto contains(double x, double y, double z) const noexcept -> bool {
        return (lowerX <= x) && (x < lowerX + cell_size) 
            && (lowerY <= y) && (y < lowerY + cell_size)
            && (lowerZ <= z) && (z < lowerZ + cell_size); 
    }
    constexpr auto contains(double const* points) const noexcept -> bool {
        return contains(points[0], points[1], points[2]); 
    }

    auto shiftTo(GridIndex index) & noexcept -> GridCellBounds & {
        lowerX = index.xIndex * cell_size;
        lowerY = index.yIndex * cell_size;
        lowerZ = index.zIndex * cell_size;
        return *this;
    } 
    auto shiftTo(GridIndex index) && noexcept -> GridCellBounds && {
        return std::move(*this);
    }
    auto shiftTo(double x, double y, double z) & noexcept -> GridCellBounds& {
        double scale = 1.0 / cell_size; 
        lowerX = GridHash::ifloor(x * scale) * cell_size;
        lowerY = GridHash::ifloor(y * scale) * cell_size;
        lowerZ = GridHash::ifloor(z * scale) * cell_size; 
        return *this; 
    }
    auto shiftTo(double x, double y, double z) && noexcept -> GridCellBounds&& {
        return std::move(shiftTo(x, y, z)); 
    }

    auto moveTo(GridIndex index) noexcept -> GridCellBounds& {
        lowerX = index.xIndex * cell_size;
        lowerY = index.yIndex * cell_size;
        lowerZ = index.zIndex * cell_size;
        return *this;
    }
};
