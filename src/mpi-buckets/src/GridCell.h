#pragma once 
#include "GridHash.h"
#include <utility>

struct GridCell {
    double lowerX, lowerY, lowerZ, cell_size;

    auto isInside(double x, double y, double z) const noexcept -> bool {
        return (lowerX <= x) && (x < lowerX + cell_size) 
            && (lowerY <= y) && (y < lowerY + cell_size)
            && (lowerZ <= z) && (z < lowerZ + cell_size); 
    }
    auto isInside(double const* points) const noexcept -> bool {
        return isInside(points[0], points[1], points[2]); 
    }

    auto shiftTo(double x, double y, double z) & noexcept -> GridCell& {
        double scale = 1.0 / cell_size; 
        lowerX = GridHash::ifloor(x * scale) * cell_size;
        lowerY = GridHash::ifloor(y * scale) * cell_size;
        lowerZ = GridHash::ifloor(z * scale) * cell_size; 
        return *this; 
    }
    auto shiftTo(double x, double y, double z) && noexcept -> GridCell&& {
        double scale = 1.0 / cell_size; 
        lowerX = GridHash::ifloor(x * scale) * cell_size;
        lowerY = GridHash::ifloor(y * scale) * cell_size;
        lowerZ = GridHash::ifloor(z * scale) * cell_size; 
        return std::move(*this); 
    }
};
