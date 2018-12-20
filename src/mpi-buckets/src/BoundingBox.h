#pragma once
#include "GridHash.h"
#include "defines.h"

struct BoundingBox {
    /**
     * @brief The minimum x, y, and z coordinates of the bounding box
     *
     */
    double x, y, z;
    /**
     * @brief The dimensions of the box
     *
     */
    double xSize, ySize, zSize;

    template <class F>
    CONSTEXPR_CPP14 void getHashPoints(GridHash const hash, F&& sendHashes) const noexcept {
        int64_t min_x = hash.index(x), 
                min_y = hash.index(y), 
                min_z = hash.index(z);
        int64_t max_x = hash.index(x + xSize), 
                max_y = hash.index(y + ySize),
                max_z = hash.index(z + zSize);
        for(int64_t x_index = min_x; x_index <= max_x; ++x_index) {
            for(int64_t y_index = min_y; y_index <= max_y; ++y_index) {
                for(int64_t z_index = min_z; z_index <= max_z; ++z_index) {
                    sendHashes(hash(x_index, y_index, z_index)); 
                }
            }
        }
    }
};
