#pragma once
#include "DataChannel.h"
#include "WaitAnyResult.h"
#include "RequestQueue.h"
#include "GridHash.h"
#include "GridCell.h"
#include "GridIndex.h"
#include "BoundingBox.h"
#include <cstddef>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <string>

auto operator<<(std::ostream& out, GridIndex index) -> std::ostream& {
    return out << index.xIndex << " " << index.yIndex << " " << index.zIndex;
}
template<class T>
auto operator<<(std::ostream& out, NumericRange<T> range) -> std::ostream& {
    return out << "[" << range.min << ", " << range.max << "]"; 
}
auto operator<<(std::ostream& out, BoundingBox const& box) -> std::ostream& {
    return out << "(" << box.xBounds << " " << box.yBounds << " " << box.zBounds << ")"; 
}
void printIndexBoxMap(std::unordered_map<GridIndex, BoundingBox> const & boxmap) {
    std::stringstream stream; 
    for(auto & grid_box_pair : boxmap) {
        auto index = grid_box_pair.first;
        auto box = grid_box_pair.second;

        stream << "Cell " << index << ": " << box << "\n"; 
    }
    std::cerr << stream.str();
}



void processCollision(double* ptr, size_t count, double const cell_size) {
    if(count < 3) return; 
    double const scale = 1.0 / cell_size;

    auto current_index = GridIndex(ptr[0], ptr[1], ptr[2], scale); 
    auto current_cell = GridCellBounds(current_index, cell_size); 

    using Map = std::unordered_map<GridIndex, BoundingBox>;

    auto boxes = Map{};
    auto *current_box = & boxes[current_index]; 
    *current_box = BoundingBox(ptr[0], ptr[1], ptr[2]); 


    for(size_t i = 0; i < count; i += 3) {
        double x = ptr[i]; 
        double y = ptr[i + 1];
        double z = ptr[i + 2];

        if(!current_cell.contains(x, y, z)) {
            current_index = GridIndex(scale, x, y, z);
            current_cell = current_index;
            auto iter = boxes.find(current_index);
            if(iter != boxes.end()) {
                current_box = &iter->second; 
            } else {
                current_box = &boxes[current_index]; 
            }
        }
        current_box->include(x, y, z); 
    }

    printIndexBoxMap(boxes); 
    //We have a complete map containing all bounding boxes
}