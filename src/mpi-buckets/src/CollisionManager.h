#pragma once
#include "DataChannel.h"
#include "WaitAnyResult.h"
#include "RequestQueue.h"
#include "GridHash.h"
#include "GridCell.h"
#include "GridIndex.h"
#include "BoundingBox.h"
#include "meta.h"
#include <cstddef>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <sstream>
#include <string>

auto operator<<(std::ostream& out, GridIndex index) -> std::ostream& {
    return out << index.xIndex << " " << index.yIndex << " " << index.zIndex;
}
template < class T >
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

template <
    class View, 
    class Map = std::unordered_map < GridIndex, BoundingBox >
> 
auto addPointsToBoundingBoxMap(View&& pointSource, double const cell_size, Map&& boundingBoxLookup = Map{}) 
    -> decltype(boundingBoxLookup)
{
    size_t const count = pointSource.extent(0); 
    if(count > 0) 
    {
        double const scale = 1.0 / cell_size;

        double x = pointSource(0, 0);
        double y = pointSource(0, 1);
        double z = pointSource(0, 2);
        
        auto current_cell_bounds = GridCellBounds(current_index, cell_size); 
        auto current_box = BoundingBox(x, y, z); 
        auto current_index = GridIndex(x, y, z, scale); 

        for(size_t i = 1; i < count; i++) {
            x = pointSource(i, 0);
            y = pointSource(i, 1);
            z = pointSource(i, 2);

            if(current_cell_bounds.contains(x, y, z)) 
            {
                current_box.include(x, y, z); 
            }
            else 
            {
                //Include the box in the map
                current_box.includeInMap(boundingBoxLookup, current_index);
                //Recenter the bounding box
                current_box.centerAt(x, y, z);
                //Recalculate the current index 
                current_index = GridIndex(x, y, z, scale);
                //Move the cell bounds to the new index
                current_cell_bounds.moveTo(current_index);
            }
        }
        //Ensure that the most recent version of the current box is included in the map
        current_box.includeInMap(boundingBoxLookup, current_index);
    }
    return std::forward<decltype(boundingBoxLookup)>(boundingBoxLookup);
}
template<class View, class Map>
auto addPointsToBoundingBoxMapConst(
    View&& pointSource, 
    double const cell_size, 
    double const scale,
    size_t const index, 
    size_t const count, 
    GridCellBounds const current_cell_bounds,
    BoundingBox const current_box, 
    GridIndex const current_index,
    Map&& boundingBoxLookup = Map{}) 
    -> decltype(boundingBoxLookup)
{
    if(index == count) {
        current_cell_bounds.includeInMap(boundingBoxLookup, current_index);
        return std::forward<decltype(boundingBoxLookup)>(boundingBoxLookup);
    } else {
        double const x = pointSource(0, 0);
        double const y = pointSource(0, 1);
        double const z = pointSource(0, 2);
        if (current_cell_bounds.contains(x, y, z)) {
            return addPointsToBoundingBoxMapConst(
                pointSource, 
                cell_size,
                scale,
                index + 1, 
                count,
                current_cell_bounds,
                current_box.include(x, y, z),
                current_index,
                boundingBoxLookup
            );
        } else {
            current_cell_bounds.includeInMap(boundingBoxLookup, current_index);
            return addPointsToBoundingBoxMap(
                pointSource,
                cell_size,
                scale,
                index + 1,
                count, 
                GridCellBounds(GridIndex(x, y, z, scale), cell_size),
                BoundingBox(x, y, z),
                GridIndex(x, y, z, scale),
                boundingBoxLookup
            );
        }
    }
}


template < 
    class Range,
    class Func,
    class Destination = std::unordered_map <
        DecayedOutputType < 
            Func, 
            ElemType < Range >
        >, 
        std::vector < 
            ElemType < Range > 
        >
    >
>
auto GatherBy(Range&& range, Func func, Destination&& dest = Destination{}) 
    -> decltype(dest) 
{
    for(auto&& elem : range)
        dest[func(elem)].push_back(elem);
    
    return std::forward<decltype(dest)>(dest);
}