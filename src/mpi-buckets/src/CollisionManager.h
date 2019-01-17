#pragma once
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include "BoundingBox.h"
#include "DataChannel.h"
#include "GridCell.h"
#include "GridHash.h"
#include "GridIndex.h"
#include "RequestQueue.h"
#include "WaitAnyResult.h"
#include "meta.h"

auto operator<<(std::ostream& out, GridIndex index) -> std::ostream&
{
    return out << index.x << " " << index.y << " " << index.z;
}
template <class T>
auto operator<<(std::ostream& out, NumericRange<T> range) -> std::ostream&
{
    return out << "[" << range.min << ", " << range.max << "]";
}
auto operator<<(std::ostream& out, BoundingBox const& box) -> std::ostream&
{
    return out << "(" << box.xBounds << " " << box.yBounds << " " << box.zBounds
               << ")";
}
void printIndexBoxMap(std::unordered_map<GridIndex, BoundingBox> const& boxmap)
{
    std::stringstream stream;
    for (auto& grid_box_pair : boxmap)
    {
        auto index = grid_box_pair.first;
        auto box   = grid_box_pair.second;

        stream << "Cell " << index << ": " << box << "\n";
    }
    std::cerr << stream.str();
}

using BoundingBoxMap = std::unordered_map<GridIndex, BoundingBox>;
template <class Range, class Func, template <class...> class map = std::unordered_map>
using GatherByMap_t
    = map<DecayedOutputType<Func, ElemType<Range>>, std::vector<ElemType<Range>>>;

template <class View, class Map = BoundingBoxMap>
auto addPointsToBoundingBoxMap(View&&       values,
                               double const cell_size,
                               Map&&        boundingBoxLookup = Map{})
    -> decltype(boundingBoxLookup)
{
    size_t const count = values.extent(0);
    if (count > 0)
    {
        double const scale = 1.0 / cell_size;

        Point3d point{values(0, 0), values(0, 1), values(0, 2)};

        auto current_index       = GridIndex(point, scale);
        auto current_cell_bounds = GridCellBounds(current_index, cell_size);
        auto current_box         = BoundingBox(point);

        for (size_t i = 1; i < count; i++)
        {
            point = {values(i, 0), values(i, 1), values(i, 2)};

            if (current_cell_bounds.contains(point, cell_size))
            {
                current_box.include(point);
            }
            else
            {
                // We include the current bounding box in the map,
                // then calculate a new bounding box, index, and cell
                current_box.includeInMap(boundingBoxLookup, current_index);
                current_box.centerAt(point);
                current_index = GridIndex(point, scale);
                current_cell_bounds.moveTo(current_index, cell_size);
            }
        }
        // Ensure that the most recent version of the current box is included in
        // the map
        current_box.includeInMap(boundingBoxLookup, current_index);
    }
    return std::forward<decltype(boundingBoxLookup)>(boundingBoxLookup);
}

template <class Range, class Func, class Destination = GatherByMap_t<Range, Func>>
auto GatherBy(Range&& range, Func func, Destination&& dest = Destination{})
    -> decltype(dest)
{
    for (auto&& elem : range)
        dest[func(elem)].push_back(elem);

    return std::forward<decltype(dest)>(dest);
}

template <class RankMap>
void SendBoundingBoxData(RankMap&& rank_data_map, DataChannel channel)
{
    RequestQueue requests;
    for (auto&& rank_data_pair : rank_data_map)
    {
        auto  destRank = rank_data_pair.first;
        auto& message  = rank_data_pair.second;

        requests.push(channel.Isend(message, destRank));
    }

    int const my_rank = channel.commRank();

    DataChannel on_completion_channel{channel.comm, channel.tag + 1};

    // We expect the rank above us to tell us when to stop waiting for recieves
    // from the main channel The rank above us is basically floor(my_rank /  2),
    // that is, my_rank >> 1
    if (my_rank != 0)
    {
        requests.push(on_completion_channel.Iawait(my_rank >> 1));
    }

    do {
        
        auto reply = requests.pop(); 
        if(reply.status.MPI_TAG == channel.tag) {

        }
    } while(requests.has());
}

template <class View>
void handleCollisions(View&& kokkos_view, double const cell_size, DataChannel channel)
{
    HashFunction hash;
    int const    n_ranks = channel.commSize();

    using rank_data_type = ElemType<View>;
    // clang-format off
    SendBoundingBoxData(
        GatherBy(
            addPointsToBoundingBoxMap(
                kokkos_view, 
                cell_size, 
                BoundingBoxMap()
            ),
            [=](rank_data_type const& message) {
                auto index = message.first;
                return hash(index.x, index.y, index.z) % n_ranks;
            }
        ),
        channel
    );
    // clang-format on
}
