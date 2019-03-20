#pragma once
#include "BarrierTree.h"
#include "BoundingBox.h"
#include "DataChannel.h"
#include "GridCell.h"
#include "GridHash.h"
#include "GridIndex.h"
#include "RequestQueue.h"
#include "WaitAnyResult.h"
#include "meta.h"

#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "text_conversion.h"



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
template <class PacketMap>
void SyncSendMessageSizes(PacketMap&&   packets,
                          DataChannel   channel,
                          RequestQueue& queue,
                          size_t*       sizeStorage)
{
    for (auto&& packet : packets)
    {
        auto  destination = packet.first;
        auto& message     = packet.second;

        *sizeStorage = message.size();
        queue.push(channel.Issend(sizeStorage, 1, destination));
        sizeStorage += 1;
    }
}

template <class PacketMap, class SourceIdentifiedCallback, class SendProcessedCallback>
auto IdentifySources(PacketMap&&              packets,
                     DataChannel              channel,
                     DataChannel              barrierChannel,
                     SourceIdentifiedCallback notifySourceIdentified,
                     SendProcessedCallback    notifySendProcessed) -> void
{
    using namespace std;

    auto queue                 = RequestQueue();
    auto messageCounts         = vector<size_t>(packets.size());
    auto active_outgoing_count = size_t(messageCounts.size());
    auto barrier               = BarrierTree(barrierChannel);
    SyncSendMessageSizes(packets, channel, queue, messageCounts.data());
    barrier.enqueueBarrier(queue);

    auto incoming_size       = size_t();
    auto active_recv_request = channel.Irecv(&incoming_size, 1, MPI_ANY_SOURCE);
    queue.push(active_recv_request);

    do
    {
        auto reply = queue.pop();
        if (reply.request == active_recv_request)
        {
            notifySourceIdentified(reply.status.MPI_SOURCE, incoming_size);

            active_recv_request = channel.Irecv(&incoming_size, 1, MPI_ANY_SOURCE);
            queue.push(active_recv_request);
        }
        else if (reply.status.MPI_TAG == channel.tag)
        {
            notifySendProcessed(reply.status.MPI_SOURCE);

            active_outgoing_count--;
            if (active_outgoing_count == 0)
            {
                barrier.markComplete();
            }
        }
        else if (reply.status.MPI_TAG == barrier.tag())
        {
            barrier.processStatus(reply.status);
        }
    } while (!barrier.test());
}

template <class PacketMap,
          class DataT = ElemType<decltype(std::begin(iof<PacketMap&>())->second)>>
auto ExchangeData(PacketMap&& packets, MPI_Comm comm, int tag1, int tag2, int tag3)
    -> std::vector<std::pair<int, std::vector<DataT>>>
{
    using namespace std;

    auto incomingData = vector<pair<int, vector<DataT>>>();
    incomingData.reserve(10);

    auto         messageChannel = DataChannel({comm, tag3});
    RequestQueue queue;

    // clang-format off
    IdentifySources(
        packets,
        DataChannel({comm, tag1}),
        DataChannel({comm, tag2}),
        [&](int source, size_t incoming_size) {
            incomingData.emplace_back();
            auto& back   = incomingData.back();
            back.first   = source;
            auto& buffer = back.second;
            buffer.resize(incoming_size);
            queue.push(messageChannel.Irecv(buffer.data(), incoming_size, source));
        },
        [&](int dest) {
            auto& message = packets[dest];
            queue.push(messageChannel.Isend(message, dest));
        }
    );
    // clang-format on
    queue.wait_all();

    return incomingData;
}

auto anyIntersect(std::vector<BoundingBox> const& boxes1,
                  std::vector<BoundingBox> const& boxes2) -> bool
{
    for (auto b1 : boxes1)
    {
        for (auto const& b2 : boxes2)
        {
            if (b1.intersects(b2))
                return true;
        }
    }
    return false;
}

using BBPacket     = std::pair<int, std::vector<BoundingBox>>;
using BBPacketList = std::vector<BBPacket>;
/* returns a list of the ranks that this rank intersects with */
auto notifyRanksOfIntersection(BBPacketList const&     rankBoxInfo,
                               DataChannel             dataChannel,
                               DataChannel             sizeChannel,
                               std::vector<int> const& sourceRanks)
    -> std::vector<int>
{
    using namespace std;
    RequestQueue                 sizeRecvQueue;
    std::vector<size_t>          incomingSizes(sourceRanks.size());
    std::unordered_map<int, int> inverseSourceRanks;
    vector<vector<int>>          incomingRanks(sourceRanks.size());

    inverseSourceRanks.reserve(sourceRanks.size());

    // Prepare to recieve incoming messages
    {
        size_t index = 0;
        for (int rank : sourceRanks)
        {
            inverseSourceRanks[rank] = (int)index;
            sizeRecvQueue.push(sizeChannel.Irecv(&incomingSizes[index], 1, rank));
            ++index;
        }
    }

    RequestQueue sendQueue;
    auto         sizes             = vector<size_t>(rankBoxInfo.size());
    auto         intersectingRanks = vector<vector<int>>(rankBoxInfo.size());

    {
        size_t index = 0;
        for (auto& destInfo : rankBoxInfo)
        {
            int   destRank    = destInfo.first;
            auto& listOfRanks = intersectingRanks[index];

            for (auto& otherInfo : rankBoxInfo)
            {
                if (&destInfo == &otherInfo)
                    continue;
                if (anyIntersect(destInfo.second, otherInfo.second))
                    listOfRanks.push_back(otherInfo.first);
            }

            sendQueue.push(sizeChannel.Isend(&sizes[index], 1, destRank));
            sendQueue.push(dataChannel.Isend(listOfRanks, destRank));
            index++;
        }
    }


    RequestQueue rankRecvQueue;
    size_t       total_incoming_ranks = 0;
    while (sizeRecvQueue.has())
    {
        WaitAnyResult result = sizeRecvQueue.pop();
        if (result.status.MPI_TAG == sizeChannel.tag)
        {
            int   source       = result.status.MPI_SOURCE;
            int   source_index = inverseSourceRanks.at(source);
            auto& vect         = incomingRanks[source_index];
            vect.resize(incomingSizes[source_index]);
            rankRecvQueue.push(dataChannel.Irecv(vect.data(), vect.size(), source));
            total_incoming_ranks += vect.size();
        }
    }
    rankRecvQueue.wait_all();


    std::unordered_set<int> ranks;
    ranks.reserve(total_incoming_ranks);
    for (auto& rankList : incomingRanks)
    {
        ranks.insert(rankList.begin(), rankList.end());
    }
    return std::vector<int>(ranks.begin(), ranks.end());
}
template <class View>
auto getExchangeMembers(View&&       kokkos_view,
                        double const cell_size,
                        DataChannel  channel) -> std::vector<int>
{
    using namespace std;
    HashFunction hash;
    int const    n_ranks = channel.commSize();

    // clang-format off
    auto packets = GatherBy(
        addPointsToBoundingBoxMap(
            kokkos_view, 
            cell_size, 
            BoundingBoxMap()
        ),
        [=](pair<const GridIndex, BoundingBox> const& message) {
            auto index = message.first;
            return hash(index.x, index.y, index.z) % n_ranks;
        }
    ); 
    std::vector<int> source_ranks; 
    for(auto& packet : packets) {
        source_ranks.push_back(packet.first); 
    }
    auto rankBoxInfo = std::vector<std::pair<int, std::vector<BoundingBox>>>(
        ExchangeData(
            packets,
            channel.comm,
            channel.tag, 
            channel.tag + 1,
            channel.tag + 2
        )
    );
    // clang-format on
    return notifyRanksOfIntersection(rankBoxInfo,
                                     channel,
                                     DataChannel{channel.comm, channel.tag + 1},
                                     source_ranks);
}
